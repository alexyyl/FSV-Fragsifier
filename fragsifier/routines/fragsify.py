#!/usr/bin/env python
#
# [ Project Fragsifier ]
# Fragsifier sequence extraction algorithm
#
# STR sequence fragment classifier
#
# Alexander YY Liu | yliu575@aucklanduni.ac.nz

from Bio import pairwise2

import warnings
import numpy as np
align = pairwise2.align.localms

# Import routines based on if main script is run as package vs script
if '.'.join(__name__.split('.')[:-1]) == 'routines':
    from routines.string_functions import *
    from routines.preprocessing import *
else:
    from fragsifier.routines.string_functions import *
    from fragsifier.routines.preprocessing import *

# Hidden settings
allele_num = 0
known_strs_present = list()
show_table = 0

# Make percentage flankshow dictionary
def make_flank_threshold_dict(min_percentage_flank_aligned=0.2, per_flank_threshold=7):
    return dict([[x[0], (max(len(x[1][0]) * min_percentage_flank_aligned, per_flank_threshold),
                         max(len(x[1][1]) * min_percentage_flank_aligned, per_flank_threshold))] for x in
                 flank_dict.items()])

flank_alignment_threshold_dict = make_flank_threshold_dict()

# Flank alignment filter includes both per_flank_threshold and min_percentage_flank_aligned thresholds
def flank_alignment_filter(score, locus, orient):
    if score >= flank_alignment_threshold_dict[locus][orient]:
        return score
    else:
        return 0


def fragsify(sequence, flank_threshold=10, seq_threshold=0.5,
             flank_search_length=30, show_table=show_table, inwards_offset=4,
             known_strs_present=[], n_longest_repeat_stretches=3,
             max_repeat_stretches=11, min_tandem_repeat=2, allele_num=allele_num):

    # Calculate sequence length
    seq_length = len(sequence)

    # Find all repeat sites and indexes
    repeat_elements = [x for x in find_rus(sequence) if x[3] >= min_tandem_repeat]

    if len(repeat_elements) != 0 and len(repeat_elements[0]) == 0:
        return [""]

    # Calculate all combinations of repeat stretches
    # Repeat combination selection goes here (eg. proposed STR must be present in repeat_stretches)
    top_repeat_stretches = [y[0] for y in sorted([(i,x[-1]) for i,x in enumerate(repeat_elements)], key=lambda x: x[-1])[-n_longest_repeat_stretches:]]
    # max_STR_complexity limits the complexity of STRs/the max number of repeat stretches in the proposed sequence
    # Prioritize longest uninterrupted stretch (repeat_stretches)
    repeat_ranges = [[(left, right), (i,j) ] for i, left in enumerate([x[1] for x in repeat_elements])
                     for j,right in enumerate([x[2] for x in repeat_elements[i:i+max_repeat_stretches+1]], start=i)
                     if any([i <= y and j >= y for y in top_repeat_stretches])]

    try:
        base_ranges, repeat_ranges = list(map(list, zip(*repeat_ranges)))
    except ValueError:
        return [""]

    ### Sequence classification
    ## Prepare k-gram (n-gram/k-mer) combination for input into sequence classifier
    ## Instead of performing k-gram and one-hot encoding for each sequence by iteration, do once for input sequence and
    ## sample from the post-converted sequence

    # Perform k-gram conversion and encoding
    str_sequences_encoded = seq_vectorizer.transform(sequence[max(0,x-seq_ref_flank_length):min(seq_length, y-k+seq_ref_flank_length)] for x,y in base_ranges)

    ## Sequence prediction
    preds = STRSSC_predict_proba(sequence_model, str_sequences_encoded)
    preds = [list(zip(*preds))][0]

    if known_strs_present != []:
        preds[0] = tuple([x if x in known_strs_present else 'Neg' for x in preds[0]])

    ## Flanking sequence alignment

    # Compile alignment positions
    required_alignments = []
    # Side of repeat, locus, base position, repeat number
    for locus, ranges, pred_proba, stretch in zip(preds[0], base_ranges, preds[1], repeat_ranges):
        if pred_proba >= seq_threshold:
            required_alignments += [(i, locus, ranges[i], stretch[i]) for i in range(2)]

    required_alignments = set(x for x in required_alignments if 'Neg' not in x[1])

    # Perform alignment in those positions
    alignment_results_dict = dict()
    for side, locus, position, stretch in required_alignments:
        # Initialize dict entry
        if stretch not in alignment_results_dict:
            alignment_results_dict[stretch] = {0:{},1:{}}
            alignment_results_dict[stretch][side] = dict()
        try:
            if side:
                alignment_results_dict[stretch][side][locus] = align(sequence[position - inwards_offset:position + flank_search_length],
                      flank_dict[locus][1], 1, -2, -2, -2, penalize_end_gaps=False, one_alignment_only=True)[0]
            else:
                alignment_results_dict[stretch][side][locus] = align(sequence[max(position-flank_search_length,0):position+inwards_offset],
                                  flank_dict[locus][0], 1, -2, -2, -2, penalize_end_gaps=False, one_alignment_only=True)[0]
        except IndexError:
            alignment_results_dict[stretch][side][locus] = [[]]

    # Retrieve alignment scores
    left_align_score = []
    right_align_score = []

    for locus, stretch in zip(preds[0], repeat_ranges):
        try:
            left_align_score += [flank_alignment_filter(alignment_results_dict[stretch[0]][0][locus][2], locus, 0)]
        except (KeyError, IndexError):
            left_align_score += [0]
        try:
            right_align_score += [flank_alignment_filter(alignment_results_dict[stretch[1]][1][locus][2], locus, 1)]
        except (KeyError, IndexError):
            right_align_score += [0]

    alignment_score = [x+y for x,y in zip(left_align_score, right_align_score)]

    # Assemble decision frame
    decision_table = pd.DataFrame(preds+[left_align_score]+[right_align_score]+[alignment_score]).T
    decision_table.columns = ['sequence_pred_locus','sequence_pred_proba', 'left_align_score', 'right_align_score', 'alignment_score']

    if show_table == 1:
        decision_table['stretch_boundaries'] = base_ranges
        decision_table['stretches_contained'] = repeat_ranges
        decision_table = decision_table[decision_table.alignment_score != 0]
        boundaries = []
        for locus, str_index in zip(decision_table.sequence_pred_locus, decision_table.index):
            try:
                boundaries += [(((len(alignment_results_dict[repeat_ranges[str_index][0]][0][locus][1].rstrip('-')) -
                                  len(alignment_results_dict[repeat_ranges[str_index][0]][0][locus][0]) + inwards_offset +
                                  base_ranges[str_index][0])),
                                ((len(alignment_results_dict[repeat_ranges[str_index][1]][1][locus][0].lstrip('-')) -
                                  len(alignment_results_dict[repeat_ranges[str_index][1]][1][locus][1].lstrip(
                                      '-')) - inwards_offset + base_ranges[str_index][1])))]
            except KeyError:
                boundaries += [(0,0)]
        decision_table['boundaries'] = boundaries
        decision_table['allele'] = [':'+str(calculate_allele_num(decision_table['sequence_pred_locus'][x].split(':')[0],
                             decision_table['boundaries'][x][1] - decision_table['boundaries'][x][0])) for x in decision_table.index.tolist()]
        print(decision_table.sort_values(by=['alignment_score','sequence_pred_proba'], ascending=[False,False]))

    # Apply flank score filtering
    # Filter out entries that does not fulfill requirement
    # Flank classifier prediction scores must be above threshold
    decision_table = decision_table[(decision_table['alignment_score'] >= flank_threshold) &
                                        (decision_table['sequence_pred_proba'] >= seq_threshold)]

    # Sort values
    decision_table = decision_table.sort_values(by=['alignment_score','sequence_pred_proba'], ascending=[False,False])

    if len(decision_table) == 0:
        return ['']

    decision_table = decision_table.groupby('sequence_pred_locus').head(1)

    # Calculate range of selected STR proposals
    boundaries = []
    for locus, str_index in zip(decision_table.sequence_pred_locus, decision_table.index):
        boundaries += [(((len(alignment_results_dict[repeat_ranges[str_index][0]][0][locus][1].rstrip('-')) -
                          len(alignment_results_dict[repeat_ranges[str_index][0]][0][locus][0]) + inwards_offset + base_ranges[str_index][0])),
                       ((len(alignment_results_dict[repeat_ranges[str_index][1]][1][locus][0].lstrip('-')) -
                         len(alignment_results_dict[repeat_ranges[str_index][1]][1][locus][1].lstrip('-')) - inwards_offset + base_ranges[str_index][1] )))]
    decision_table['boundaries'] = boundaries

    # Produces locus:allele call:sequence:score
    out_seqs = []

    # Calculate allele number
    if allele_num == 1:
        decision_table['allele'] = [':'+str(calculate_allele_num(decision_table['sequence_pred_locus'][x].split(':')[0],
                             decision_table['boundaries'][x][1] - decision_table['boundaries'][x][0])) for x in decision_table.index.tolist()]
    else:
        decision_table['allele'] = ''

    # Format out sequence
    for x in decision_table.index.tolist():
        out_seqs += [':'.join([decision_table['sequence_pred_locus'][x]+decision_table['allele'][x],
                                sequence[decision_table['boundaries'][x][0]:decision_table['boundaries'][x][1]],
                                str(decision_table['alignment_score'][x])])]

    return out_seqs