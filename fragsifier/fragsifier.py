#!/usr/bin/env python
#
# [ Project Fragsifier ]
# Fragsifier main executable
#
# Main control flow script for data preprocessing and sequence extraction
#
# Alexander YY Liu | yliu575@aucklanduni.ac.nz


from __future__ import absolute_import
from tqdm import tqdm
from pathlib2 import Path
from multiprocessing import Pool
import pandas as pd
import argparse
import os
from termcolor import cprint

#if __name__ == '__main__':
if True:
    """
    Run this code block regardless of 
    """
    cprint(""" #######  #### ###  ##
         ###   ###  ##
 ####### ###   ###  ##
 ##      ###    #####
 ##   #####      ###

 == Project Fragsifier ==
 STR Read Fragment Classifier: Accurate forensic STR detection
     """, 'yellow')

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", help='input labelled sequences file containing untrimmed 351 bp ForenSeq reads')
    parser.add_argument("outdir", help='output directory the labelled reads')
    parser.add_argument("-fsl", "--flank_search_length", help="The number of bases adjacent to the repeat stretch to search for flanking sequences. An integer (default: %(default)s)", default=30)
    parser.add_argument("-ft", "--flank_threshold", help="The minimum alignment scores for a pair of flanking sequence. A float (default: %(default)s)", default=18)
    parser.add_argument("-io", "--inwards_offset", help="The number of bases into the repeat stretch to search for flanks. An integer (default: %(default)s)", default=4)
    parser.add_argument("-mSc", "--max_STR_complexity", help="The maximum number of repeat stretches allowed in a possible STR. An integer (default: %(default)s)", default=11)
    parser.add_argument("-mpfa", "--min_percentage_flank_aligned", help="The minimum percentage that a reference flanking sequence must be aligned. A float (default: %(default)s)", default=0.2)
    parser.add_argument("-mtr", "--min_tandem_repeat", help="How many times the motif needs to repeat to count as a repeat stretch. A float (default: %(default)s)", default=2)
    parser.add_argument("-nlrs", "--n_longest_repeat_stretches", help="How many of the longest stretches to consider as the major stretch of a possible STR. (default: %(default)s)", default=3)
    parser.add_argument("-nc", "--num_cores", help="Number of parallel processes to initailze. (default: %(default)s)", default=1)
    parser.add_argument("-pft", "--per_flank_threshold", help="The minimum alignment scores for each flanking sequence. A float (default: %(default)s)", default=2)
    parser.add_argument("-st", "--seq_threshold", help="The minimum sequence classification prediction probability. A float (default: %(default)s)", default=0.5)
    args = parser.parse_args()
    
    infile = args.inputfile
    outdir = Path(args.outdir)

    num_cores = int(args.num_cores)

    # Create output dir if it doesnt exist
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)

    # Import routines based on if main script is run as package vs script
    if __package__:
        from fragsifier.routines.string_functions import *
        from fragsifier.routines.preprocessing import *
        from fragsifier.routines.fragsify import *
    else:
        from routines.string_functions import *
        from routines.preprocessing import *
        from routines.fragsify import *

    # Added 10/06/18
    seq_ref_flank_length = 25

    # Load query sequences
    # Check if input file is FASTQ((k, *vals) for k, vals in library.items())
    # Use only the sequence lines
    with open(infile, 'r') as f:
        test_sequences = f.readlines()
        try:
            if test_sequences[0][0] == '@':
                print('input file %s is FASTQ format with %s sequences' % (infile, len(test_sequences)))
                test_sequences = [x.rstrip() for i, x in enumerate(test_sequences) if i % 4 == 1]
            else:
                print('input file %s is raw sequences format with %s sequences' % (infile, len(test_sequences)))
                test_sequences = [x.rstrip() for x in test_sequences]
        except IndexError:
            print('Input file is empty!')
            quit()

    flank_alignment_threshold_dict = make_flank_threshold_dict(min_percentage_flank_aligned=float(args.min_percentage_flank_aligned),
                                                               per_flank_threshold=int(args.per_flank_threshold))

    def fragsify_wrapper(input_sequence):
        return fragsify(input_sequence,
                        flank_search_length=int(args.flank_search_length),
                        inwards_offset=int(args.inwards_offset),
                        min_tandem_repeat=int(args.min_tandem_repeat),
                        n_longest_repeat_stretches=int(args.n_longest_repeat_stretches),
                        flank_threshold=int(args.flank_threshold),
                        seq_threshold=float(args.seq_threshold))

def main():
    '''Section code so that it runs in both script mode and package mode'''

    cprint('Performing sequence extraction with {} cores'.format(num_cores), 'yellow')

    predictions = []
    with Pool(num_cores) as p:
        predictions = list(tqdm(p.imap(fragsify_wrapper, test_sequences), total=len(test_sequences)))

    with open(outdir / (infile.split('/')[-1].split('.')[0] + '.extractions'), 'w') as f:
        f.writelines(['\t'.join(i) + '\n' for i in predictions])

        # Make ckseqs formatted output file

        # Count alleles and produce ckseqs standard input
    alleles = [':'.join(y.split(':')[:-1]) for x in predictions for y in x if y != '']

    extracted_sequences_dict = {}
    # Make sequence dictionary
    # sequence: F counts|R counts
    for result in list(Counter(alleles).items()):
        locus, orientation, seq = result[0].split(':')[:3]
        counts = result[1]
        if orientation == 'R':
            seq = reverse_complement(seq)
        try:
            extracted_sequences_dict[seq + ':' + locus][{'F': 0, 'R': 1}[orientation]] = counts
            # Calculate allele length
            # Column 4 is allele length
            extracted_sequences_dict[seq + ':' + locus][3] = calculate_allele_num(locus, len(seq))
        except KeyError:
            extracted_sequences_dict[seq + ':' + locus] = [0, 0, 0, 0, locus]
            extracted_sequences_dict[seq + ':' + locus][{'F': 0, 'R': 1}[orientation]] = counts
            # Calculate allele number
            extracted_sequences_dict[seq + ':' + locus][3] = calculate_allele_num(locus, len(seq))

    # Append to skseqs output file
    ckseqs_df = pd.DataFrame.from_dict(extracted_sequences_dict, orient='index').reset_index()
    ckseqs_df[2] = ckseqs_df[0] + ckseqs_df[1]
    ckseqs_df['index'] = [x.split(':')[0] for x in ckseqs_df['index']]
    ckseqs_df.columns = [1, 2, 3, 4, 5, 0]
    ckseqs_df = ckseqs_df.sort_values(by=0, ascending=True)
    ckseqs_df[[0, 1, 2, 3, 4, 5]].to_csv(outdir / (infile.split('/')[-1].split('.')[0] + '.ckseqs'), header=None,
                                         index=None)

    cprint('Processing complete!\n', 'green')

if __name__ == '__main__':
    main()
