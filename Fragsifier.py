#!/usr/bin/env python
"""
##### Fragsifier exectuable

Base script to load input data file and process with Fragsifier algorithm

# Optimized parameters
# Fixed allele join problem

Outputs results in the format:
[0:locus, 1:sequence, 2:forward_counts, 3:reverse_counts, 4:base length, 5:allele call]

Alexander YY Liu | yliu575@aucklanduni.ac.nz
"""

from tqdm import tqdm
from multiprocessing import Pool
import pandas as pd
from FSV_string_functions import *
from FSV_preprocessing import *
import argparse
import os

if __name__ == '__main__':

    print(""" #######  #### ###  ##
         ###   ###  ##
 ####### ###   ###  ##
 ##      ###    #####
 ##   #####      ###

 == Project Fragsifier==
 STR Read Fragment Classifier: Accurate forensic STR detection
     """)

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("seqs", help='input labelled sequences file containing untrimmed 351 bp ForenSeq reads')
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
    infile = args.seqs
    outdir = args.outdir

    num_cores = int(args.num_cores)

    parameters = "flank_search_length = {}\ninwards_offset = {}\nmin_tandem_repeat = {}\nmax_STR_complexity = {}\nn_longest_repeat_stretches = {}\nmin_percentage_flank_aligned = {}\nper_flank_threshold = {}\nnum_cores = {}\nflank_threshold = {}\nseq_threshold = {}".format(
        int(args.flank_search_length), int(args.inwards_offset), int(args.min_tandem_repeat),
        int(args.max_STR_complexity), int(args.n_longest_repeat_stretches), float(args.min_percentage_flank_aligned),
        int(args.per_flank_threshold), int(args.num_cores), int(args.flank_threshold), float(args.seq_threshold))


    with open("FSV_fragsify.py") as f:
        lines = f.readlines()

    parameter_block_bounds = [i for i, x in enumerate(lines) if x == '#@\n']

    with open("TEMP_fragsify_runtime_params.py", "w") as f:
        f.write("".join(lines[:parameter_block_bounds[0]] + [parameters] + lines[parameter_block_bounds[1]+1:]))

    # Added 10/06/18
    seq_ref_flank_length = 25

    # Load query sequences
    # Check if input file is FASTQ
    # Use only the sequence lines
    print('processing file:', infile)
    with open(infile, 'r') as f:
        test_sequences = f.readlines()
        if test_sequences[0][0] == '@':
            print('input file is FASTQ')
            test_sequences = [x.rstrip() for i, x in enumerate(test_sequences) if i % 4 == 1]
        else:
            print('input file contains sequences')
            test_sequences = [x.rstrip() for x in test_sequences]
    
    print('found {} sequences'.format(len(test_sequences)))
    print('Begin processing with {} cores'.format(num_cores))

    # Use multiprocessing
    predictions = []

    with Pool(num_cores) as p:
        from TEMP_fragsify_runtime_params import *
        predictions = list(tqdm(p.imap(fragsify, test_sequences), total=len(test_sequences)))

    #print(outdir + infile.split('/')[-1].split('.')[0] + '.extractions')

    with open(outdir + infile.split('/')[-1].split('.')[0] + '.extractions', 'w') as f:
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
    ckseqs_df[2] = ckseqs_df[0]+ckseqs_df[1]
    ckseqs_df['index'] = [x.split(':')[0] for x in ckseqs_df['index']]
    ckseqs_df.columns = [1, 2, 3, 4, 5, 0]
    ckseqs_df = ckseqs_df.sort_values(by=0, ascending=True)
    ckseqs_df[[0, 1, 2, 3, 4, 5]].to_csv(outdir + infile.split('/')[-1].split('.')[0] + '.ckseqs', header=None, index=None)

    print('Processing complete!\n')

    os.remove('TEMP_fragsify_runtime_params.py')
