#!/usr/bin/env python
"""
##### FSV string functions
Core functions for string manipulation

# 19/06/18 Update

Alexander YY Liu | yliu575@aucklanduni.ac.nz
"""
import numpy as np
import pandas as pd

# Reverse complement from Stackoverflow
alt_map = {'ins': '0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    for k, v in alt_map.items():
        seq = seq.replace(k, v)
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    for k, v in alt_map.items():
        bases = bases.replace(v, k)
    return bases

# makes an alphabetical key
def alphabetical_transform(strin):
    key_template = strin + strin
    suffix_array = []
    for i in range(0, len(strin)):
        suffix_array += [key_template[i:]]
    suffix_array.sort()
    return suffix_array[0][0:len(strin)]

# Method to find all repeat units in a string
# returns a list of repeat stats, such as longest run number-
# and a nested list of rus found in form of tuples: (repeat unit, start pos, end pos, number of repeat cycles)
# changed deque processing to pure list processing
# waaaay faster (13 secs > 0.06 sec)

def find_rus(seq, k_range=range(3, 7)):
    found_repeats = []
    in_run_ru = 0
    run_counter = 0
    run_start_index = 0
    for i in range(len(seq)):
        if in_run_ru > 0:
            if seq[i:i + in_run_ru] == seq[i + in_run_ru:i + in_run_ru * 2]:
                # If already in run
                run_counter += 1
            else:
                found_repeats += [(seq[run_start_index:run_start_index + in_run_ru], run_start_index, i+(in_run_ru*2)-1, ((i+(in_run_ru*2)-1) - run_start_index) / in_run_ru)]
                in_run_ru = 0
                run_counter = 0
                run_start_index = 0
        else:
            # Search for new run
            for k in k_range:
                if seq[i:i + k] == seq[i+k:i+(k*2)]:
                    # A new run
                    in_run_ru = k
                    run_counter += 1
                    run_start_index = i
                    break

    # Terminate all runs at end of seq
    if in_run_ru == 1:
        # Run terminates
        found_repeats += [(seq[i:i + in_run_ru], run_start_index, len(seq), (i + run_start_index) / in_run_ru)]

    return found_repeats

