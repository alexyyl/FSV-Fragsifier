#!/usr/bin/env python
"""
##### FSV string functions
Functions for string manipulation

Alexander YY Liu | yliu575@aucklanduni.ac.nz
"""
import numpy as np
import pandas as pd

# Reverse complement function
complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    return ''.join([complement_dict.get(x,x) for x in reversed(seq)])

# Method to find all repeat units in a string
# returns a list of repeat stats, such as longest run number-
# and a nested list of rus found in form of tuples: (repeat unit, start pos, end pos, number of repeat cycles)

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

