#!/usr/bin/env python
"""
##### FSV sequence classification preprocessing

* Preprocess training data
* Training sequence classification model

Upon first fun: train the sequence classifier

# 22/05/19 Update
# 19/06/18 Update

Alexander YY Liu | yliu575@aucklanduni.ac.nz
"""

import pickle
import os
from collections import Counter
import random
from FSV_string_functions import *

# Set k value
k = 7 # Used 7-mers as sequence features
train_folder = 'train/'

seq_ref_flank_length = 25
num_sequences_to_augment_per_locus = 5

pd.set_option('display.expand_frame_repr', False)

# sklearn's one hot does not work with sequences with uneven lengths
# Make our own label ohe encoder
def labels_onehot(labels):
    encoded = np.zeros(len(feature_le.classes_))
    for i in labels:
        encoded[i] += 1
    return encoded

# 14/06/18 Added

def STRSSC_predict_proba(model, X_query):

    predictions = model.predict_proba(X_query)
    index_of_best_predictions = [np.argmax(x) for x in predictions]

    return [[class_le.inverse_transform(index), pred[index]] for pred, index in zip(predictions, index_of_best_predictions)]


## Load in preinitiated configurations if exists
## Or else do preprocess and build new models

## Load allele definitions
print('>','Loading allele definitions')

if os.path.isfile(train_folder + 'flank_dict.pkl') and os.path.isfile(train_folder+'allele_len_conversion_dict') and os.path.isfile(train_folder+'allele_num_calculation_dict'):
    flank_dict = pickle.load(open(train_folder+'flank_dict.pkl', 'rb'))
    allele_num_calculation_dict = pickle.load(open(train_folder + 'allele_num_calculation_dict.pkl', 'rb'))
    allele_len_conversion_dict = pickle.load(open(train_folder + 'allele_len_conversion_dict.pkl', 'rb'))
else:
    # Load in locus flanking sequence database
    locus_ref = pd.read_csv('locus.config', sep='\t', header=None)

    # Remove Amel and SNPs if they exist
    locus_ref = locus_ref.iloc[
        [i for i, x in enumerate(locus_ref[0].tolist()) if x[:2] not in ['rs', 'N2', 'mh', 'Am'] if
         x[-2:] != 'FS']].reset_index(drop=True)
    locus_ref[0] = [x.split('.')[0] for x in locus_ref[0]]

    # Build flanking sequence dictionary
    all_loci = [x + y for x in locus_ref[0].values for y in [':F', ':R']]
    all_flanks = [list(locus_ref.loc[x, y[0]:y[1]]) for x in range(len(locus_ref)) for y in [(2, 3), (4, 5)]]
    flank_dict = dict(zip(all_loci, all_flanks))
    pickle.dump(flank_dict, open(train_folder+'flank_dict.pkl', 'wb'))

    # Build locus dict contain dict of length to allele num conversions
    allele_len_conversion_dict = {}
    for i in range(len(locus_ref)):
        locus = locus_ref[0][i]
        list_allelenum, list_length = zip(*[x.split(':')[1].split('=') for x in locus_ref[7][i].split(';')[:-1]])
        conversion_dict = dict(zip([int(x) for x in list_length], list_allelenum))
        allele_len_conversion_dict[locus] = conversion_dict

    # 21/06/18
    # Build allele length conversion formulae to calculate allele number in cases where allele not recorded
    allele_num_calculation_dict = {}
    for locus in allele_len_conversion_dict.keys():
        # Figure out number of bases in repeat
        complete_alleles = [[x[0],int(x[1])] for x in allele_len_conversion_dict[locus].items() if '.' not in x[1]]
        # Find first consecutive pair
        for i in range(len(complete_alleles)):
            allele = complete_alleles[i]
            allele_next = complete_alleles[i+1]
            if allele[1] + 1 == allele_next[1]:
                motif_size = allele_next[0] - allele[0]
                offset = allele[0] - motif_size*allele[1]
                allele_num_calculation_dict[locus] = (motif_size, offset)
                break

    pickle.dump(allele_len_conversion_dict, open(train_folder + 'allele_len_conversion_dict.pkl', 'wb'))
    pickle.dump(allele_num_calculation_dict, open(train_folder + 'allele_num_calculation_dict.pkl', 'wb'))

print('Completed!')

## Load reference sequences and models
print('>','Loading sequence models')
if os.path.isfile(train_folder + 'sequence_model.pkl') and os.path.isfile(train_folder + 'key_loci_dict.pkl'):
    sequence_model = pickle.load(open(train_folder + 'sequence_model.pkl', 'rb'))
    key_loci_dict = pickle.load(open(train_folder + 'key_loci_dict.pkl', 'rb'))
    seq_vectorizer = pickle.load(open(train_folder + 'seq_vectorizer.pkl', 'rb'))
    class_le = pickle.load(open(train_folder + 'class_le.pkl', 'rb'))
else:
    print('>', 'Creating new model')
    print('>', 'Loading training data')

    # Load reference sequences and build new sequence models
    train_files = [x for x in os.listdir(train_folder) if 'csv' in x]
    sequence_reference = pd.DataFrame()

    # Read all training files except the negative examples
    for file in train_files:
        if file != 'FSV_training_negative_examples.csv':
            sequence_reference = sequence_reference.append(
                pd.read_csv(train_folder + file, header=None).sort_values(by=0))

    # Remove negatives sequences from training examples and use trailing noise instead
    sequence_reference.columns = ['locus', 'orientation', 'left_bound', 'right_bound', 'read']
    # Check that sequence must be present
    sequence_reference = sequence_reference[sequence_reference.read != '']
    # Remove Amel and SNPs
    total_loci = [x for x in set(sequence_reference.locus) if x[:2] not in ['rs', 'N2', 'mh', 'Am']]
    sequence_reference = sequence_reference[sequence_reference['locus'].isin(total_loci)]

    ## Calculate offset from last repeat
    print('>', 'Calibrating optimum flanking sequence offsets')

    for locus in total_loci:
        calibration_df = sequence_reference[sequence_reference.locus == locus].reset_index()
        calibration_df = calibration_df[calibration_df.read != '']
        calibration_df['read'] = [s[l:r] for s, l, r in
                                  zip(calibration_df.read, calibration_df.left_bound, calibration_df.right_bound)]
        # [0, 48, 100, 'TCTGTCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA', 882, 882]
        temp_read_info = calibration_df.groupby(
            ['orientation', 'left_bound', 'right_bound', 'read']).count().sort_values(by='locus',
                                                                                      ascending=False).reset_index().iloc[
                         0, :].values.tolist()
        # Record: left offset, right offset, number of reads, orientation
        orient = {0: ':F', 1: ':R'}[temp_read_info[0]]

    # Trim reads to make training sequences
    sequence_reference['orientation'] = [{0: ':F', 1: ':R'}[x] for x in sequence_reference['orientation']]
    trimmed_sequences = []
    for loc, o, l, r, s in sequence_reference.values.tolist():
        trimmed_sequences += [s[max(0, l - seq_ref_flank_length):r + seq_ref_flank_length]]
    sequence_reference['sequence'] = trimmed_sequences
    sequence_reference = sequence_reference[['locus', 'orientation', 'sequence', ]]
    sequence_reference = sequence_reference[sequence_reference.sequence != '']

    # Sample n sequences from each classsort=False
    reads_to_sample = 300
    sampling_df = pd.DataFrame()
    for locus in set(sequence_reference.locus):
        sampling_df = sampling_df.append(sequence_reference[sequence_reference.locus == locus].sample(reads_to_sample, replace=True), sort=False)
    sequence_reference = sampling_df

    ## Perform sequence augmentation
    print('>', 'Performing sequence augmentation')
    augmented_sequences = []
    augment_steps = int(seq_ref_flank_length/2)
    augment_step_size = 2   # Remove 2 bases each step

    for locus in total_loci:

        temp_augmented_sequences = []
        try:
            temp_sequences = sequence_reference[sequence_reference.locus == locus].sample(num_sequences_to_augment_per_locus).values.tolist()
        except ValueError:
            print(locus)

        for info in temp_sequences:
            # Remove m bp from left (5') flank
            # Remove m bp from right (3') flank
            # Remove m bp from both flanks
            orient, seq = info[1:]
            for i in range(1,augment_steps):
                temp_augmented_sequences += [[locus, orient, seq[i * augment_step_size:]]]
                temp_augmented_sequences += [[locus, orient, seq[:-i * augment_step_size]]]
                temp_augmented_sequences += [[locus, orient, seq[i * augment_step_size:-i * augment_step_size]]]

        augmented_sequences += temp_augmented_sequences

    augmented_sequences = pd.DataFrame(augmented_sequences)
    augmented_sequences = augmented_sequences[augmented_sequences[2] != '']
    augmented_sequences.columns = ['locus', 'orientation', 'sequence',]

    # Append to dataset
    sequence_reference = sequence_reference.append(augmented_sequences, sort=False)

    # Load in trailing noise, format and append to training data
    print('>', 'Adding noise sequences')
    negatives_df = pd.read_csv(train_folder + 'FSV_training_negative_examples.csv', header=None).sort_values(by=0)
    negatives_df.columns = ['sequence']
    negatives_df['locus'] = 'Neg'
    negatives_df['orientation'] = ':F'

    # Trim negative noise
    trimmed_noise = []
    for seq in negatives_df.sequence:
        random_left = random.randint(0,int(len(seq)/2))
        random_right = random.randint(random_left,len(seq))
        trimmed_noise += [seq[random_left:random_right]]
    negatives_df['sequence'] = trimmed_noise

    # Add 3000 noise sequences to training sequences
    sequence_reference = sequence_reference.append(negatives_df.sample(3000), sort=False)

    # Data curation
    training_sequences = sequence_reference['sequence'].values.tolist()
    locus_label = sequence_reference['locus'].values.tolist()
    orientation_flags = sequence_reference['orientation'].values.tolist()

    # Enumerate reverse complement sequences
    print('>', 'Enumerating reverse orientation sequences')
    training_sequences += [reverse_complement(x) for x in training_sequences]
    orientation_flags += [{':R':':F',':F':':R'}[x] for x in orientation_flags]

    # Attach locus label orientations
    locus_label = [x + y for x, y in zip(locus_label+locus_label, orientation_flags)]

    sequence_reference = pd.DataFrame({'locus': locus_label, 'sequence': training_sequences})

    # make key to loci dictionary
    # find and record the 3 largest repeat units for each locus

    locus_key_dict = {}
    key_loci_dict = {}

    ### Model building
    # Use Countvectorizer to make k-mers/onehot encode
    print('>', 'Vectorizing sequences')

    from sklearn.feature_extraction.text import CountVectorizer

    # Character-level k-gram model method
    seq_vectorizer = CountVectorizer(analyzer="char", token_pattern=".", tokenizer=None, preprocessor=None, stop_words=None,
                                 ngram_range=(k, k))

    # Turn all sequences into kmers and encode
    X = seq_vectorizer.fit_transform(sequence_reference.sequence)

    # Label encode k-mers
    # Note that kgram used to train the encoder is k-mers extracted from reference database

    from sklearn.preprocessing import LabelEncoder

    # Encode class labels
    class_le = LabelEncoder()
    classes_encoded = class_le.fit_transform(sequence_reference.locus)
    y = classes_encoded

    # Train models
    # Removed decision tree and KNN
    print('>', 'Building model')

    #from sklearn.linear_model import LogisticRegression
    #sequence_model = LogisticRegression()

    #from sklearn.tree import DecisionTreeClassifier
    #sequence_model = DecisionTreeClassifier(random_state=0)

    from sklearn.ensemble import RandomForestClassifier
    sequence_model = RandomForestClassifier(n_estimators=30, random_state=0)

    sequence_model.fit(X, y)

    pickle.dump(sequence_model, open(train_folder + 'sequence_model.pkl', 'wb'))
    pickle.dump(key_loci_dict, open(train_folder + 'key_loci_dict.pkl', 'wb'))
    pickle.dump(seq_vectorizer, open(train_folder + 'seq_vectorizer.pkl', 'wb'))
    pickle.dump(class_le, open(train_folder + 'class_le.pkl', 'wb'))

print('Completed!')

# Classifier function
def STRSSC_predict(model, X_query):
    return class_le.inverse_transform(model.predict(X_query))


def calculate_allele_num(locus, allele_len):
    try:
        allele_num = allele_len_conversion_dict[locus][allele_len]
    except KeyError:
        if locus[:2] in ['rs', 'N2', 'mh', 'Am']:
            return 0
        else:
            conversion_rules = allele_num_calculation_dict[locus]
            if (allele_len - conversion_rules[1]) % conversion_rules[0] == 0:
                allele_num = str((allele_len - conversion_rules[1]) // conversion_rules[0])
            else:
                allele_num = '.'.join([str((allele_len - conversion_rules[1]) // conversion_rules[0]),
                                       str((allele_len - conversion_rules[1]) % conversion_rules[0])])
    return allele_num



