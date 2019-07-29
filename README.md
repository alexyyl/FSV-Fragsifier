```
 #######  #### ###  ##
         ###   ###  ##
 ####### ###   ###  ##
 ##      ###    #####
 ##   #####      ###            
```
# Project Fragsifier: STR Sequence Fragment Classifier

We present a machine learning approach to short tandem repeat (STR) sequence detection and extraction called Fragsifier. Using massively parallel sequencing data, STRs are detected by first locating the longest repeat stretches on a read, then by locus prediction using k-mers with a machine learning sequence model, followed by reference flanking sequence alignment to determine precise STR boundaries. 

This tool uses the STR flanking sequence information of STRait Razor v2s, which is accessible at: https://www.unthsc.edu/graduate-school-of-biomedical-sciences/laboratory-faculty-and-staff/strait-razor/

# Algorithm overview
The tool performs the following steps in sequence: 
* Detect all repeat stretches in a read
* Create a list all possible STR sequences
* Use sequence classifiers to predict locus of possible STR
* Perform flanking sequence alignment to locate STR boundaries and generate alignment score
* Assign highest-scoring sequence as a legitimate STR

# Usage
Install requirements using the requirements.txt
Run Fragsifier.py and type --help for a list of adjustable parameters
Anaconda python is recommended.

## Using custom data to retrain
If you wish to retrain the sequence prediction model:
* Delete all files in the 'train' folder
* Add STR training sequences to the FSV_reference_sequences_examples.csv file 
* Add negative (non-STR) training sequences to the FSV_training_negative_examples.csv file 
* Make sure to update locus.config with new STR loci information
* All prediction models will be created upon next run

# Known issues
Tested with Anaconda with Python 3.6
Fragsifier currently only perform extraction of STR sequence.

To be updated...

---
 Alexander YY Liu | yliu575@aucklanduni.ac.nz





