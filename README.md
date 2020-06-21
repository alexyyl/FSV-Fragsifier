```
 #######  #### ###  ##
         ###   ###  ##
 ####### ###   ###  ##
 ##      ###    #####
 ##   #####      ###            
```

# Project Fragsifier: STR Sequence Fragment Classifier

Fragsifier is a forensic short tandem repeat (STR) sequence extraction tool for massively parallel sequencing data.

> Original article: [Forensic STR allele extraction using a machine learning paradigm](https://doi.org/10.1016/j.fsigen.2019.102194)

Fragsifier uses a machine learning approach to classify potential STRs in addition to the classic flanking sequence alignment technique.


# Algorithm overview
To find STRs on a read, the tool performs the following procedures: 

* Detect all repeat stretches present on the read
* Create a list of STR candidates containing all possible STR sequences (combinations of repeat stretches)
* Use k-mer/n-gram based sequence classifiers to predict the locus of each candidate
* Perform flanking sequence alignment on each STR candidate to locate the STR boundaries and generate an alignment score
* Candidate STRs with congruent sequence classification and flanking sequence alignment results are returned as a legitimate STR

Note: Fragsifier is designed to analyze repeat stretches so will skip SNP sequences.

# Installation
Fragsifier was developed in Python 3.6. It is recommened that the program be installed in a [Conda virtual environment](https://www.anaconda.com/products/individual). 

Fragsifier is now available on PyPI and can be installed via pip:

`pip install fragsifier`


Installing using the pip command will also install the required packages automatically. 


Alternatively, Fragsifier can be installed from the download source folder by running:

```
pip install .	# run in setup.py directory
```

Fragsifier can also be run directly from the source folder by running:
```
pip install requirements.txt	# install requirements first
python fragsifier.py		# type --help for a list of adjustable parameters
```

# Usage
Fragsifier works out of the box with the pretrained models as described in the paper. 
The tool readily extracts STR loci described in [STRait Razor v2s](https://doi.org/10.1016/j.fsigen.2017.03.013).


### Sequence model retraining
The sequence model can be retrained using custom data following the below procedures:

* Delete all the files in the `models` folder
* Add training sequence files to the `models` folder:
	* A `FSV_reference_sequences_examples.csv` file containing example sequence from each STR locus
	* A `FSV_training_negative_examples.csv` file negative (non-STR, noise) sequences
	* Examples of above files can be found in the models folder (prior to being deleted) 
* Make sure to update the `locus.config` file in the `routines` folder with new STR loci information
* Machine learning prediction models will be rebuilt automatically when the tool is next run. A 'building models' message will be displayed. 

Tip: The models can be retrained in the source folder when downloaded, and then installed to the system using `pip install .` afterwards.

# Release notes

## v1.0.1
* Reorganized code
* Code can now be run as a script in the local folder or be installed as a package
* Installation via pip now installs dependencies
* Fragsifier now available as a package on PyPI

---
 Alexander YY Liu | yliu575@aucklanduni.ac.nz




e 
