# # 001_sequencing_data_processing: 
# 
# This script performes the quality control with FastQC for the Nextera NovaSeq data. It then merges the reads and concatenated the lanes.
# 
# This script assumes the following folder layout
# - Home folder
#     - Scripts folder with this script
#     - Data/Sequencing_data folder with the raw sequencing data
#     - Data/uni_amplicon_sequences with the indexes needed to demultiplex the data
#     - Data/Analysis_NovaSeq for the results 
#     
# This script assumes the following programs are installed (paths to the executables need to be specified below):
# - FastQC
# - Pandaseq

import os
import subprocess

import pandas as pd
print(pd.__name__, pd.__version__)

import numpy as np
print(np.__name__, np.__version__)

import matplotlib.pyplot as plt
import matplotlib
print(matplotlib.__name__, matplotlib.__version__)

import scipy.stats as stats
import scipy
print(scipy.__name__, scipy.__version__)

import re
print(re.__name__, re.__version__)

from collections import Counter

import sys
print(sys.version)

import seaborn as sns
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec

from Bio import SeqIO

## Define FastQC path 
fastqc_path = '/path/to/fastqc'

## Define some helper functions that will help us for the further analysis
def reverse_complement(dna):
    """ function that reverse complements DNA
    dna: input dna sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def get_dict_of_seq(fasta_file):
    """ function that converts a fasta file to a dictionnary of sequences
    fasta_file: the input fasta file
    """
    
    file_fasta_dict = {}
    # output dict of imported seqs
    
    with open(fasta_file, 'r') as fasta:    
        for line in fasta:
            # loops through the file

            if line.startswith('>') == True:
                seq_info = line.strip('>').strip('\n').split('\t')[0]
                file_fasta_dict[seq_info] = ''
                # checks if seq header

            else:
                file_fasta_dict[seq_info] += line.strip('\n')
                # If not, append nt to current seq
                
    return file_fasta_dict

### This directory must be the Novaseq directory
os.chdir("/path/Novaseq_data")

## Import the indexes to check the demultiplex

universal_seqs_fasta = './Data/uni_amplicon_sequences/universal_amplicon_novaseq.fa'
universal_seqs = get_dict_of_seq(universal_seqs_fasta)
print(universal_seqs)

indexes_fasta = './Data/uni_amplicon_sequences/indexes_novaseq.fa'
indexes = get_dict_of_seq(indexes_fasta)
print(indexes)

indexes_F = './Data/uni_amplicon_sequences/indexes_novaseq_F.fa'
indexes_F = get_dict_of_seq(indexes_F)
print(indexes_F)
indexes_R = './Data/uni_amplicon_sequences/indexes_novaseq_R.fa'
indexes_R = get_dict_of_seq(indexes_R)
print(indexes_R)

degen_seq = 'AAAAA'


# # Start the Quality Control with FASTQC

### Define the path for the raw data
path_to_R1_file = '/path/sequencing/data/*R1.fastq.gz'
path_to_R2_file = '/path/sequencing/data/*R2.fastq.gz'

# In[ ]:


## FastQC quality control for R1
subprocess.check_output(fastqc_path + ' ' + path_to_R1_file, shell=True)


## FastQC quality control for R2
subprocess.check_output(fastqc_path + ' ' +path_to_R2_file, shell=True)


# # Start merging with Pandaseq

### Initialize intermediate folders for the analysis
amplicon_sequences_path = "./Data/Analysis_Novaseq/amplicon_sequences"
os.makedirs(amplicon_sequences_path, exist_ok = True)

temp_path = './Data/Analysis_Novaseq/temp'
os.makedirs(temp_path, exist_ok = True)

### Define the path for the raw data
path_to_files="./Data/Sequencing_data/"

### Before I start the merge, I need to unzip the fastq files using gunzip to go from .fasta.gz to .fastq:
subprocess.check_output("gunzip ./Data/Sequencing_data/*.gz", shell=True)


### Create the directory where the merged files are going to be saved
os.makedirs("./Data/Analysis_Novaseq/merged_reads", exist_ok = True)
path_to_merged_files="./Data/Analysis_Novaseq/merged_reads"


### Prepare a dictionary to identify the reads in the directory for the merging

dict_path_to_files = {}

for filename in os.listdir(path_to_files):
    f = os.path.join(path_to_files, filename)
    f1 = f.rstrip("_fastqc")
    
    f3 = f1.rstrip("R1.")
    
    f4 = f3.rstrip("R2.")
    
    f5 = f4.rstrip("L001")

    sample_name=f5.lstrip("./Data/Sequencing_data/")
    
    sample_name2=sample_name.rstrip("_S")
        
    dict_path_to_files[sample_name2] = f5
    
        
dict_path_to_files


### Merge R1 and R2 

merged_file_list = []

# Loop through the files
for filepath in dict_path_to_files.values():

    filepath_for = filepath+'R1.fastq'
    filepath_rev = filepath+'R2.fastq'

    sample_read_count = 0

    with open(filepath_for, 'r') as source:

        for line in source:

            if line.startswith('@A'):

                sample_read_count += 1

    if sample_read_count >= 1:

        # Derive the name of the output file based on the input file
        filepath_out_prefix = './Data/Analysis_Novaseq/merged_reads/'
        
        filepath_out = os.path.join(filepath + '_merged.fasta')

        # Add to the list of merged files so we can keep working with it
        merged_file_list.append(filepath_out)
        
        ## Pandaseq arguments:
        # -f input file with forward reads
        # -r input file with reverse reads
        # -L maximum length for a sequence
        # -O maximum overlap for a sequence
        # -k kmers
        # -B allow unbarcoded sequences
        # -N eliminate all sequences with unknown nucleotides
        # -t threshold (minimum probability a sequence must have to assemble)
        # -T threads
        # -w output file in fasta.bz2 format
        panda_seq_call = 'pandaseq -f '+filepath_for+' -r '+filepath_rev+ ' -L 550 -O 400 -k 4 -B -N -t 0.5 -T 6 -w '+ filepath_out

        subprocess.check_output(panda_seq_call, shell=True)
        
        print('--------')


# # CONTINUE IN THE 2_AGGandALIGN NOTEBOOK
