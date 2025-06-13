
# coding: utf-8

# # 002_sequencing_data_processing: 
# 
# This script aggregates the reads into the diferent unique variants counting how many reads cover each variant. After aggregation the script alignes the reads to the reference sequence.
# 
# This script assumes the following folder layout
# - Home folder
#     - Scripts folder with this script
#     - Data/uni_amplicon_sequences with the reference amplicon fasta file
#     - Data/Analysis_NovaSeq/merged_reads folder with the merged and concatenated data files with the reads obtained in the previous script
#     - Data/Analysis_NovaSeq folder for the results 
#     
# This script assumes the following programs are installed (paths to the executables need to be specified below):
# - Cutadapt
# - vsearch
# - Needle


### Import modules and check versions
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


os.chdir("/path/Novaseq_data")


## Make a dictionary to identify the different merged files

path_to_merged_files = "./Data/Analysis_Novaseq/merged_reads"

dict_path_to_merged_files = {}

for filename in os.listdir(path_to_merged_files):
    
    f = os.path.join(path_to_merged_files, filename)
    
    f1 = f.rstrip("_merged.fasta")
       
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/merged_reads")
           
    dict_path_to_merged_files[sample_name] = f1
    
dict_path_to_merged_files


# # Trimming Using cutadapt to keep only the sequence on which we performed the DMS from position 10 to position 78. For this I use the primers as adaptors AGTAGCAATGAAGTCAGT...GTTTAAACGGTCTCCAGC

### Create a directoty to save the trimmed files

outfile_trimmed_merged_prefix = './Data/Analysis_Novaseq/trimmed_merged_reads/'

os.makedirs(outfile_trimmed_merged_prefix, exist_ok = True)

## I will use cutadapt to trim the sequences using anchored adapters to remove all the sequence that are before the coding sequence. 

### Trim the first batch with 
trimmed_merged_list = []

## Loop through the merged files
for file in dict_path_to_merged_files.values():
    
    filename = file + '_merged.fasta'

    filename_out = file + '_trimmed_merged.fa'
    
    test2=filename_out.lstrip("./Data/Analysis_Novaseq/merged_reads/")
    
    outfile_trimmed_merged = os.path.join(outfile_trimmed_merged_prefix, test2)
    
    trimmed_merged_list.append(outfile_trimmed_merged)
    
    cutadapt_trim_call = 'cutadapt -a AGTAGCAATGAAGTCAGT...GTTTAAACGGTCTCCAGC -j 4 -o ' + outfile_trimmed_merged + " " + filename
            
    subprocess.check_output(cutadapt_trim_call, shell=True)
    
    print(test2)
    print('--------')
    


# # Aggregate the reads with Vsearch

aggregate_list = []

aggregate_prefix = './Data/Analysis_Novaseq/aggregate_reads/'
os.makedirs(aggregate_prefix, exist_ok = True)

## Loop through the trimmed merged files and aggregate
for trimmed_merged_file in trimmed_merged_list:
    
    
    test = trimmed_merged_file.rstrip("trimmed_merged.fa")
        
    filename_out = test + '_aggregate.fa'
        
    test2=filename_out.lstrip("./Data/Analysis_Novaseq/trimmed_merged_reads/")
        
    outfile_aggregate = os.path.join(aggregate_prefix, test2)
    
    aggregate_list.append(outfile_aggregate)
        
    ## vsearch parameters:
    # --derep_fulllength dereplicate sequences in the given fasta file
    # --relabel add a given prefix
    # --output output file
    # -- sizeout include abundance information

    vsearch_aggregate_call = 'vsearch --derep_fulllength ' + trimmed_merged_file + ' --relabel seq --sizeout --output '
    vsearch_aggregate_call += outfile_aggregate

    subprocess.check_output(vsearch_aggregate_call, shell=True)
    
    print(test2)
    print('-----')


# # Define and aling to the reference amplicon with Needle

os.mkdir('./Data/Analysis_Novaseq/amplicon_align')

path_to_amplicons = './Data/uni_amplicon_sequences/Novaseq_R67.fasta'

amplicon_info_dict = get_dict_of_seq(path_to_amplicons)

amplicon_length_dict = {}

amplicon_dict = {}

for amplicon in amplicon_info_dict:
    
    amplicon_name = amplicon.split('|')[0]

    amplicon_dict[amplicon_name] = amplicon_info_dict[amplicon]
    
    amplicon_fasta_path = './Data/Analysis_Novaseq/amplicon_sequences/'+amplicon_name+'.fasta'
    
    with open(amplicon_fasta_path, 'w') as dest:
        
        seq_ID = '>'+amplicon_name+'\n'
        seq = amplicon_dict[amplicon_name]+'\n'
        
        dest.write(seq_ID)
        dest.write(seq)
    

print (amplicon_info_dict.keys())
print (amplicon_dict)

## Define a function to do the alignments with needle
def needle_align_on_ref(ref_orf, filepath):
    
    """ function that creates the command and executes it in the terminal to align the aggregated reads to the reference sequence using needle 
    ref_orf: The reference sequence name to which you want to align as it apears in the key of your dictionary
    filepath: The path to the aggregated file or the list that has the path for all the aggregate files
    """
                   
    ref_seq = amplicon_dict[ref_orf]
    
    ref_fasta_path = './Data/Analysis_Novaseq/amplicon_sequences/'+ref_orf+'.fasta '
    
    ## Get the name of the samples
    
    test = filepath.rstrip("aggregate.fa")      

    test2=test.lstrip("./Data/Analysis_Novaseq/aggregate_reads/") 
    
    test3 = test2.rstrip("_")
        
    ## Needle arguments:
    # -auto Turn of prompts
    # -gapopen penalty for opening a gap
    # -asequence one of the two input sequences to align
    # -bsequence the second input sequence to align
    # -aformat output format for the alignment (markx10 is an easily parseable format http://emboss.sourceforge.net/docs/themes/alnformats/align.markx10 )
    # -outfile path to the output file
    
    # For the analysis_wt folder
    needle_out = './Data/Analysis_Novaseq/amplicon_align/'+ ref_orf + '_' + test3 +'.needle'
    
    needle_call = 'needle -auto -gapopen 50 -asequence '+ ref_fasta_path
    
    needle_call += '-bsequence '+ filepath +' -aformat3 markx10 -outfile '+needle_out
    
    print(test3)
        
    subprocess.check_output(needle_call, shell = True)
    
    print("----done----")
          
    return

### Run the alignments looping through the files of size filtered aggregate reads
for filepath in aggregate_list:
    needle_align_on_ref('Novaseq_R67_bacteria', filepath)


# # CONTINUE IN THE 3_Count_variants NOTEBOOK
