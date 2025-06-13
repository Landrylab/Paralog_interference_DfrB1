
# coding: utf-8

# # 004 Sequencing data processing and calculations for selection coeffitient: 
# 
# This script calculates the selection coefficient for each codon in each library according to the following equation (Cisneros et al. 2023):
# 
#     selection_coefficient= (log2((Nmut[t10]/median_Nwt[t10])/(Nmut[t0]/median_Nwt[t0])))/k
#     
#     N=number of reads (abundance that stand for the counts of the mutant by the total number of counts for all the variants in the sample)
#     k=number of generations after t (=10)
# 
# This script assumes the following folder layout
# - Home folder
#     - Scripts folder with this script
#     - Data/Complete_codons dataframes that containing the dataframes with the codon counts and abundances for each sample in each run
# 
# 


## Import modules and check versions

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


# ## CHANGE THE PATH TO YOUR LOCAL FOLDER TO MAKE THE WHOLE NOTEBOOK WORK

os.chdir("/path/Novaseq_data")

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


### Create a directoty to save the selection coeffitient files

sel_coff_prefix = './Analysis/Sel_coff/'
temp_files_prefix= "./intermediate_temporary_dataframes/"

os.makedirs(sel_coff_prefix, exist_ok = True)
os.makedirs(temp_files_prefix, exist_ok = True)

# Define the genetic code
codons = [
  'ATA', 'ATC', 'ATT', 'ATG',
  'ACA', 'ACC', 'ACG', 'ACT',
  'AAC', 'AAT', 'AAA', 'AAG',
  'AGC', 'AGT', 'AGA', 'AGG',
  'CTA', 'CTC', 'CTG', 'CTT',
  'CCA', 'CCC', 'CCG', 'CCT',
  'CAC', 'CAT', 'CAA', 'CAG',
  'CGA', 'CGC', 'CGG', 'CGT',
  'GTA', 'GTC', 'GTG', 'GTT',
  'GCA', 'GCC', 'GCG', 'GCT',
  'GAC', 'GAT', 'GAA', 'GAG',
  'GGA', 'GGC', 'GGG', 'GGT',
  'TCA', 'TCC', 'TCG', 'TCT',
  'TTC', 'TTT', 'TTA', 'TTG',
  'TAC', 'TAT', 'TAA', 'TAG',
  'TGC', 'TGT', 'TGA', 'TGG'
]

residues = [
  'I', 'I', 'I', 'M',
  'T', 'T', 'T', 'T',
  'N', 'N', 'K', 'K',
  'S', 'S', 'R', 'R',
  'L', 'L', 'L', 'L',
  'P', 'P', 'P', 'P',
  'H', 'H', 'Q', 'Q',
  'R', 'R', 'R', 'R',
  'V', 'V', 'V', 'V',
  'A', 'A', 'A', 'A',
  'D', 'D', 'E', 'E',
  'G', 'G', 'G', 'G',
  'S', 'S', 'S', 'S',
  'F', 'F', 'L', 'L',
  'Y', 'Y', '*', '*',
  'C', 'C', '*', 'W'
]

gen_code = {'Codon': codons, 'Residue': residues}
genetic_code = pd.DataFrame(data=gen_code)

### I will try to merge the codon and residues dataframe

### First I need to define the genetic code
codontable_standard = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

### Then I define a function that converts the codon in a cell for the aa that is encoded

# Create a function for translation

def translate_codon(codon):
    try:
        return codontable_standard.get(codon, 'Unknown')
    except:
        return 'Error'
    
    
### Define also a dictionary for the WT sequence of DfrB1

wt_codons = {10: 'AAT',11: 'CCA',12: 'GTT',13: 'GCT',14: 'GGC',15: 'AAT',16: 'TTT',17: 'GTA',18: 'TTC',19: 'CCA',
             20: 'TCG',21: 'GAC',22: 'GCC',23: 'ACG',24: 'TTT',25: 'GGT',26: 'ATG',27: 'GGA',28: 'GAT',29: 'CGC',30: 'GTG',31: 'CGC',32: 'AAG',33: 'AAA',34: 'TCC',
             35: 'GGC',36: 'GCC',37: 'GCC',38: 'TGG',39: 'CAA',40: 'GGT',41: 'CAG',42: 'ATT',43: 'GTC',44: 'GGG',45: 'TGG',46: 'TAC',47: 'TGC',48: 'ACA',49: 'AAT',
             50: 'TTG',51: 'ACC',52: 'CCC',53: 'GAA',54: 'GGC',55: 'TAC',56: 'GCC',57: 'GTC',58: 'GAG',59: 'TCT',60: 'GAG',61: 'GCT',62: 'CAC',63: 'CCA',64: 'GGC',
             65: 'TCA',66: 'GTA',67: 'CAG',68: 'ATT',69: 'TAT',70: 'CCT',71: 'GTT',72: 'GCG',73: 'GCG',74: 'CTT',75: 'GAA',76: 'CGC',77: 'ATC',78: 'AAC'}

### Load the data. Here I can add as many runs I've done for the project 
Run1=pd.read_csv("./Data/Complete_codons_DF_run1.csv")
Run2=pd.read_csv("./Data/Complete_codons_DF_run2.csv")

### Concatenate all the run datafranes into a single one 
Novaseq_data = pd.concat([Run1, Run2]) ## Add here all the runs that have been done in the project
Novaseq_data = Novaseq_data.reset_index(drop=True)

# Apply translation to the dataframe column for library encoded codons
Novaseq_data['Residue'] = Novaseq_data['Codon'].apply(translate_codon)
    
# Apply translation to the dataframe column for WT codons
Novaseq_data['WT_Residue'] = Novaseq_data['WT_Codon'].apply(translate_codon)

print("The length of Novaseq Data " + str(len(Novaseq_data)) + " should be equal to " + str(len(Run1)+len(Run2)))


# Update the data by removing the stop codon at the last position
Novaseq_data = Novaseq_data[Novaseq_data['Position'] != 70]

# Add +9 to the 'Position' column to start from position 10 where actually my DMS starts
Novaseq_data['Position'] = Novaseq_data['Position'] + 9

# Rename a specific column
Novaseq_data = Novaseq_data.rename(columns={'Count': 'read_counts'})
# #Calculate the log2 for the codon counts
Novaseq_data["log2_counts"]=np.log2(Novaseq_data["read_counts"]+1)

# Test if the calculated abundances are correct. The sum of abundances per sample should be equal to 1

### I will test it for t0 because it is easier and all is done in the same way
Sum=Novaseq_data[Novaseq_data['Timepoint'] == 0]

### Remove also lane 2
Sum=Sum[Sum['Lane'] == 1]

Sum=Sum["read_abundance"].sum()

### But the sum has to be divided by all the samples in my dataframe, that in this case is 18 (13 for Run1 and 5 for Run2)

Sum=Sum/18

print("The mean abundance per sample " + str(Sum) + " should be equal to 1. IF NOT, eliminate the read abundance and recalculate it later on.")

# Remove read abunance column 
Novaseq_data=Novaseq_data.drop(columns= ['read_abundance'])

# Save the file 
Novaseq_data.to_csv("./All_run_masterfile.csv", sep = ',', index=False)


Novaseq_data

### I want to do a wide format dataframe separating all the variable columns bu timepoint 0 and 10

# Drop the 'ID' column
Novaseq_data_1 = Novaseq_data.drop('ID', axis=1)

## Split the data into t0 and t10
Novaseq_data_t0 = Novaseq_data_1[Novaseq_data_1['Timepoint'] == 0]
Novaseq_data_t10 = Novaseq_data_1[Novaseq_data_1['Timepoint'] == 10]

# Drop the 'Timepoint' column for each dataframe
Novaseq_data_t0_1 = Novaseq_data_t0.drop(['Timepoint', 'Replicate'], axis=1)
Novaseq_data_t10_1 = Novaseq_data_t10.drop('Timepoint', axis=1)

## Rejoin the dataframe on the common columns 
Novaseq_data_wide= pd.merge(Novaseq_data_t0_1, Novaseq_data_t10_1, on=['Sample', 'Run', 'Date', 'Lane',
       'Expected_Reads', 'Position', 'Codon', 'WT_Codon', 'Residue', 'WT_Residue'], suffixes=('_t0', '_t10'))

# Reset the index
Novaseq_data_wide = Novaseq_data_wide.reset_index(drop=True)

# Reorder the columns
Novaseq_data_wide=Novaseq_data_wide[['Sample', 'Replicate', 'Run', 'Date', 'Lane', 'Expected_Reads', 
                                     'Raw_Reads_t0', 'ng_PCR_for_library_t0', 'Merged_reads_count_R1_R2_t0', 'Percentage_merged_t0', 
                                     'Raw_Reads_t10', 'ng_PCR_for_library_t10', 'Merged_reads_count_R1_R2_t10', 'Percentage_merged_t10',
                                     'Position', 'WT_Codon', 'Codon', 'read_counts_t0','read_counts_t10', 'log2_counts_t0', 'log2_counts_t10',
                                     'WT_Residue','Residue']]


# Save the file 
Novaseq_data_wide.to_csv("./All_run_masterfile_wide.csv", sep = ',', index=False)


### I keep only the data for Lane 1, because in case of having more lanes in the run they have been concatenated

Data = Novaseq_data_wide[Novaseq_data_wide['Lane'] != 2]

print("The lenght of Novaseq Data wide " + str(len(Data)) + " should be equal to " + str((69*64*3*18)))
Data.to_csv("./Working_masterfile_wide.csv", sep = ',', index=False)


### Load Working Masterfile
Data=pd.read_csv("./Working_masterfile_wide.csv")
Novaseq_data=pd.read_csv("./All_run_masterfile.csv")


### I do the filtering at t0 removing all the codons that have less than 100 reads to explore the data and decide where to set the threshold

Data_t0_with_less_100_counts=Data[Data.read_counts_t0 < 100]
Data_t0_with_more_100_counts=Data[Data.read_counts_t0 >= 100]

# Filter the resulting DataFrame for samples with Run == 1 and Run == 2
Data_t0_with_less_100_counts_run1 = Data_t0_with_less_100_counts[Data_t0_with_less_100_counts.Run == 1]
Data_t0_with_less_100_counts_run2 = Data_t0_with_less_100_counts[Data_t0_with_less_100_counts.Run == 2]
Data_t0_with_more_100_counts_run1 = Data_t0_with_more_100_counts[Data_t0_with_more_100_counts.Run == 1]
Data_t0_with_more_100_counts_run2 = Data_t0_with_more_100_counts[Data_t0_with_more_100_counts.Run == 2]

Data_run1 = Data[Data.Run == 1]
Data_run2 = Data[Data.Run == 2]

### Print how much we lose and how much we keep 

print("For Run 1 we lose " + str(len(Data_t0_with_less_100_counts_run1)) + " out of " + str(len(Data_run1)) + " (" + str((len(Data_t0_with_less_100_counts_run1)*100/len(Data_run1)))+ " %) variants")
print("For Run 2 we lose " + str(len(Data_t0_with_less_100_counts_run2)) + " out of " + str(len(Data_run2)) + " (" + str((len(Data_t0_with_less_100_counts_run2)*100/len(Data_run2)))+ " %) variants")


print("In Run 1 the lost reads have on average "+str(Data_t0_with_less_100_counts_run1["read_counts_t0"].mean())+" +/- " + str(Data_t0_with_less_100_counts_run1["read_counts_t0"].std()) + " counts per variant" )
print("In Run 2 the lost reads have on average "+str(Data_t0_with_less_100_counts_run2["read_counts_t0"].mean())+" +/- " + str(Data_t0_with_less_100_counts_run2["read_counts_t0"].std()) + " counts per variant" )


### See also how many variants have more than 75 counts
more75_run1=Data_t0_with_less_100_counts_run1[Data_t0_with_less_100_counts_run1.read_counts_t0 > 75]
more75_run2=Data_t0_with_less_100_counts_run2[Data_t0_with_less_100_counts_run2.read_counts_t0 > 75]

print("For Run 1 " + str(len(more75_run1)) + " out of " + str(len(Data_t0_with_less_100_counts_run1)) + " (" + str((len(more75_run1)*100/len(Data_t0_with_less_100_counts_run1)))+ " %) variants have >75 counts")
print("For Run 2 " + str(len(more75_run2)) + " out of " + str(len(Data_t0_with_less_100_counts_run2)) + " (" + str((len(more75_run2)*100/len(Data_t0_with_less_100_counts_run2)))+ " %) variants have >75 counts")


# # To also help decide where to set the threshold calculate also the Log2 fold change for each sample


# Calculate FC
Data["Log2FC"]=np.log2(Data["read_counts_t10"]/Data["read_counts_t0"])

# Compute the absolute values
Data['Log2FC_abs'] = Data['Log2FC'].abs()

test=Data.copy()

## Divide by run
test_run1 = test[test['Run'] == 1]
test_run2 = test[test['Run'] == 2]

## Separate just one sample
test_run1_CDS = test_run1[test_run1['Sample'] == "CDS_"]

## Separate also for replicate because I want to do first the t0
test_run1_CDS_rep1 = test_run1_CDS[test_run1_CDS['Replicate'] == 1]

print ("The lengt of test_run1_CDS_rep1 " + str(len(test_run1_CDS_rep1)) + " should be equal to " + str(69*64))

print("-----Continue------")


## Calculate all the counts that are in my sample for normalize 
total_counts_t0_CDS=test_run1_CDS_rep1["read_counts_t0"].sum()

# Create a new column with the abundance as the result of read_counts_t0 / total_counts_t0_CDS
test_run1_CDS_rep1['read_abundance_t0'] = test_run1_CDS_rep1['read_counts_t0'] / total_counts_t0_CDS

Sum=test_run1_CDS_rep1["read_abundance_t0"].sum()


print("The abundance per sample " + str(Sum) + " should be equal to 1")
print("-----COMPLETE-----")

## Adding 1 to all the values in the column to avoid NAs

Data["read_counts_t0"]=Data["read_counts_t0"]+1
Data["read_counts_t10"]=Data["read_counts_t10"]+1

print("After correction my dataset has now a minimum value of " + str(Data["read_counts_t0"].min()) + " for t0 and " + str(Data["read_counts_t10"].min()) + " for t10.")

# Function to calculate read_abundance_t0 for each sample and replicate
def calculate_read_abundance_t0(df):
    # Calculate total counts for each combination of Sample and Replicate
    df['total_counts_t0'] = df['read_counts_t0'].sum()

    # Calculate read abundance for each combination of Sample and Replicate
    df['read_abundance_t0'] = df['read_counts_t0'] / df['total_counts_t0']

    return df

# Apply the function to each unique combination of Sample and Replicate
Data = Data.groupby(['Sample', 'Replicate', "Run"]).apply(calculate_read_abundance_t0)

Sum=Data["read_abundance_t0"].sum()

print("The abundance per sample " + str((Sum/3)/18) + " should be equal to 1")
print("-----COMPLETE-----")

Data.to_csv("./All_run_masterfile_wide.csv", sep = ',', index=False)

# Function to calculate also the read_abundance_t0 for each sample and replicate
def calculate_read_abundance_t10(df):
    # Calculate total counts for each combination of Sample and Replicate
    df['total_counts_t10'] = df['read_counts_t10'].sum()

    # Calculate read abundance for each combination of Sample and Replicate
    df['read_abundance_t10'] = df['read_counts_t10'] / df['total_counts_t10']

    return df

# Apply the function to each unique combination of Sample and Replicate
Data = Data.groupby(['Sample', 'Replicate', "Run"]).apply(calculate_read_abundance_t10)

Sum=Data["read_abundance_t10"].sum()

print("The abundance per sample " + str((Sum/3)/18) + " should be equal to 1")
print("-----COMPLETE-----")

Data.to_csv("./All_run_masterfile_wide.csv", sep = ',', index=False)


# # Filter the dataframe on a threshold of a minimum number of counts at t0 to remove all those that are below the threshold to avoid bias because of the limit of foldchange detection. 

# If you have a very low nomber at the initial time point the negative foldchange you can achieve is smaller than if you have a bigger number because the limit of detection is 0.
# If you have a variant that goes from 10 to 1 you will assume the same effect than if you go from 100 to 10, because of the limit of detection, but the first could have a lower or a higher effect but we are just not able to capture it with such a low nomber of input reads at the initial timepoint. Also random sampling effects and non-selective variation within the population can falsely enhance the fitness effect of a low count variant, because the probability of getting lost very early is higher than in high number counts.

### I will set a thershold of at least 100 counts for each codon to be keept in the further processing based on the previous analysis. This limit can be changed at any time.

filtered_data_more_than_100counts_at_t0=Data[Data.read_counts_t0 >= 100]

### Define some functions for ploting the heatmaps to explore the raw data

def generate_codon_heatmap(data):
    data.index = pd.CategoricalIndex(data.index, categories=codon_order)
    data.sort_index(level=0, inplace=True)
    data = data.fillna(0)
    fig, ax = plt.subplots(figsize=(22, 18))
    sns.heatmap(np.log2(data + 1), cmap='viridis', cbar_kws={'label': 'n reads'}, vmax=10, vmin=4, ax=ax)
    ax.set(xlabel='R67 / DfrB1 codon number', ylabel='Codon')
    ax.xaxis.label.set_fontsize(16)
    ax.yaxis.label.set_fontsize(16)
    plt.xticks(rotation=90)
    return fig, ax

def add_wildtype_codons(ax, data, codon_order, wt_codons):
    for codon in range(len(data.columns)):
        x_pos = codon + 0.5
        codon_n = int(list(data.columns)[codon])
        seq = wt_codons[codon_n]
        y_pos = list(data.index).index(seq) + 0.5
        ax.plot(x_pos, y_pos, 'ro')

        
        
def add_wildtype_residues(ax, data, codon_order, wt_codons, aa_order):
    for codon in range(len(data.columns)):
        x_pos = codon + 0.5
        codon_n = int(list(data.columns)[codon])
        seq = wt_codons[codon_n]
        y_pos = aa_order.index(codontable_standard[seq]) + 0.5
        ax.plot(x_pos, y_pos, 'ro')
        
# Arrange codons
codon_order = ['TAA', 'TAG', 'TGA', # *
            'GGA', 'GGC', 'GGG', 'GGT', # G
            'GCA', 'GCC', 'GCG', 'GCT', # A
            'GTA', 'GTC', 'GTG', 'GTT', # v
            'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG', # L
            'ATA', 'ATC', 'ATT', # I
            'ATG', # M
            'TGC', 'TGT', # C
            'CCA', 'CCC', 'CCG', 'CCT', # P
            'TGG', # W
            'TTC', 'TTT', # F
            'TAC', 'TAT', # Y
            'AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT',  # S
            'ACA', 'ACC', 'ACG', 'ACT', # T
            'AAC', 'AAT', # N
            'CAA', 'CAG', # Q
            'CAC', 'CAT', # H
            'AAA', 'AAG', # K
            'AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT', # R
            'GAC', 'GAT', # D
            'GAA', 'GAG' # E
             ]       

aa_order=['*', 'G','A','V', 'L', 'I', 'M', 'C', 'P', 'W', 'F', 'Y', 'S', 'T', 'N', 'Q', 'H', 'K', 'R', 'D', 'E']


Run1_figures = './Analysis/Figures/Run1/'
Run2_figures= "./Analysis/Figures/Run2/"

os.makedirs(Run1_figures, exist_ok = True)
os.makedirs(Run2_figures, exist_ok = True)

### Do the t0 Heatmaps for Run 1 
Run1 = filtered_data_more_than_100counts_at_t0[filtered_data_more_than_100counts_at_t0.Run == 1]

# Get unique sample names from the DataFrame
unique_samples = Run1['Sample'].unique()

# Loop through each sample and apply the code
for sample_name in unique_samples:
    # Prepare the data and plot the heatmaps for the current sample
    data_pivot = Run1[(Run1.Sample == sample_name) & (Run1.Replicate == 1)]
    data_pivot = data_pivot.pivot(index="Codon", columns='Position', values='read_counts_t0')
    
    # Generate the heatmap
    heatmap_fig, ax1 = generate_codon_heatmap(data_pivot)
    
    # Add wildtype codons
    add_wildtype_codons(ax1, data_pivot, codon_order, wt_codons)
    
    # Set the title of the heatmap to the current sample name
    ax1.set_title(sample_name)
    
    # Save the heatmap for each sample
    heatmap_fig.savefig(Run1_figures+f'{sample_name}_heatmap.png')

### Do the t0 Heatmaps for Run 2 
Run2 = filtered_data_more_than_100counts_at_t0[filtered_data_more_than_100counts_at_t0.Run == 2]

# Get unique sample names from the DataFrame
unique_samples = Run2['Sample'].unique()

# Loop through each sample and apply the code
for sample_name in unique_samples:
    # Prepare the data and plot the heatmaps for the current sample
    data_pivot = Run2[(Run2.Sample == sample_name) & (Run2.Replicate == 1)]
    data_pivot = data_pivot.pivot(index="Codon", columns='Position', values='read_counts_t0')
    
    # Generate the heatmap
    heatmap_fig, ax1 = generate_codon_heatmap(data_pivot)
    
    # Add wildtype codons
    add_wildtype_codons(ax1, data_pivot, codon_order, wt_codons)
    
    # Set the title of the heatmap to the current sample name
    ax1.set_title(sample_name)
    
    # Save the heatmap for each sample
    heatmap_fig.savefig(Run2_figures+f'{sample_name}_heatmap.png')

### Save the filtered dataframe 

filtered_data_more_than_100counts_at_t0.to_csv("./Filtered_100at0_wide.csv", sep = ',', index=False)


# #  CALCULATE SELECTION COEFFICIENTS

### Load the data
filtered_data = pd.read_csv("./Filtered_100at0_wide.csv")

# Keep all the rows that on Codon are different from TAG. I will filter out this codon. 
### Stop codon TAG was removed from all datasets because it has been shown to have a lower termination efficiency than the other stop codons in E.coli       
filtered_data = filtered_data[filtered_data.Codon != 'TAG']


filtered_data.to_csv(sel_coff_prefix + "Filtered_100at0_wide_no_TAG.csv", sep = ',', index=False)
filtered_data.columns

### Get the median for the WT

## Keep only the rows for those the Codon is the WT_Codon

WT_medians = filtered_data[filtered_data.Codon == filtered_data.WT_Codon]

WT_medians=WT_medians[['Sample', 'Replicate', 'Run', 'Position', 'WT_Codon', 'Codon', 
                       'read_counts_t0','log2_counts_t0', 'read_abundance_t0',
                       'read_counts_t10', 'log2_counts_t10', 'read_abundance_t10',  
                       'WT_Residue', 'Residue']]
WT_medians=WT_medians.sort_values(by=['Position'])

WT_medians

WT_medians = WT_medians.groupby(['Sample', 'Replicate', "Run"]).agg({'read_abundance_t0': ['median'],'read_abundance_t10': ['median']})

WT_medians.reset_index()

# Rename the columns
WT_medians.columns = pd.MultiIndex.from_tuples([('read_abundance_t0', 'WT_median_t0'),('read_abundance_t10', 'WT_median_t10'),])

WT_medians.columns = WT_medians.columns.droplevel()

filtered_data_WTmedian=pd.merge(filtered_data, WT_medians, on=['Sample', 'Replicate', 'Run'])

# complete_codons_with_WT_no_TAG_WTmedian.columns

print("The length of filtered_data_WTmedian " + str(len(filtered_data_WTmedian)) + " should be equal to " + str(len(filtered_data)))

### Divide the normalized number of reads by the median of WT reads for each position for t0 amd t10
filtered_data_WTmedian["mut_wt_ratio_t0"]=filtered_data_WTmedian["read_abundance_t0"].div(filtered_data_WTmedian["WT_median_t0"])
filtered_data_WTmedian["mut_wt_ratio_t10"]=filtered_data_WTmedian["read_abundance_t10"].div(filtered_data_WTmedian["WT_median_t10"])


## Add a column for the timepoint 
filtered_data_WTmedian["Timepoint"]=10

## Calculate the selection coefficient
filtered_data_WTmedian["Sel_coeff"]=(np.log2(filtered_data_WTmedian["mut_wt_ratio_t10"]/filtered_data_WTmedian["mut_wt_ratio_t0"]))/filtered_data_WTmedian["Timepoint"]

print("The length of filtered_data_WTmedian " + str(len(filtered_data_WTmedian)) + " should be equal to " + str(len(filtered_data)))

## Save the new file with the selection coefficients
filtered_data_WTmedian.to_csv(sel_coff_prefix + "Filtered_100at0_wide_no_TAG_sel_coeff.csv", sep = ',', index=False)

filtered_data_WTmedian = pd.read_csv("./Analysis/Sel_coff/Filtered_100at0_wide_no_TAG_sel_coeff.csv")

