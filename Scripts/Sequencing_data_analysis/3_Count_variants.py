
# coding: utf-8

# # 003_sequencing_data_processing: 
# 
# This script finds the mutations in the alignment and count the numers of times each mutation shows up in each library.
# 
# This script assumes the following folder layout
# - Home folder
#     - Scripts folder with this script
#     - Data/Analysis_NovaSeq/amplicon_sequences folder with the reference amplicon fasta file
#     - Data/Analysis_NovaSeq/aggregate_reads folder for aggregated sequences from the done in the previous step
#     - Data/Analysis_NovaSeq/amplicon_align folder for aligned sequences from the done in the previous step
#     - Data/Analysis_NovaSeq folder for the results 

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

def convert_dict_to_aa(codon_dict):
    
    aa_dict_of_dicts = {}
    
    for pos in list(codon_dict.keys()):
        
        aa_dict_of_dicts[pos] = {}
        
        for codon in list(codon_dict[pos].keys()):            
            
            if np.isnan(codon_dict[pos][codon]):
                continue
            
            aa = codontable_standard[codon]
            
            if aa in list(aa_dict_of_dicts[pos].keys()):           
                aa_dict_of_dicts[pos][aa] += codon_dict[pos][codon]
                
            else:
                aa_dict_of_dicts[pos][aa] = codon_dict[pos][codon]
                
    return aa_dict_of_dicts

os.chdir("/path/Novaseq_data")

### Define the reference sequence 

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


# # Define the functions to find mutations and count the variants at each position

def parse_needle_output(needle_align_path):
    
    n_aligns = 0
    align_seqs_dict = {}       
   
    with open(needle_align_path, 'r') as source:

        current_align = ''
        current_qseq = ''
        current_sseq = ''
        qseq_done = 0

        for line in source:

            if line.startswith('>>>') == True:

                n_aligns +=1
                align_name = line.strip('>>>')

                if n_aligns != 1:

                    align_seqs_dict[current_align] = [current_qseq, current_sseq]
                    current_align = align_name
                    current_qseq = ''
                    current_sseq = ''
                    qseq_done = 0

                else:

                    current_align = align_name

            elif line.startswith(';') == False and line.startswith('>') == False and line.startswith('\n') == False and line.startswith('#') == False:

                if qseq_done == 1:
                    current_sseq += line.strip('\n').upper()

                else:
                    current_qseq += line.strip('\n').upper()

            elif line.startswith('#--') == True:

                align_seqs_dict[align_name] = [current_qseq, current_sseq]

            else:

                if qseq_done == 0 and current_qseq != '':
                    qseq_done =1            
                             
    return align_seqs_dict, n_aligns

def find_mutations(path, ref_orf):
    
    allele_dict = {}

    align_dict, align_count = parse_needle_output(path)

    for entry in list(align_dict.keys()):

        read_var_list = []

        query_seq = align_dict[entry][1]
        # aligned prot sequence of the strain

        align_ref = align_dict[entry][0]
        # aligned prot sequence of the reference

        gap_adjust = 0
        # value used to adjust the protein sequence index for the presence of insertions in the strain sequence vs the 
        # reference strain

        backtrack_adjust = 0

        temp_var = None
        # temporary variable to hold the sequence of an insertion or deletion as a string. When the gap ends, annotation 
        # will be added to strain_var_list

        indel_start = 0
        # position start of the indel annotation in the reference sequence, with adjustment for gap presence

        ref_seq_no_gaps = align_ref.replace('-','')
        # Make a copy of the reference sequence without gaps 
        
        align_start = (amplicon_dict[ref_orf].upper().index(ref_seq_no_gaps))+1
        # Look for the starting position of the gapless part of the reference sequence
        # This helps remove gaps at the start and the end of the alignment
        
        query_seq_no_gaps = len(query_seq.replace('-',''))
        # Make a copy of the query sequence without gaps
        
        for nt in range(0, len(align_ref)):
            # iterates through the entire alignment of the strain prot sequence

            if query_seq[nt] == '-':
                # detect a deletion variant

                # logic for indel detection/annotation:
                #
                # suppose we have this alignment  
                #
                # 1 2 3 4 5 6 7 8 9
                # A T - - A A A T G    strain variant: del gaps are indexed because the aa index is based on reference
                # A T K P A - - T G
                # 1 2 3 4 5     6 7    reference: insert gaps not indexed because aa positions do (actually don't?) exist in reference
                #
                # following this logic, every time an insertion is detected and annotated, the gap_adjust value is 
                # incremented by the length of the gap and used to adjust the variant mapping to make it match the 
                # reference index values. The indel aa postion is the first residue detected as part of the indel


                if indel_start == 0:
                    # checks if the character is the start or the continuation of a gap in the alignment

                    temp_var = 'del'+ align_ref[nt]
                    indel_start = (nt+1-gap_adjust)
                    # if it is, starts a new annotation entry with a start position compensated for previous insertions
                    # (if any)

                    backtrack_adjust += 1

                else:

                    temp_var += align_ref[nt]
                    # if it is not, adds the following aa to the deletion annotation

                    backtrack_adjust += 1


            elif align_ref[nt] == '-':
                # detects an insertion variant

                if indel_start == 0:
                    # checks if the character is the start or the continuation of a gap in the alignment

                    temp_var = 'ins'+ query_seq[nt]

                    indel_start = (nt+1-gap_adjust)
                    # if it is, starts a new annotation entry with a start position compensated for previous insertions
                    # (if any)

                    gap_adjust += 1
                    # increments the gap adjust for the this added aa in the strain sequence                   

                else:

                    temp_var += query_seq[nt]
                    # if it is not, adds the following aa to the insertion annotation

                    gap_adjust += 1
                    # increments the gap adjust for the this added aa in the strain sequence

            elif query_seq[nt] != align_ref[nt]:
                # detects a mismatch between the strain sequence and the reference

                variant = align_ref[nt]+'|'+str((nt+1-gap_adjust))+'|'+query_seq[nt]
                read_var_list.append(variant)
                # creates an annotation for the strain-reference aa mismatch and appends it to the list of 
                # annotations

            else:

                 if indel_start != 0:
                    # detects if there is currently an open gap entry. If there is, then the detected mismatch means 
                    # that it has now concluded

                    read_var_list.append(str((indel_start))+temp_var)
                    temp_var = None
                    indel_start = 0
                    # adds the indel annotation to the strain variant list and resets temporary variables for the next 
                    # indel entry

        if len(read_var_list)<25:  
            allele_dict[entry] = read_var_list, align_start
                           
    return allele_dict 

def get_variant_count_1(mutation_set, ref_seq, frag_start, codon_start, n_aa):
    
    variant_abundance_dict ={}
    variants = list(mutation_set.keys())
    codon_groups = {}
    codon = 0
    wt_count =0
    valid_seq=0
    
    # Calculate the end of the fragment with the frag start and the number of residues
    
    frag_end = frag_start + ((n_aa -1 )* 3)
    
    for nt in range(0, (n_aa - 1)*3):
        pos = nt + frag_start
        
        if nt % 3 == 0:
            codon += 1
            
        codon_groups[pos] = codon
        
        variant_abundance_dict[codon] = {}
        
    wt_codons = {}
    
    ref = amplicon_dict[ref_seq].upper()
    
    # Set the abundance of WT codons to not a number (nan) since they would be overrepresented
    
    for aa in range(0, n_aa - 1): 
        
        offset = frag_start - 1
        start = offset+(aa*3)
        wt_codon=ref[start:(start+3)]
        wt_codons[(aa+codon_start)] = wt_codon
        variant_abundance_dict[aa+codon_start][wt_codon]=np.nan
        
    for variant in variants:
        
        var_info = variant.split(',')
        var_count =int(var_info[1].split(';')[1].strip('size='))      
        mut_list = mutation_set[variant][0]
        filtered_list = []

        
        # Only keep variants that appeared in more than 20 reads
        if var_count>=20:
        
            for mutation in mut_list:

                if 'del' in mutation:
                    mut_info = mutation.split('del')
                    mut_pos = int(mut_info[0])
                    
                    if mut_pos == 1:
                        mut_list.remove(mutation)
        
            # If not an indel
            if 'ins' not in str(mut_list) and 'del' not in str(mut_list):
                
                if len(mut_list) ==0:
                    wt_count += var_count
                    
                else:
                    mut_nt_list = []
                    out_list = []
                    
                    for mutation in mut_list:
                        
                        mut_pos = int(mutation.split('|')[1])
                        
                        if mut_pos >= frag_start and mut_pos <= frag_end:
                        
                            # mut_nt_list contains the list of mutated positions inside the coding sequence
                            mut_nt_list.append(codon_groups[mut_pos])
                            
                        else:
                            # out_list contains the list of mutated positions outside the coding sequence
                            out_list.append(mut_pos)
                        
                    if len(set(mut_nt_list)) == 1:
                        valid_seq+=var_count
                        
                        codon = int(list(set(mut_nt_list))[0])
                        
                        wt_seq = wt_codons[codon]
                                                                      
                        new_seq = [x for x in wt_seq]
                        
                        for mutation in mut_list:
                            mut_pos = int(mutation.split('|')[1])
                            mutation = mutation.split('|')[2]
                            
                            codon_pos = (mut_pos-1)%3
                            
                            new_seq[codon_pos] = mutation
                            
                        new_codon = ''.join(new_seq)
                        
                        if new_codon in list(variant_abundance_dict[codon].keys()):

                            variant_abundance_dict[codon][new_codon]+=var_count
                            
                        else:
                            variant_abundance_dict[codon][new_codon]=var_count
  
                    # This is for mutations outside of the coding sequence
                    elif len(set(mut_nt_list)) == 0 and len(out_list)>=1:
                        wt_count+=var_count            
            
        ## For low abundance variants            
        elif var_count < 20:

            for mutation in mut_list:

                if 'del' in mutation:
                    mut_info = mutation.split('del')
                    mut_pos = int(mut_info[0])

                    if mut_pos == 1:
                        mut_list.remove(mutation)
                  
            if len(mut_list) <=3 and 'ins' not in str(mut_list) and 'del' not in str(mut_list):

                if len(mut_list) ==0:
                    wt_count += var_count

                else:
                    #print(mut_list)
                    mut_nt_list = []

                    out_list = []

                    for mutation in mut_list:

                        mut_pos = int(mutation.split('|')[1])

                        if mut_pos >= frag_start and mut_pos <= frag_end:

                            mut_nt_list.append(codon_groups[mut_pos])

                        else:
                            out_list.append(mut_pos)

                    if len(set(mut_nt_list)) == 1:
                        valid_seq+=var_count

                        codon = int(list(set(mut_nt_list))[0])

                        wt_seq = wt_codons[codon]

                        new_seq = [x for x in wt_seq]

                        for mutation in mut_list:
                            mut_pos = int(mutation.split('|')[1])
                            mutation = mutation.split('|')[2]

                            codon_pos = (mut_pos-1)%3

                            new_seq[codon_pos] = mutation

                        new_codon = ''.join(new_seq)
                        
                        if new_codon in list(variant_abundance_dict[codon].keys()):
                            
                            variant_abundance_dict[codon][new_codon]+=var_count
  
                        else:
                            variant_abundance_dict[codon][new_codon]=var_count
                    
    return variant_abundance_dict, wt_count, wt_codons


# # Define a new function to count the variants normalizing by the total number of reads in each library


def get_timepoint_fraction_df(needle_file, Sample_ID):
    
    # Find mutations (this function parses the needle file)
    pool_1_muts = find_mutations(needle_file, 'Novaseq_R67_bacteria')
    
    frag_1_dict = get_variant_count_1(pool_1_muts, 'Novaseq_R67_bacteria', 1, 1, 71)
    mut_dict = frag_1_dict[0]
    wt_count = frag_1_dict[1]
        
    variant_df = pd.DataFrame(mut_dict)
    
    array_size=len(variant_df.to_numpy().flatten())
    
    # Calculate the total number of reads before adding the WT count to the matrix
    read_total = variant_df.sum().sum()+wt_count+array_size
    print(Sample_ID, read_total, wt_count, array_size)
    
    # Extract the dictionary with the WT codons
    wt_codons = frag1_dict[2]
    
    # Add a loop to fill in the WT codons
    for position, wt_codon in wt_codons.items():
        mut_dict[position][wt_codon] = wt_count
    
    # Update variant_df
    variant_df = pd.DataFrame(mut_dict)
    variant_df_no_NaN = variant_df.fillna(0)
    
    variant_df_no_NaN = variant_df_no_NaN + 1
        
    read_fraction_df = variant_df_no_NaN/(read_total)
    print(read_total, 1/read_total)
    
    read_fraction_df.rename_axis(read_total, inplace=True)
    
    df_out_path='./Data/Analysis_Novaseq/read_abundances/Codons/'+str(Sample_ID)+'_read_frac.csv'
    
    read_fraction_df.to_csv(df_out_path, sep=',')
    
    # Convert dataframe to the residue level
    aa_frag_1 = convert_dict_to_aa(mut_dict)
    variant_df_aa = pd.DataFrame(aa_frag_1)
    
    # Get proportion of reads at the residue level
    array_size=len(variant_df_aa.to_numpy().flatten())
    read_total = variant_df_aa.sum().sum() + wt_count + array_size
    
    print(Sample_ID, read_total, wt_count, array_size)
    
    variant_df_aa_no_NaN = variant_df_aa.fillna(0)
    variant_df_aa_no_NaN = variant_df_aa_no_NaN + 1
    
    # Normalize by the total of reads
    read_fraction_df_aa = variant_df_aa_no_NaN/(read_total)
    print(read_total, 1/read_total)
    
    read_fraction_df_aa.rename_axis(read_total, inplace=True)
    
    # Save file
    df_out_path='./Data/Analysis_Novaseq/read_abundances/Residues/'+str(Sample_ID)+'_read_frac.csv'
    read_fraction_df_aa.to_csv(df_out_path, sep=',')
    
    # variant_df has the raw counts for each codon
    # read_fraction_df has counts normalized by the total of reads
    return read_fraction_df, variant_df

os.makedirs('./Data/Analysis_Novaseq/read_abundances/', exist_ok = True)
os.makedirs('./Data/Analysis_Novaseq/read_abundances/Codons/', exist_ok = True)
os.makedirs('./Data/Analysis_Novaseq/read_abundances/Residues/', exist_ok = True)

# The loop to calculate read proportions in each sample

for file in dict_path_to_align_files.values():
    needle_file=file

    # Set the sample name
    test = file.rstrip(".needle")
    sample_name=test.lstrip("./Data/Analysis_Novaseq/amplicon_align/Novaseq_R67_bacteria/")
    
    print(sample_name)
    
    read_fraction_df, variant_df = get_timepoint_fraction_df(needle_file, sample_name)
        
    print('-----done----')


# # Generate a master dataframe that contains all the information


# Define the genetic code and get the sample information from the excel sheet
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


### Get the sample information from the excel sheet
metada= pd.read_excel('./Data/Novaseq_Sample_Info_df.xlsx', index_col=1) 
metada=metada.drop(columns=['Sample'])
metada.head()

### Create a directoty to save the selection coeffitient files

sel_coff_prefix = './Data/Analysis_Novaseq/Sel_coff/'
os.makedirs(sel_coff_prefix, exist_ok = True)


## Make a dictionary to identify the different files

path_to_codon_read_abundances = "./Data/Analysis_Novaseq/read_abundances/Codons/"

dict_path_to_codons_read_abundances = {}

path_to_residues_read_abundances = "./Data/Analysis_Novaseq/read_abundances/Residues/"

dict_path_to_residues_read_abundances = {}


for filename in os.listdir(path_to_codon_read_abundances):
    
    f = os.path.join(path_to_codon_read_abundances, filename)
    
    f1 = f.rstrip("_read_frac.csv")
    ### I want the sample name to be the same as the Sample.1 in the excel dataframe
    
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/read_abundances/")
    
    sample_name2=sample_name.lstrip("Codons")
    
    sample_name3=sample_name2.lstrip("/")
    
    dict_path_to_codons_read_abundances[sample_name3] = f

dict_path_to_codons_read_abundances

for filename in os.listdir(path_to_residues_read_abundances):
    
    f = os.path.join(path_to_residues_read_abundances, filename)
    
    f1 = f.rstrip("_read_frac.csv")
       
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/read_abundances/")
    
    sample_name2=sample_name.lstrip("Residues")
    
    sample_name3=sample_name2.lstrip("/")
    
    dict_path_to_residues_read_abundances[sample_name3] = f
    

## Make also a dictionary to identify the different files with the counts instead of the fractions

path_to_codon_counts = "./Data/Analysis_Novaseq/aggregate_dataframes/Codons/"

dict_path_to_codons_counts = {}

path_to_residues_counts = "./Data/Analysis_Novaseq/aggregate_dataframes/Residues/"

dict_path_to_residues_counts = {}


for filename in os.listdir(path_to_codon_counts):
    
    f = os.path.join(path_to_codon_counts, filename)
    
    f1 = f.rstrip("_codons.txt")
    ### I want the sample name to be the same as the Sample.1 in the excel dataframe
    
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/aggregate_dataframes/")
    
    sample_name2=sample_name.lstrip("Codons")
    
    sample_name3=sample_name2.lstrip("/Novaseq_R67_bacteria_")
    
    dict_path_to_codons_counts[sample_name3] = f
    
for filename in os.listdir(path_to_residues_counts):
    
    f = os.path.join(path_to_residues_counts, filename)
    
    f1 = f.rstrip("_residues.txt")
       
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/aggregate_dataframes/")
    
    sample_name2=sample_name.lstrip("Residues")
    
    sample_name3=sample_name2.lstrip("/Novaseq_R67_bacteria_")
    
    dict_path_to_residues_counts[sample_name3] = f
    


### Get the codon abundances for each sample and format the dataframe into the correct 

codon_abundances = './Data/Analysis_Novaseq/Sel_coff/Codons/'
os.makedirs(codon_abundances, exist_ok = True)

all_codon_data=pd.DataFrame()

for file in dict_path_to_codons_read_abundances.values():
    
    ### Import the file
    codon_df=pd.read_csv(file)
    
    ### Change the first column name
    codon_df.rename(columns = {list(codon_df)[0]:'Codon'}, inplace=True)
    
    ### Unpivot the codon_df to a longer dataframe in wich each row is the read_abundance for each codon at each position
    new_codon= pd.melt(codon_df, id_vars="Codon")
    new_codon.rename({'variable': 'Position', 'value': 'read_abundance'}, axis=1, inplace=True)

    ### Get the sample ID as it is in metadata
    f1 = file.rstrip("_read_frac.csv")
        
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/read_abundances/")
         
    sample_name2=sample_name.lstrip("Codons")
    
    sample_name3=sample_name2.lstrip("/")

    
    ### Add a new column with the sample ID to the codon df
    
    new_codon['Sample.1'] = sample_name3
    
    
    ### Save the file
    new_codon.to_csv(codon_abundances + sample_name3 + '_codons.csv', sep = ',', index=False)

    ### I need to concatenate all the dataframes into the same master dataframe "all_codon_data"
    all_codon_data=all_codon_data.append(new_codon,ignore_index=True)
    
### Save the file master file 
all_codon_data.to_csv(sel_coff_prefix + "All_codon_abundance.csv", sep = ',', index=False)    

### Get the residues abundances for each sample and format the dataframe into the correct format

Residues_abundances = './Data/Analysis_Novaseq/Sel_coff/Residues/'
os.makedirs(Residues_abundances, exist_ok = True)

all_residues_data=pd.DataFrame()

for file in dict_path_to_residues_read_abundances.values():
    
    ### Import the file
    residues_df=pd.read_csv(file)
    
    ### Change the first column name
    residues_df.rename(columns = {list(residues_df)[0]:'Residue'}, inplace=True)
    
    ### Unpivot the residues_df to a longer dataframe in wich each row is the read_abundance for each codon at each position
    new_residue= pd.melt(residues_df, id_vars="Residue")
    new_residue.rename({'variable': 'Position', 'value': 'read_abundance'}, axis=1, inplace=True)

    ### Get the sample ID as it is in metadata
    f1 = file.rstrip("_read_frac.csv")
        
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/read_abundances/")
         
    sample_name2=sample_name.lstrip("Residues")
    
    sample_name3=sample_name2.lstrip("/")
    
    ### Add a new column with the sample ID to the codon df
    
    new_residue['Sample.1'] = sample_name3
    
    
    ### Save the file
    new_residue.to_csv(Residues_abundances + sample_name3 + '_residues.csv', sep = ',', index=False)

    ### I need to concatenate all the dataframes into the same master dataframe "all_codon_data"
    all_residues_data=all_residues_data.append(new_residue,ignore_index=True)
    
### Save the file master file 
all_residues_data.to_csv(sel_coff_prefix + "All_residues_abundance.csv", sep = ',', index=False)    


### Get the codon counts and the residues counts for each sample and format the dataframe into the correct f

codon_abundances = './Data/Analysis_Novaseq/Sel_coff/Codons/'
os.makedirs(codon_abundances, exist_ok = True)

all_codon_counts=pd.DataFrame()

for file in dict_path_to_codons_counts.values():
    
    ### Import the file
    codon_df=pd.read_csv(file, sep="\t")
    
    ### Change the first column name
    codon_df.rename(columns = {list(codon_df)[0]:'Codon'}, inplace=True)
    
    ### Unpivot the codon_df to a longer dataframe in wich each row is the read_abundance for each codon at each position
    new_codon= pd.melt(codon_df, id_vars="Codon")
    new_codon.rename({'variable': 'Position', 'value': 'Count'}, axis=1, inplace=True)

      
    ### Get the sample ID as it is in metadata
    
    f1 = file.rstrip("_codons.txt")
   
    
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/aggregate_dataframes/")
    
    sample_name2=sample_name.lstrip("Codons")
    
    sample_name3=sample_name2.lstrip("/Novaseq_R67_bacteria_")
    
    ### Add a new column with the sample ID to the codon df
    
    new_codon['Sample.1'] = sample_name3
    
    
    ### Save the file
    new_codon.to_csv(codon_abundances + sample_name3 + '_codons.csv', sep = ',', index=False)

    ### I need to concatenate all the dataframes into the same master dataframe "all_codon_data"
    all_codon_counts=all_codon_counts.append(new_codon,ignore_index=True)

    
    
all_residues_counts=pd.DataFrame()

for file in dict_path_to_residues_counts.values():
    
    ### Import the file
    residues_df=pd.read_csv(file, sep="\t")
    
    ### Change the first column name
    residues_df.rename(columns = {list(residues_df)[0]:'Residue'}, inplace=True)
    
    ### Unpivot the residues_df to a longer dataframe in wich each row is the read_abundance for each codon at each position
    new_residue= pd.melt(residues_df, id_vars="Residue")
    new_residue.rename({'variable': 'Position', 'value': 'Count'}, axis=1, inplace=True)

    ### Get the sample ID as it is in metadata
    f1 = file.rstrip("_residues.txt")
        
    sample_name=f1.lstrip("./Data/Analysis_Novaseq/aggregate_dataframes/")
         
    sample_name2=sample_name.lstrip("Residues")
    
    sample_name3=sample_name2.lstrip("/Novaseq_R67_bacteria_")

    ### Add a new column with the sample ID to the codon df
    
    new_residue['Sample.1'] = sample_name3
    
    
#     ### Save the file
    new_residue.to_csv(Residues_abundances + sample_name3 + '_residues.csv', sep = ',', index=False)

    ### I need to concatenate all the dataframes into the same master dataframe "all_codon_data"
    all_residues_counts=all_residues_counts.append(new_residue,ignore_index=True)
# ### Save the file master file 
all_residues_counts.to_csv(sel_coff_prefix + "all_residues_counts.csv", sep = ',', index=False)    
all_codon_counts.to_csv(sel_coff_prefix + "all_codon_counts.csv", sep = ',', index=False)    

### Merge the counts and the abundance tables into one dataframe 

all_codons_counts_and_abundances=pd.merge(all_codon_data, all_codon_counts, on=['Codon', 'Position', 'Sample.1'])
all_residues_counts_and_abundances=pd.merge(all_residues_data, all_residues_counts, on=['Residue', 'Position', 'Sample.1'])

### Save the master file 
all_residues_counts_and_abundances.to_csv(sel_coff_prefix + "all_residues_counts_and_abundances.csv", sep = ',', index=False)    
all_codons_counts_and_abundances.to_csv(sel_coff_prefix + "all_codons_counts_and_abundances.csv", sep = ',', index=False) 

### Merge the dataframes with the metadata to include the sample information to the read abundance information

complete_codons=pd.merge(all_codons_counts_and_abundances, metada, on="Sample.1")
complete_residues=pd.merge(all_residues_counts_and_abundances, metada, on="Sample.1")

### Add the WT information for each position 
WT_matrix=pd.read_csv("./Data/Analysis_Novaseq/Sel_coff/WT_sequence_table.txt", sep="\t")

### Convert the position column for all the data frames into a numeric object
WT_matrix["Position"] = pd.to_numeric(WT_matrix["Position"])

complete_codons["Position"] = pd.to_numeric(complete_codons["Position"])

complete_residues["Position"] = pd.to_numeric(complete_residues["Position"])

### Merge the files 
complete_codons_with_WT=pd.merge(complete_codons, WT_matrix, on="Position")
complete_residues_with_WT=pd.merge(complete_residues, WT_matrix, on="Position")

complete_codons_with_WT['Sample'] = complete_codons_with_WT['ID'].str.slice(0,4)
complete_residues_with_WT['Sample'] = complete_residues_with_WT['ID'].str.slice(0,4)

### Reorder the columns to 

codon_cols=["ID", 'Sample', 'Run', 'Date', 'Timepoint', 'Replicate', 'Lane', 'Expected_Reads', 'Raw_Reads', 'ng_PCR_for_library',
       'Merged_reads_count_R1_R2', 'Percentage_merged', 'Position','Codon', 'WT_Codon', 'Count', 'read_abundance']
complete_codons_with_WT = complete_codons_with_WT.reindex(columns=codon_cols)

res_cols=["ID", 'Sample', 'Run', 'Date', 'Timepoint', 'Replicate', 'Lane', 'Expected_Reads', 'Raw_Reads', 'ng_PCR_for_library',
       'Merged_reads_count_R1_R2', 'Percentage_merged', 'Position','Residue', 'WT_Residue', 'Count', 'read_abundance']
complete_residues_with_WT = complete_residues_with_WT.reindex(columns=res_cols)

### Save the dataframes
complete_codons_with_WT.to_csv(sel_coff_prefix + "Complete_codons_DF.csv", sep = ',', index=False)
complete_residues_with_WT.to_csv(sel_coff_prefix + "Complete_residues_DF.csv", sep = ',', index=False)

### Save the dataframe also in the masterfolder wher I will put together all the runs and continue the analysis

complete_codons_with_WT.to_csv("../Novaseq_combined/Data/Complete_codons_DF_run2.csv", sep = ',', index=False)
complete_residues_with_WT.to_csv("../Novaseq_combined/Data/Complete_residues_DF_run2.csv", sep = ',', index=False)


# # CONTINUE IN home/Novaseq_combined/ 4_Calculate_sel_coff NOTEBOOK 
