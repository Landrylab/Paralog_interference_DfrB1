
# # Organize Rosetta output
# 
# This script will format the output of flexddG for our analysis script.

# Load libraries
import numpy as np
import pandas as pd
import matplotlib
import csv
import os
import sys
import subprocess
from Bio.PDB import PDBParser, PDBIO
import glob
from Bio import SeqIO
import re
import json
from collections import OrderedDict
import seaborn as sns


# ## Format the DMS data for NADPH

## Load file
dms_data_nadph = pd.read_csv('/path/results/NADPH', 
                       sep = ',')

dms_data_nadph


all_data_ddg_dms_NADPH = []

## Make sure the mutation label column indicates which subunits are a WT copy
new_list = []

## Loop through the rows
for index, row in dms_data_nadph.iterrows():

    mut_dict = {'A': 'WT', 'B':'WT', 'C':'WT', 'D':'WT'}

    mutation = row['mutation']    
    mut_list = mutation.split(':')

    for new_mut in mut_list:
        mut_dict[new_mut[0]] = ''.join(new_mut.split('.')[1:])

    mut_label = mut_dict['A'] + '_' + mut_dict['B'] + '_' + mut_dict['C'] + '_' + mut_dict['D']

    new_row = [mutation, mut_label] + list(row[2:])
    new_list.append(new_row)

all_data_ddg_dms_NADPH = pd.DataFrame(new_list, columns = dms_data_nadph.columns)
all_data_ddg_dms_NADPH


all_data_ddg_dms_NADPH.to_csv('../../Data/Mutational_effects/all_results_DMS_NADPH_formatted.tsv', 
                             index = False, sep = '\t')


# ## Format the DMS data for DHF binding after NADPH

## Load file
dms_data_nadph_dhf = pd.read_csv('/path/flexddG/NADPH+DHF'), 
                       sep = ',')

dms_data_nadph_dhf

all_data_ddg_dms_NADPH_DHF = []

## Make sure the mutation label column indicates which subunits are a WT copy
new_list = []

## Loop through the rows
for index, row in dms_data_nadph_dhf.iterrows():

    mut_dict = {'A': 'WT', 'B':'WT', 'C':'WT', 'D':'WT'}

    mutation = row['mutation']    
    mut_list = mutation.split(':')

    for new_mut in mut_list:
        mut_dict[new_mut[0]] = ''.join(new_mut.split('.')[1:])

    mut_label = mut_dict['A'] + '_' + mut_dict['B'] + '_' + mut_dict['C'] + '_' + mut_dict['D']

    new_row = [mutation, mut_label] + list(row[2:])
    new_list.append(new_row)

all_data_ddg_dms_NADPH_DHF= pd.DataFrame(new_list, columns = dms_data_nadph_dhf.columns)
all_data_ddg_dms_NADPH_DHF

all_data_ddg_dms_NADPH_DHF.to_csv('../../Data/Mutational_effects/all_results_DMS_NADPH_DHF_formatted.tsv', 
                             index = False, sep = '\t')

