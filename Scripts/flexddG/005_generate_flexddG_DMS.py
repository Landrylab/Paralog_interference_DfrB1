
# # Generate flex_ddG DMS files
# 
# This script produces the input mutation files for a flexddG DMS run.


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
import shutil
import json


# ## Prepare systematic DMS runs for all possible mutations

# Load the data with all the mutations we can simulate
foldx_effects = pd.read_csv('../../Data/Mutational_effects/sciadv.add9109_table_s3.tsv', 
                               sep = '\t')
foldx_effects

## Remove positions from the disordered region and remove unnecessary columns
foldx_effects_noNA = foldx_effects[['Position', 'WT_Residue', 'Residue', 'Mean_ddG_int_HM_A_C']].dropna(subset = ['Mean_ddG_int_HM_A_C'])
foldx_effects_noNA


# Start with the run for NADPH
already_seen = []

## New way to organize the output
main_folder = '/path/output/folder'
os.makedirs(main_folder, exist_ok = True)

tetramer_path = main_folder

chains_move_tetramer = 'A,B,C,D'
count_folders = 0

for index, line in foldx_effects_noNA.iterrows():

    ## Read the mutation
    mut1 = line['WT_Residue'] + str(line['Position']) + line['Residue']
    
    mut1_wt = mut1[0]
    mut1_pos = mut1[1:-1]
    mut1_res = mut1[-1]
        
    ## Prepare files for mutation 1 if it has not been done already and it is not in the disordered region
    if not mut1 in already_seen and int(mut1_pos) >= 21:
    
        count_folders += 1
        
        ## Create the directories for this pair of mutations
        out_path_tetramers_mut = os.path.join(tetramer_path, mut1)
        if not os.path.exists(out_path_tetramers_mut):
            os.makedirs(out_path_tetramers_mut)

        outfile_tetramers_handle = open(os.path.join(out_path_tetramers_mut, 'mutation_list.txt'), 'w')

        ## Rosetta flex ddG input format: 
        # A.F.118.W,B.F.118.W A
        # (chain,WT_res,pos,new_res,chain,WT_res,pos,new_res,...[space][chains_to_move])

        #### Generate the tetramer combinations ####
        ## Tetramer 1 (mut1 in all four chains)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['B', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['C', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['D', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')

        ## Tetramer 2 (mut1 in A; WT in B,C,D)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')

        ## Tetramer 3 (mut1 in B,C,D; WT in A)
        new_line =  '.'.join(['B', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['C', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['D', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')
        
        ## Tetramer 4 (mut1 in A,C; WT in B,D)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['C', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')

        ## Tetramer 5 (mut1 in A,D; WT in B,C)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['D', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')
        
        ## Tetramer 6 (mut1 in A,B; WT in C,D)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['B', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')
        
        ## Tetramer 7 would be WT in all four subunits, which we don't need
        
        ## Add mutations to already_seen
        already_seen.append(mut1)

        ## Close the handles
        outfile_tetramers_handle.close()


# Repeat for DHF binding to the (DfrB1 tetramer + NADPH)

already_seen = []

## New way to organize the output
main_folder = '/path/output/folder'
os.makedirs(main_folder, exist_ok = True)

tetramer_path = main_folder

chains_move_tetramer = 'A,B,C,D,E'

count_folders = 0

for index, line in foldx_effects_noNA.iterrows():

    ## Read the mutation
    mut1 = line['WT_Residue'] + str(line['Position']) + line['Residue']
    
    mut1_wt = mut1[0]
    mut1_pos = mut1[1:-1]
    mut1_res = mut1[-1]
        
    ## Prepare files for mutation 1 if it has not been done already and it is not in the disordered region
    if not mut1 in already_seen and int(mut1_pos) >= 21:
    
        # print(mut1, mut1_wt, mut1_pos, mut1_res)
        
        count_folders += 1

        ## Create the directories for this pair of mutations
        out_path_tetramers_mut = os.path.join(tetramer_path, mut1)
        if not os.path.exists(out_path_tetramers_mut):
            os.makedirs(out_path_tetramers_mut)

        outfile_tetramers_handle = open(os.path.join(out_path_tetramers_mut, 'mutation_list.txt'), 'w')

        ## Rosetta flex ddG input format: 
        # A.F.118.W,B.F.118.W A
        # (chain,WT_res,pos,new_res,chain,WT_res,pos,new_res,...[space][chains_to_move])

        #### Generate the tetramer combinations ####
        ## Tetramer 1 (mut1 in all four chains)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['B', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['C', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['D', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')

        ## Tetramer 2 (mut1 in A; WT in B,C,D)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')

        ## Tetramer 3 (mut1 in B,C,D; WT in A)
        new_line =  '.'.join(['B', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['C', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['D', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')
        
        ## Tetramer 4 (mut1 in A,C; WT in B,D)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['C', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')

        ## Tetramer 5 (mut1 in A,D; WT in B,C)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['D', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')
        
        ## Tetramer 6 (mut1 in A,B; WT in C,D)
        new_line =  '.'.join(['A', mut1_wt, mut1_pos, mut1_res]) + ','
        new_line += '.'.join(['B', mut1_wt, mut1_pos, mut1_res])
        new_line += ' ' + chains_move_tetramer
        outfile_tetramers_handle.write(new_line + '\n')
        
        ## Tetramer 7 would be WT in all four subunits, which we don't need
        
        ## Add mutations to already_seen
        already_seen.append(mut1)

        ## Close the handles
        outfile_tetramers_handle.close()

    

