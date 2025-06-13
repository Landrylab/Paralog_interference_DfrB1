
import os
import pandas as pd
import plotnine as p9
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import DataFrame


## Calculate the relative growth

## Load the data and make a new dataframe with a subset of data removing all unnecesary columns
results_clean = pd.read_csv("../../Data/Growth_recovery_duplication/Data_FigS2B.csv", delimiter = ",", header=0)
data = results_clean[["sample","REP", "conc_ara", "tmp", "auc_l"]]
data.to_csv('./small_data_set_GR_calc.csv') 

## Group the data in each condiction and calculate the medians and the standart deviations

test = data.groupby(["sample", "conc_ara", "tmp"])["auc_l"].agg([np.median, np.std])
test.rename(columns = {'median':'AUC_median','std':'AUC_std' }, inplace = True)
test.to_csv('./test.csv') 

##Load the we trimmed and group them by sample

test = pd.read_csv('./test.csv')

group_samples = test.groupby('sample')

df_sample_list = [] ## Create an empty to make a data frame with all the variables we need in our equation

auc_10_0 = 0
auc_0_0 = 0


## Define our demominatior, which will be constant for each sample
for sample_name, df_sample in group_samples:
    for index, row in df_sample.iterrows():
        if row['conc_ara'] == 0.000 and row['tmp'] == 0:
            auc_0_0 = row['AUC_median']
        elif row['conc_ara'] == 0.000 and row['tmp'] == 10:
            auc_10_0 = row['AUC_median']
        else:
            pass
        
    df_sample['AUC_10_0'] = [auc_10_0]*len(df_sample.index)
    df_sample['AUC_0_0'] = [auc_0_0]*len(df_sample.index)
    
    df_sample_list.append(df_sample)
    
df_sample_final = pd.concat(df_sample_list)


df_sample_conc_list = []

group_samples_conc = df_sample_final.groupby(['sample', 'conc_ara'])

for sample_name, df_sample in group_samples_conc:
    for index, row in df_sample.iterrows():
        if row['tmp'] == 0:
            auc_0_0 = row['AUC_median']
        elif row['tmp'] == 10:
            auc_10_0 = row['AUC_median']
        else:
            pass
        
    df_sample['AUC_10_ara'] = [auc_10_0]*len(df_sample.index)
    df_sample['AUC_0_ara'] = [auc_0_0]*len(df_sample.index)
    
    df_sample_conc_list.append(df_sample)

df_sample_conc_final = pd.concat(df_sample_conc_list)

df_final = df_sample_conc_final[df_sample_conc_final['tmp'] == 0].copy() 
df_final = df_final[['sample', 'conc_ara',  'AUC_median', 'AUC_10_0', 'AUC_0_0', 'AUC_10_ara', 'AUC_0_ara']].copy()

df_final['recovery_growth'] = 100 * (1 - ((df_final['AUC_10_ara'] - df_final['AUC_0_ara']) / (df_final['AUC_10_0'] - df_final['AUC_0_0'])))

### This function defines how to label each point in the scatter plot

def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+0.0005, point['y'], str(point['val']))


g = sns.catplot(data=df_final[df_final["sample"].isin(["WT&CDS"])], 
          x="conc_ara", y="recovery_growth", col = "sample", hue="conc_ara", kind="point", errorbar="AUC_std")

plt.xticks(rotation = 45)
g.set(xlabel='% Arabinose', ylabel='Growth Recovery', title="WT & CDS")
plt.axhline(50, color = "b", linestyle='--')
plt.axhline(100, color = "black", linestyle='--')
plt.grid(alpha=0.5,linestyle='--' )

plt.savefig(r"../Figures/Supp_figures/FigS2B", format='pdf', dpi=300)


g = sns.catplot(data=df_final[df_final["sample"].isin(["WT&EMPTY"])], 
                x="conc_ara", y="recovery_growth", col = "sample", hue="conc_ara", kind="point", errorbar="AUC_std")

plt.xticks(rotation = 45)
g.set(xlabel='% Arabinose', ylabel='Growth Recovery', title="WT & EMPTY")
plt.axhline(50, color = "b", linestyle='--')
plt.axhline(100, color = "black", linestyle='--')
plt.grid(alpha=0.5,linestyle='--' )

# plt.savefig(r"./Figures/Supp_figures/growth_recovery_WT_Singleton.pdf", format='pdf', dpi=300)

