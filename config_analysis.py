import os
import sys
import pandas as pd
import numpy as np


## -----------------
## Enviroment Config
## -----------------

# set working enviroment -------------------------------

project_data_dir = '/gcm-lfs1/pablo/data/rna_silencing'
project_dir = '/home/pmonteagudo/workspace/silencing_project'

raw_dir = os.path.join(project_data_dir, 'raw')
data_dir = os.path.join(project_data_dir, 'seq_data')

# needed to import all (parastous) code
scripts_dir = os.path.join(project_dir, 'pyRNAdeg')
if scripts_dir not in sys.path:
    sys.path.append(scripts_dir)

# result directories ----------------------------------

#data_results_dir = os.path.join(project_data_dir, 'old_results')
data_results_dir = os.path.join(project_data_dir, 'results')
# results WITH repeat counts normalized by number of locations reads maps to
#data_results_dir = os.path.join(project_data_dir, 'results_with_nh-norm')
# results WITHOUT repeat counts normalized by number of locations reads maps to
#data_results_dir = os.path.join(project_data_dir, 'results_wo_nh-norm')

rna_dir = os.path.join(data_results_dir, 'xp_data/RNA')
chip_dir = os.path.join(data_results_dir, 'xp_data/ChIP')

ratios_dir = os.path.join(data_results_dir, 'Ratios')
if not os.path.exists(ratios_dir):    
    os.makedirs(ratios_dir)

#plots_dir = os.path.join(ratios_dir, 'Plots')
plots_dir = os.path.join(ratios_dir, 'Plots_16-07-21')

if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)
    

## -----------------
## Pipeline Config
## -----------------

# Ignore Data Sets ----------------------------------

#ignore_datasets = []
ignore_datasets = [
    '510_S2-ChIP_1', #?
    '491_S2-ChIP_3', # correlation with other replicates ~0.73 < 0.8
]
#all_samples_df = all_samples_df[~all_samples_df.id.isin(ignore_datasets)]

# Pipeline map dictionaries -------------------------

seq_category = {
    'S2-ChIP':'S2-ChIP', 
    'S5-ChIP':'S5-ChIP',
    'S2-ChIP-INPUT':'INPUT',
    'S2-ChIP-OIN':'INPUT',
    'S2-RIP':'S2-RIP',
    'S5-RIP':'S5-RIP',
    'pA-RNA':'pA-RNA', 
    'total-RNA':'total-RNA',
    'simulated-data':'simulated-data',
    'H3K9me2':'H3K9me2'
}

pipeline_type = {
    'S2-ChIP':'ChIP', 
    'S5-ChIP':'ChIP',
    'S2-ChIP-INPUT':'ChIP',
    'S2-ChIP-OIN':'ChIP',
    'S2-RIP':'RNA', 
    'S5-RIP':'RNA', 
    'pA-RNA':'RNA', 
    'total-RNA':'RNA',
    'simulated-data':'simulated-data',
    'H3K9me2':'ChIP'
}


## -------------
## Mutants Info
## -------------

mut_dict = {
    '1168':'caf1d*ccr4d*',
    '1022':'mot2d',
    '1023':'mot2d',
    '301':'swi6d',
    '302':'clr3d',
    '324':'chp2d',
    '491':'mit1d',
    '504':'rrp6d',
    '510':'caf1d',
    #'523':'unknown_1',
    '523':'unknown',
    #'524':'unknown_2',
    '524':'unknown',
    '530':'exo2d',
    '544':'ccr4d',
    '591':'caf1d',
    '638':'ago1d',
    '80':'clr4d',
    '63':'wt', 
    '65':'wt',
    'WT':'wt',
    'fake-reads':'fake-reads'
}

# it's not bijective
#inv_mut_dict = {v: k for k, v in mut_dict.items()}
inv_mut_dict = {
    'caf1d*ccr4d*':'1168',
    'mot2d':'1022',
    #'mot2d':'1023',
    'swi6d':'301',
    'clr3d':'302',
    'chp2d':'324',
    'mit1d':'491',
    'rrp6d':'504',
    'caf1d':'510', 
    'unknown':'523',
    #'unknown_1':'523',
    #'unknown_2':'524',
    'exo2d':'530',
    'ccr4d':'544',
    #'caf1d':'591',
    'ago1d':'638',
    'clr4d':'80',
    'wt':'WT',
    'fake-reads':'fake-reads'
}

# - Mutant group 1
mutant_group_1 = ["clr4d", "ago1d", "swi6d", "chp2d", "mit1d", "clr3d"]
#[inv_mut_dict.get(key) for key in mutant_group_1]

# - Mutant group 2
mutant_group_2 = ["rrp6d", "exo2d", "caf1d"]
#[inv_mut_dict.get(key) for key in mutant_group_2]

# - Mutant group 3
mutant_group_3 = ["caf1d", "ccr4d", "mot2d"]
#[inv_mut_dict.get(key) for key in mutant_group_3]

# - Mutant group 4 (New sequencing)
mutant_group_4 = ["caf1d", "ccr4d", "mot2d", "caf1d*ccr4d*"]
#[inv_mut_dict.get(key) for key in mutant_group_4]


## -------------------------------------
## Investigate **Heterochromatic genes**
## -------------------------------------

import viz_strands ## get deg1, deg2 and non_degraded

# centromeric genes: `deg1`
old_deg1 = ['dh1', 'dg1']
deg1 = viz_strands.deg1

# subtelomeric genes: `deg2`
old_deg2 = ['SPAC212.11', 'SPAC212.10']
deg2 = viz_strands.deg2

# Mating type region (MTR) gene counts
#deg3 = ['MAT2', 'MAT3', 'MAT1']
deg3 = viz_strands.deg3

# rest of Heterochromatic genes, including mat: `deg3`
non_degraded = viz_strands.non_degraded

all_htc_genes = deg1 + deg2 + non_degraded
htc_genes = deg1 + deg2 + deg3
het_genes = deg1 + deg2 + deg3 + non_degraded


## -------------------------
## Visualization Parameters
## -------------------------

# - Decide wether to include explicitly the **MAT locus genes** in the visualization  
include_mat_locus_vis = True
#include_mat_locus_vis = False

#annotate_plots = True
annotate_plots = False
