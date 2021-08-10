import os
import sys
import pandas as pd
import numpy as np

## ------
## Config
## ------

project_data_dir = '/gcm-lfs1/pablo/data/rna_silencing'
project_dir = '/home/pmonteagudo/workspace/silencing_project'

scripts_dir = os.path.join(project_dir, 'pyRNAdeg')
if scripts_dir not in sys.path: 
    sys.path.append(scripts_dir)

# result directories ----------------------------------

#data_results_dir = os.path.join(project_data_dir, 'old_results')
data_results_dir = os.path.join(project_data_dir, 'results')

rna_dir = os.path.join(data_results_dir, 'xp_data/RNA')
chip_dir = os.path.join(data_results_dir, 'xp_data/ChIP')

ratios_dir = os.path.join(data_results_dir, 'Ratios')

#plots_dir = os.path.join(ratios_dir, 'Plots')
plots_dir = os.path.join(ratios_dir, 'Plots_v2')

if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

## -----------
## Parameters:
## -----------

# - Decide wether to include explicitly the **MAT locus genes** in the visualization  
include_mat_locus_vis = True
#include_mat_locus_vis = False

## -------------------------------------
## Investigate **Heterochromatic genes**
## -------------------------------------

import viz_strands ## get deg1, deg2 and non_degraded

## centromeric genes: `deg1`
old_deg1 = ['dh1', 'dg1']
deg1 = viz_strands.deg1

## subtelomeric genes: `deg2`
old_deg2 = ['SPAC212.11', 'SPAC212.10']
deg2 = viz_strands.deg2

# Mating type region (MTR) gene counts
#deg3 = ['MAT2', 'MAT3', 'MAT1']
deg3 = viz_strands.deg3

## rest of Heterochromatic genes, including mat: `deg3`
non_degraded = viz_strands.non_degraded

all_htc_genes = deg1 + deg2 + non_degraded
htc_genes = deg1 + deg2 + deg3
het_genes = deg1 + deg2 + deg3 + non_degraded
