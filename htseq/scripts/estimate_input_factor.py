#!/usr/bin/env python

import itertools
import time
import os
import sys
import argparse

import numpy as np
import pandas as pd

import HTSeq
import pysam

import matplotlib.pyplot as plt
import seaborn as sns

import collections

from sklearn.linear_model import LinearRegression
from sklearn import metrics

import coverage as cov
import gene_counts as gc


## -------------------
## Auxiliary Functions
## -------------------

def get_args():

    parser = argparse.ArgumentParser()
    
    # from 'coverages'
    parser.add_argument(
        '--in_bam', '-f',
        type=str, default=None,
        help='Input bam file for which coverage is to be computed.')
        
    parser.add_argument(
        '--in_chip_input', '-i',
        type=str, default=None,
        help='Input ChIP INPUT bam file which coverage scaled by a norm factor needs to be substracted from coverage.')
    
    # # from 'counts'
    # # => use counts if have already been computed
    # parser.add_argument(
    #     '--in_count',
    #     type=str, default=None,
    #     help='Count file for ChIP data, used to estimate the norm factor (lambda) for INPUT subtraction.')
    # 
    # # => use counts if have already been computed
    # parser.add_argument(
    #     '--in_chip_input_count',
    #     type=str, default=None,
    #     help='Count file for ChIP INPUT data, used to estimate the norm factor (lambda) for INPUT subtraction. INPUT counts scaled by a norm factor need to be subtracted from ChIP counts.')

    parser.add_argument(
        '--in_gdf', '-d', type=str,
        help='Gene information table (gdf), containing gene_lengths, etc.', required=True)
    
    parser.add_argument(
        '--verbose', '-v',
        type=bool, default=True,
        help='Verbose (default: True).')
    
    parser.add_argument(
        '--out_dir', '-o',
        type=str, default='.',
        help='Save results in (default: current working directory).')
        
    parser.add_argument(
        '--feature_type',
        default=[],
        action='append', help="Feature or Features to be parse from GFF and used for counting - currently, only used for naming the output_files")

    parser.add_argument(
        '--max_nh',  type=int,
        default=16,
        help="Potentially allow multimapped-reads with less than 'max_nh' (default:16) hits - Note: we only allow this if the read maps to a repeat feature.")
        
    parser.add_argument(
        "--count_mode", "-m", dest="count_mode",
        choices=("union", "intersection-strict", "intersection-nonempty"),
        default="union",
        help="Mode to handle reads overlapping more than one feature " +
             "(choices: union, intersection-strict, intersection-nonempty; default: union).")

    parser.add_argument(
        "--ambiguous_assignment_mode", dest="ambiguous_assignment_mode", type=str,
        choices=("none", "all", "fractional"),
        default="none",
        help="Mode to handle ambiguous reads (reads compatible with several overlapping features)" +
             "(choices: none, all, fractional; default: none).")

    parser.add_argument(
        "--multimapped_mode", dest="multimapped_mode", type=str,
        choices=("none", "all", "fractional", "ignore_secondary"),
        default="none",
        help="Mode to handle multimapped read alignments (reads that are not uniquely aligned to one location)" +
             "(choices: none, all, fractional, ignore_secondary; default: none).")



    args = parser.parse_args()

    return args


## ---------
## Functions
## ---------

def raw_counts_centromer(bam_file, in_gdf, features_list, stranded=False, count_mode="union", ambiguous_assignment_mode="none", multimapped_mode="none", max_nh=16):
  
  # get sample_file
  sample_file = os.path.basename(bam_file)
  sample_id = sample_file.split(".")[0]
  
  # -------------------------------------------------------------
  # ---------    Feature Annotation Data: in_gdf      -----------
  # -------------------------------------------------------------
  
  # Import Gene Data Table (gdf)
  gdf = pd.read_csv(in_gdf, sep='\t')
  
  ## - Get 'repeat_ids' for using during counting
  repeat_ids = gdf[gdf['category'] == 'repeat']['gene_id'].tolist()
  
  ## Filter for 'features' of interest
  # TODO: too specific
  feature_ids = [xx['id'] for xx in features_list]
  gdf = gdf[gdf['gene_id'].isin(feature_ids)].reset_index(drop=True)
  #gdf = gdf[gdf['category'] == 'ip_region'].reset_index(drop=True)
  
  feature_type = gdf['type'].unique()
  
  print("\nGetting `raw_gene_counts` using feature: {}".format(feature_type))
  print("Using `stranded` mode: {}".format(stranded))
  print("Using Count Parameters: \n     - `count_mode`: {}\n     - `ambiguous_assignment_mode`: {}\n     - `multimapped_mode`: {}\n".format(count_mode, ambiguous_assignment_mode, multimapped_mode))
  
  ## Gene Data Info: (currently, seems more practical to just use this)
  # => features: as `GenomicArrayOfSets` htseq object.
  features, genes_dict = gc.get_features_from_gdf(gdf, stranded=stranded) ## use same feature as Parastou
  
  # -------------------------------------------------------------
  # -----------         Get feature counts          -------------
  # -------------------------------------------------------------
  
  count_df = []
  # Parse `.bam` file and count 'reads' falling into each `feature`s genomic location
  for kk in genes_dict.keys():
    
    df = gc.get_counts(
        bam_file,
        features,
        repeat_ids,
        filter_region = genes_dict[kk],
        #stranded = False, # taken care in features `GenomicArrayOfSets`
        count_mode = count_mode,
        ambiguous_assignment_mode = ambiguous_assignment_mode,
        multimapped_mode = multimapped_mode,
        max_nh = max_nh
        )
    # remove summary entries from `count_df`
    df = df[~df['gene_id'].str.startswith('_')]
    count_df.append(df)
  
  count_df = pd.concat(count_df).reset_index(drop=True)
  
  # raw Gene Counts Table: (gxt)
  try:
      count_df = pd.merge(count_df, gdf[['gene_id', 'gene_name', 'transcript_length', 'gene_length', 'type', 'category']], on='gene_id', how='outer')
      #gene_counts_df = pd.merge(gene_counts_df, gdf[['gene_id', 'gene_name', 'transcript_length', 'gene_length', 'type', 'category']], on='gene_id', how='outer')
  except:
      gdf_cols = ['gene_id', 'gene_name', 'type', 'category']
      gdf_cols.extend([cc for cc in gdf.columns if 'length' in cc])
      count_df = pd.merge(count_df, gdf[gdf_cols], on='gene_id', how='outer')
      #gene_counts_df = pd.merge(gene_counts_df, gdf[gdf_cols], on='gene_id', how='outer')
  
  return count_df

    
def global_input_factor(g_i, cvg_profiles, input_cvg_profiles, verbose=True, plot_coverage=True, figsize=None, ax=None):
    
    # id region 
    chrom = g_i["id"]
    
    ## ---------------
    ## Estimate lambda: globally
    ## ----------------
    
    # From: sum(chip_cvg) - l * sum(input_cvg) = 0
    # => l = sum(chip_cvg) / sum(input_cvg)
        
    # (sumarize) pool coverages (globally)
    total_sample = sum(cvg_profiles[chrom])
    total_input = sum(input_cvg_profiles[chrom])
    
    global_lambda = total_sample / total_input
    if (verbose):
        print("Lambda: {} \n".format(global_lambda))
    
    # get INPUT subtracted profile
    normed_cvg_profile = cvg_profiles[chrom] - global_lambda * input_cvg_profiles[chrom]
    
    ## -------------
    ## Visualization: Coverage after INPUT substraction
    ## -------------
    
    if (plot_coverage):
      
        # Create our figure - handle through `plot_coverage_interval` function
        # if isinstance(ax, type(None)):
        #     fig, ax = plt.subplots(figsize=figsize)
            
        cov.plot_coverage_interval(g_i, cvg = normed_cvg_profile, l=global_lambda, figsize=figsize, ax=ax)
    
    return global_lambda, normed_cvg_profile
    
    
def from_distribution_input_factor(g_i, cvg_profiles, input_cvg_profiles, stat_param='mean', verbose=True, plot_coverage=True, plot_distributions=True, figsize=None, ax=None):
    
    # id region 
    g_id = g_i["id"]
    
    ## ---------------
    ## Estimate lambda: From Distribution of local lambdas
    ## ----------------
    
    # From: chip_cvg_i - l_i * input_cvg_i = 0
    # => l_i = chip_cvg_i / input_cvg_i
    
    # compute local lambdas
    #local_lambdas = cvg_profiles[g_id] / input_cvg_profiles[g_id]
    lambdas_df = pd.DataFrame(
        data = {
            'chip_cvg':cvg_profiles[g_id],
            'input_cvg':input_cvg_profiles[g_id]
            }
        )
    lambdas_df['lambda'] = lambdas_df['chip_cvg'] / lambdas_df['input_cvg']
    lambdas_df['position'] = g_i['start'] + lambdas_df.index
    
    if (plot_distributions):
        
        lambdas_plot = lambdas_df[['position', 'chip_cvg', 'input_cvg']].melt(
            id_vars=['position'],
            var_name = 'seq_type',
            value_name = 'cvg'
            )
            
        dist_plot = sns.displot(
            data=lambdas_plot, 
            x="cvg", 
            #hue="seq_type", 
            col="seq_type",
            #col_wrap=3,
            stat="probability"
            )
    
    region_length = lambdas_df.shape[0]
    # remove 'positions' causing problems - division zero / zero
    lambdas_df = lambdas_df[~lambdas_df.isnull().any(axis=1)]
    nan_removed_positions = region_length - lambdas_df.shape[0]

    region_length = lambdas_df.shape[0]
    # remove 'positions' causing problems - after division x / zero = inf
    lambdas_df = lambdas_df[~lambdas_df.isin([np.inf]).any(axis=1)]
    inf_removed_positions = region_length - lambdas_df.shape[0]
    
    if (plot_distributions):
    
        dist_plot = sns.displot(
            data=lambdas_df[lambdas_df['lambda'] < 1],
            x="lambda",
            #hue="method",
            #col="mutant_id",
            #col_wrap=3,
            stat="probability"
            )

    dist_lambda = lambdas_df[['lambda']].describe()
    dist_lambda = dist_lambda.loc[stat_param].values[0]
    
    if (verbose):
        print('Removed `positions` causing problems:')
        print(' - Division 0/0 (Undetermined - NaN):', nan_removed_positions)
        print(' - Division x/0 (Infinite - inf):', inf_removed_positions)
        print("Use {} as summary statistic from lambda distibution: {}\n".format(stat_param, dist_lambda))

    # get INPUT subtracted profile
    normed_cvg_profile = cvg_profiles[g_id] - dist_lambda * input_cvg_profiles[g_id]
    
    ## -------------
    ## Visualization: Coverage Profile & Regression
    ## -------------
    
    if (plot_coverage):
        
        # Create our figure - handle through `plot_coverage_interval` function
        # if isinstance(ax, type(None)):
        #     fig, ax = plt.subplots(figsize=figsize)
            
        cov.plot_coverage_interval(g_i, cvg = normed_cvg_profile, l=dist_lambda, figsize=figsize, ax=ax)
        
    return dist_lambda, normed_cvg_profile

  
def linear_regression_input_factor(g_i, cvg_profiles, input_cvg_profiles, verbose=True, plot_coverage=True, figsize=None, ax=None):
    
    # id region 
    g_id = g_i["id"]
    
    ## ---------------
    ## Estimate lambda: Linear Regression
    ## ----------------
    
    # fix intercept - b = 0
    regressor = LinearRegression(fit_intercept=False)
    # TODO: check coverage have same length
    regressor.fit(input_cvg_profiles[g_id].reshape(-1,1), cvg_profiles[g_id].reshape(-1,1))
    
    a = regressor.coef_[0][0] # retrieving the slope
    b = regressor.intercept_ # retrieve the intercept
    
    if (verbose):
        print("Slope: {}, Intercept: {} \n".format(a, b))
    lregression_lambda = a

    # get INPUT subtracted profile
    normed_cvg_profile = cvg_profiles[g_id] - lregression_lambda * input_cvg_profiles[g_id]
    
    ## -------------
    ## Visualization: Coverage Profile & Regression
    ## -------------
    
    if (plot_coverage):
        
        # Create our figure
        if isinstance(ax, type(None)):
            fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize=figsize)
    
        ## -------------------------
        ## Fig 1.a. Coverage Profile - INPUT Substraction
        ## -------------------------
        
        cov.plot_coverage_interval(g_i, cvg = normed_cvg_profile, l=lregression_lambda, figsize=figsize, ax=ax[0])
    
        ## ---------------------
        ## Fig 1.b. Scatter Plot - Regression
        ## ---------------------
    
        # add 'main plot': regression
        ax[1].scatter(input_cvg_profiles[g_id], cvg_profiles[g_id], color='gray', s=1, alpha=0.5, edgecolors='face')
        ax[1].plot(input_cvg_profiles[g_id], a * input_cvg_profiles[g_id] + b, color='red', linewidth=0.8)
        
        # add title and axis names
        ax[1].set_xlabel("INPUT coverage")
        ax[1].set_ylabel("ChIP coverage")
        
        #ax[1].set_aspect(aspect=1)
        
    return lregression_lambda, normed_cvg_profile


if __name__ == "__main__":
    
    ## -----------------
    ## Script Parameters
    ## -----------------

    params = vars(get_args())

    ## initialize args/params passed to script
    
    ## (Optional): INPUT subtraction - from BAM
    param_in_bam = params['in_bam'] # path to alignment file (.bam or .sam)
    param_in_input = params['in_chip_input'] # path to `ChIP INPUT` alignment file (.bam or .sam)
    
    # ## (Optional): INPUT subtraction - from counts
    # param_in_count = params['in_count'] # path to `ChIP` count file 
    # param_in_input_count = params['in_chip_input_count'] # path to `ChIP INPUT` count file
    
    param_in_gdf = params['in_gdf'] # gene information table (gdf), containing features and corresponding attributes: gene_lengths, etc.
    
    param_verbose = params['verbose']
    param_out_dir = params['out_dir']  # save results in (default: current working directory).'
    
    ## ----------------
    ## Count Parameters
    ## ----------------
    
    # Feature/s to be parsed from GFF and used for counting - currently, only used for naming the 'output_file's
    param_feature_type = list(set(params['feature_type']))
    
    ## Potentially allow multimapped-reads - with less than max_nh (16) hits
    ## => We only allow this if the read maps to a repeat feature.
    param_max_nh = params['max_nh']

    ## How to deal with overlap of READ and FEATURE?
    # ("union", "intersection-strict", "intersection-nonempty")
    param_count_mode = params['count_mode']
    #
    ## How to deal with reads multiple-overlapping FEATURES? (according to count_mode)
    # ("none", "all", "fractional")
    param_ambiguous_assignment_mode = params['ambiguous_assignment_mode']

    ## How to deal with multi-mapped reads? including secondary alignments
     # ("none", "all", "fractional", "ignore_secondary")
    param_multimapped_mode = params['multimapped_mode']

    out_dir = '.'
    if param_out_dir:
        out_dir = param_out_dir
        
    ## ------
    ## Config
    ## ------
    
    #fragment_size = 500
    #fragment_size = 200
    fragment_size = None
    
    # - Define `ChIP` sample
    chip_sample_file = os.path.basename(param_in_bam)
    chip_sample_id = chip_sample_file.split(".")[0]
    
    # and corresponding `INPUT` sample
    input_sample_file =  os.path.basename(param_in_input)
    input_sample_id =  input_sample_file.split(".")[0]
    
    if 'H3K9me2' in chip_sample_id:
        input_region = 'mrna'
    else:
        input_region = 'centromeric'
    
    # --------------------------
    # Workflow:
    # - Load `sample` and corresponding `INPUT` BAM files.
    # - Create coverage for region of interest.
    # - Using fragment_size: None
    #     - https://htseq.readthedocs.io/en/master/tss.html#using-the-full-coverage
    # --------------------------
    
    # -----------------------
    # A. INPUT Normalization: Centromeric Region
    # -----------------------
    
    if (input_region == 'centromeric'):
        
        # - Create coverage for region of interest - Centrometric regions:
        #     - `chr I   - 3,765,776 - 3,776,547`
        #     - `chr II  - 1,619,269 - 1,629,017`
        #     - `chr III - 1,093,002 - 1,105,987`
        
        # - Define characteristics for `coverage interval` - Centromeric regions:
        #   => between consecutive dg/dh Heterochromatin peaks in WT 65_H3k9me2 sample.
        # instead we can follow:
        # https://htseq.readthedocs.io/en/master/tss.html#using-indexed-bam-files
        i_1 = {'id':'centromer_I', 'chrom':'I', 'start': 3765776, 'end': 3776547}
        i_2 = {'id':'centromer_II', 'chrom':'II', 'start': 1619269, 'end': 1629017}
        i_3 = {'id':'centromer_III', 'chrom':'III', 'start': 1093002, 'end': 1105987}
        genomic_intervals = [i_1, i_2, i_3]
        
    # -----------------------
    # B. INPUT Normalization: Protein Coding Genes (for H3K9me2)
    # -----------------------
    
    else:
        
        # - Create coverage for region of interest - long mRNA genes:
        #     - SPAC23G3.02c (sib1): `I   - 854453   - 869704` 
        #     - SPCC737.08 (mdn1):   `III - 1900472  - 1914833`
        #     - SPAC1093.06c (dhc1): `I   - 4622166  - 4635111`
 
        ## Import Gene Data Table (gdf)
        gdf = pd.read_csv(param_in_gdf, sep='\t')
        
        # - Define characteristics for `coverage interval` - long mRNA genes:
        #   => Longest **protein coding regions** (mRNA genes) which contain no heterochromatin and are therefore noise.
        mrna_gdf = gdf[(gdf['bio_type'] == 'mRNA') & ~(gdf['category'] == 'repeat')].sort_values('gene_length', ascending=False)
        
        genomic_intervals = []
        #for gg_row in mrna_gdf.head(3).itertuples():
        for gg_row in mrna_gdf.itertuples():
            genomic_intervals.append({'id':gg_row.gene_id, 'chrom':gg_row.seqid, 'start':gg_row.start, 'end':gg_row.end})
            
    #genomic_intervals
    #len(genomic_intervals)
    
    # -------------------------------------------------------------
    # --  1. Use .bam files to get coverage in genomic regions   --
    # -------------------------------------------------------------

    coverage_plot_dir = os.path.join(out_dir, 'coverage_plots') 
    if not os.path.exists(coverage_plot_dir):
        os.makedirs(coverage_plot_dir)

    ## ---------------
    ## A. ChIP Sample:
    ## ---------------
    
    # - Load `BAM` files: instantiate BAM_Reader Object
    #bam_file = os.path.join(data_dir, data_batch, "bam", dataset_id, dataset_id + ".Aligned.sortedByCoord.out.bam")
    bam_file = param_in_bam
    #print(bam_file)
    #bam = HTSeq.BAM_Reader(bam_file)
    
    # - Instantiate a `GenomicArray` object for the **ChIP sample** coverage (`cvg`)
    #cvg = cov.coverage_genomic_intervals(bam_file, chrom, start, end, count_type="int", fragment_size=fragment_size)
    cvg = cov.coverage_genomic_intervals(bam_file, genomic_intervals, count_type="frac", fragment_size=fragment_size)
    
    fig, ax = plt.subplots(nrows = 3, ncols = 1)

    # - Visualize coverage (`cvg`) in the genomic intervals of interest
    # => return `cvg_profiles` as np.array()
    cvg_profiles = {}
    for ii, gg_i in enumerate(genomic_intervals):
        # only plot first 3 regions if there are too many
        if (ii < 3):
            cvg_profiles[gg_i["id"]] = cov.plot_coverage_interval(gg_i, cvg, ax=ax[ii])
        else:
            cvg_profiles[gg_i["id"]] = cov.plot_coverage_interval(gg_i, cvg, plot_coverage=False)
            
    #cvg_profiles
    fig.tight_layout(pad=1.0)
    plt.savefig(os.path.join(coverage_plot_dir, 'chip_sample_coverage_' + input_region + '_regions.pdf'), figsize=(10, 20))

    ## ----------------
    ## B. INPUT Sample:
    ## ----------------
    
    # - Load `BAM` files: instantiate BAM_Reader Object
    #input_bam_file = os.path.join(data_dir, data_input_batch, "bam", input_dataset_id, input_dataset_id +  ".Aligned.sortedByCoord.out.bam")
    input_bam_file = param_in_input
    #print(input_bam_file)
    #input_bam = HTSeq.BAM_Reader(input_bam_file)
    
    # - Instantiat|e a `GenomicArray` object for the **INPUT sample** coverage (`input_cvg`)
    #input_cvg = cov.coverage_genomic_intervals(input_bam_file, chrom, start, end, fragment_size=fragment_size)
    input_cvg = cov.coverage_genomic_intervals(input_bam_file, genomic_intervals, count_type="frac", fragment_size=fragment_size)
    
    fig, ax = plt.subplots(nrows = 3, ncols = 1)
    
    # - Visualize coverage (`cvg`) in the genomic intervals of interest
    # => return `input_cvg_profiles` as np.array()
    input_cvg_profiles = {}
    for ii, gg_i in enumerate(genomic_intervals):
        # only plot first 3 regions if there are too many
        if (ii < 3):
            input_cvg_profiles[gg_i["id"]] = cov.plot_coverage_interval(gg_i, input_cvg, ax=ax[ii])
        else:
            input_cvg_profiles[gg_i["id"]] = cov.plot_coverage_interval(gg_i, input_cvg, plot_coverage=False)
            
    #input_cvg_profiles
    fig.tight_layout(pad=1.0)
    plt.savefig(os.path.join(coverage_plot_dir, 'input_sample_coverage_' + input_region + '_regions.pdf'), figsize=(10, 20))

    # -------------------------------------------------------------
    # ------     2. Estimate INPUT substraction factor      -------
    # -------------------------------------------------------------
    
    lambda_estimate_plot_dir = os.path.join(out_dir, 'lambda_plots') 
    if not os.path.exists(lambda_estimate_plot_dir):
        os.makedirs(lambda_estimate_plot_dir)
    
    # ------------------------------
    # A. Estimate INPUT Norm Factor: Globally
    # ------------------------------
    
    global_lambdas = {}
    normed_cvg_profiles = {}
    
    fig, ax = plt.subplots(nrows = 3, ncols = 1)
    
    # - Loop over `genomic_intervals`
    for ii, g_i in enumerate(genomic_intervals):
        # only plot first 3 regions if there are too many
        if (ii < 3):
            global_lambdas[g_i["id"]], normed_cvg_profiles[g_i["id"]] = global_input_factor(g_i, cvg_profiles, input_cvg_profiles, ax=ax[ii])
        else:
            global_lambdas[g_i["id"]], normed_cvg_profiles[g_i["id"]] = global_input_factor(g_i, cvg_profiles, input_cvg_profiles, verbose=False, plot_coverage=False)
            
    fig.tight_layout(pad=1.0)
    plt.savefig(os.path.join(lambda_estimate_plot_dir, 'lambda_global_estimate_input_subtraction.pdf'), figsize=(10, 20))

    # ------------------------------
    # B. Estimate INPUT Norm Factor: by Linear Regression (Optimization)
    # ------------------------------
    
    # - with `polyfit` can't fix the intercept = b = 0
    #a, b = np.polyfit(window_input_cvg, window_cvg, 1, )
    #b
    
    # - with `sklearn.linear_model.LinearRegression`
    lregression_lambdas = {}
    lregression_normed_cvg_profiles = {}
    
    widths = [20, 10]
    heights = [10, 10, 10]
    ## https://matplotlib.org/3.1.0/tutorials/intermediate/gridspec.html
    gs_kw = dict(width_ratios=widths, height_ratios=heights)

    fig, ax = plt.subplots(nrows = 3, ncols = 2, gridspec_kw=gs_kw)
    
    # - Loop over `genomic_intervals`
    for ii, g_i in enumerate(genomic_intervals):
        # only plot first 3 regions if there are too many
        if (ii < 3):
            lregression_lambdas[g_i["id"]], lregression_normed_cvg_profiles[g_i["id"]] = linear_regression_input_factor(g_i, cvg_profiles, input_cvg_profiles, ax=ax[ii])
        else:
            lregression_lambdas[g_i["id"]], lregression_normed_cvg_profiles[g_i["id"]] = linear_regression_input_factor(g_i, cvg_profiles, input_cvg_profiles, verbose=False, plot_coverage=False)
            
    fig.tight_layout(pad=1.0)
    plt.savefig(os.path.join(lambda_estimate_plot_dir, 'lambda_linear_regression_estimate_input_subtraction.pdf'), figsize=(10, 30))
    
    # ------------------------------
    # C. Estimate INPUT Norm Factor: from Counts
    # ------------------------------
    
    # Compute counts for ChIP
    chip_count_df = raw_counts_centromer(
      bam_file = bam_file,
      in_gdf = param_in_gdf, 
      features_list = genomic_intervals,
      stranded = False, 
      count_mode = param_count_mode,
      ambiguous_assignment_mode = param_ambiguous_assignment_mode,
      multimapped_mode = param_multimapped_mode,
      max_nh = param_max_nh
      )
    
    # Compute counts for INPUT
    ip_count_df = raw_counts_centromer(
      bam_file = input_bam_file, 
      in_gdf = param_in_gdf, 
      features_list = genomic_intervals,
      stranded = False, 
      count_mode = param_count_mode,
      ambiguous_assignment_mode = param_ambiguous_assignment_mode,
      multimapped_mode = param_multimapped_mode,
      max_nh = param_max_nh
      )
    
    count_df = chip_count_df.merge(
        ip_count_df.rename(columns={"count": "count_ip"})[['gene_id', 'count_ip']],
        how = 'inner',
        on='gene_id'
        )
        
    count_df['input_factor'] = count_df['count'] / count_df['count_ip']
    count_df['method'] = 'count_lambda'
    
    # -------------------------
    # Agreement between lambdas
    # -------------------------
    
    # - Good agreement between both approaches and across regions:
    #lregression_lambdas
    #global_lambdas
    
    # - Convert results to Data Frame
    combined_lambdas = {key:[lregression_lambdas[key], global_lambdas[key]] for key in lregression_lambdas}
    combined_lambdas['method'] = ['linear_regression', 'global_lambda']
    
    lambda_df = pd.DataFrame.from_dict(combined_lambdas)
    lambda_df = lambda_df.melt(id_vars='method', var_name='chrom', value_name='input_factor')
    #lambda_df = lambda_df.melt(id_vars='method', var_name='region_id', value_name='input_factor')
    
    # add 'count_lambda'
    count_df = count_df.rename(columns={"gene_id": "chrom"})[['method', 'chrom', 'input_factor']]
    lambda_df = pd.concat([lambda_df, count_df]).reset_index(drop=True)

    # -Add ChIP id & INPUT id:
    lambda_df['chip_id'] = chip_sample_id
    lambda_df['input_id'] = input_sample_id

    #lambda_df['mean'] = lambda_df.mean(axis=1)
    #print(" Scaling factor (λ) used in the analysis: {:.3f}".format(mean_lambda))
    
    # - Store as .csv file
    lambda_input_factors_file = os.path.join(out_dir, "INPUT_factors.csv")
    lambda_df.to_csv(lambda_input_factors_file, sep='\t', index=None)
    
    if (input_region == 'mrna'):
        
        import seaborn as sns
        print(sns.__version__)
        plt.style.use('ggplot')
        
        # 95th Percentile
        def q95(x):
            return x.quantile(0.95)
        # 98th Percentile
        def q98(x):
            return x.quantile(0.98)
        
        # A. Summarize results for **lambda coefficients** (across genes) ------------
        lambda_summary = lambda_df.groupby('method').agg(
            # Get max of the 'input_factor' column for each group
            max_if =('input_factor', max),
            # Get min of the 'input_factor' column for each group
            min_if=('input_factor', min),
            # Get mean of the 'input_factor' column for each group
            mean_if=('input_factor', 'mean'),
            # Get median of the 'input_factor' column for each group
            median_if=('input_factor', 'median'),
            # Get standard deviation of the 'input_factor' column for each group
            std_if=('input_factor', 'std'),
            # Get 95th quantile of the 'input_factor' column for each group
            quantile_95_if=('input_factor', q95),
            # Get 98th quantile of the 'input_factor' column for each group
            quantile_98_if=('input_factor', q98),
            # Apply a lambda to date column
            #num_days=("date", lambda x: (max(x) - min(x)).days)    
            )
            
        # - Store as .csv file
        summary_lambda_input_factors_file = os.path.join(out_dir, "summary_INPUT_factors.csv")
        lambda_summary.to_csv(summary_lambda_input_factors_file, sep='\t', index=None)
        
        # Visualize lambda factor **distributions**

        lambda_max = 5
        #lambda_max = max(lambda_summary['quantile_98_if'])
        
        lambda_plot = lambda_df[lambda_df['input_factor'] < lambda_max]
        #lambda_plot = lambda_df[lambda_df['input_factor'] > lambda_max]

        sns.displot(lambda_plot, x="input_factor", hue="method", stat="probability")
        plt.savefig(os.path.join(lambda_estimate_plot_dir, 'distribution_lambda_estimate_input_subtraction.pdf'), figsize=(10, 15))
