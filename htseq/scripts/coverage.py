#!/usr/bin/env python

import itertools
import time
import os
import sys
import argparse

import numpy as np

import HTSeq

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches

import collections
import pickle

import pandas as pd
import copy


## ----------
## Parameters
## ----------

# CIGAR match characters (including alignment match, sequence match, and sequence mismatch)
com = ('M', '=', 'X')

########################################################
##  Divide Heterochromatic genes in 4 groups/regions  ##
########################################################

# ------------------------------
# A. Centromeric repeats `deg1`: dg, dh counts
# ------------------------------

#deg1 = ['dh1', 'dg1'] ## unstranded
#deg1 = ['dh1+',  'dh1-', 'dg1+', 'dg1-']
deg1 = ['dg1a', 'dg1b', 'dh1a', 'dh1b'] ## gff_v2
deg1_names = {'dg1a': 'dg1+', 'dg1b': 'dg1-',
              'dh1a': 'dh1+', 'dh1b': 'dh1-'}

# -----------------------------
# B. Subtelomeric Genes `deg2`: tlh, SPAC212.10 counts
# -----------------------------

#deg2 = ['tlh1', 'SPAC212.10'] ## unstranded
#deg2 = ['tlh1+', 'tlh1-', 'SPAC212.10']
deg2 = ['SPAC212.11', 'SPAC212.11b', 'SPAC212.10', 'SPAC212.10b'] ## gff_v2
deg2_names = {'SPAC212.11': 'tlh1-', 'SPAC212.11b': 'tlh1+',
              'SPAC212.10': 'SPAC212.10-', 'SPAC212.10b': 'SPAC212.10+'}

# ------------------------------------
# C.  Mating type region (MTR) `deg3`: MAT, other regions, etc... counts 
# ------------------------------------

deg3 = ['MAT2', 'MAT3', 'MAT1']
deg3_names = {'MAT1': 'MAT1', 'MAT2': 'MAT2', 'MAT3': 'MAT3'}

# ------------------------------------------
# D. rest Heterochromatic + mat gene counts: `deg4`
# ------------------------------------------

non_degraded = ['SPAC212.09c', 'SPNCRNA.70', 'SPAC212.08c', 'SPAC212.07c', 'SPAC212.12', 'SPAC212.06c',
                'SPAC212.04c', 'SPAC212.03', 'SPAC212.02', 'SPAC212.01c', 'SPAC977.01', 'SPAC977.18',
                'SPAC977.02', 'SPAC977.03', 'SPAC977.04', 'SPAC212.05c']
                #'MAT2', 'MAT3', 'MAT1']
    
# ---------------
# E. Other genes: `other_genes`
# ---------------

other_genes = ['SPAC212.09c', 'SPAC977.01', 'SPAC977.10', 'SPAP7G5.06']
other_gene_nanes = {'SPAC212.09c':'SPAC212.09c', 'SPAC977.01':'ftm1', 'SPAC977.10':'nhe1', 
                    'SPAP7G5.06':'per1'}

htc_genes = deg1 + deg2
#htc_genes = deg1 + deg2 + deg3 + non_degraded

visualize_genes = deg1 + deg2 + deg3 + other_genes
visualize_gene_names = {**deg1_names, **deg2_names, **deg3_names, **other_gene_nanes}


## -------------------
## Auxiliary Functions
## -------------------

def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--in_bam', '-f',
        type=str, default=None,
        help='Input bam file for which coverage is to be computed.')
    
    ## (Optional): INPUT subtraction
    parser.add_argument(
        '--in_chip_input', '-i',
        type=str, default=None,
        help='Input ChIP INPUT bam file which coverage scaled by a norm factor needs to be subtracted from coverage.')
   
    ## (Optional): INPUT subtraction
    parser.add_argument(
        '--norm_factor', '-n',
        type=str, default=None,
        help='Norm factor used to scale the ChIP INPUT coverage before subtraction from coverage.')
   
    parser.add_argument(
        '--in_gtf', '-g',
        type=str,
        help='Reference Annotation File (GTF/GFF).', required=True)
    
    parser.add_argument(
        '--verbose', '-v',
        type=bool, default=True,
        help='Verbose (default: True).')
    
    parser.add_argument(
        '--out_dir', '-o',
        type=str, default='.',
        help='Save results in (default: current working directory).')
    
    parser.add_argument(
        '--seq_type', '-s',
        type=str, default='ChIP',
        help='Sequencing technology, for now `ChIP` or `RNA` (e.g. pA-RNA, RIP, total-RNA ...), this has to do with the `strandedness` of the data.')

    args = parser.parse_args()

    return args

def invert_strand(iv):
    
    iv2 = iv.copy()
    
    if iv2.strand == "+":
        iv2.strand = "-"
        
    elif iv2.strand == "-":
        iv2.strand = "+"
        
    else:
        raise ValueError("Illegal strand")
        
    return iv2


## ----------------------------
## Calculating coverage vectors
## ----------------------------

# - Compute coverage vector for a specific `genomic interval`:
#    * Simplified version only used for ChIP!
#    => e.g. might have issues with spliced reads
def coverage_vectors(bam_file, seq_type='ChIP', typecode='i'):
    """
    By `coverage vector`, we mean a vector (one-dimensional array) of the length of a chromosome, where
    each element counts how many reads cover the corresponding base pair in their alignment.

    A `GenomicArray` can conveniently bundle the coverage vectors for all the chromosomes in a genome.

    Hence, we start by defining a `GenomicArray` for the coverage(`cvg`):
    cvg = HTSeq.GenomicArray("auto", stranded=True, typecode="i")

    Instead of listing all chromosomes, we instruct the `GenomicArray` to add chromosome vectors as needed,
    by specifiyng `"auto"`. As we set `stranded=True`, there are now two chromosome vectors for each chromosome,
    all holding integer values (`typecode="i"`).

    They all have an "infinite" length as we did not specify the actual lengths of the chromosomes.

    To fill the `coverage vectors`, we now simply iterate through all the reads and add, either:
        - value: 1
        - value: 1/NH
    depending if we are interest in the `integer` or `fractional` coverage at the interval to which each read
    was aligned to.

    :param bam_file:
    :param seq_type:
    :return:
    """

    ## initialize counter
    counts = collections.Counter()

    # select bam file: instantiate BAM_Reader Object
    bam = HTSeq.BAM_Reader(bam_file)

    if seq_type == 'ChIP':
        strandedness = False

    else:
        strandedness = True

    ## instantiate a `GenomicArray` object for the coverage(`cvg`):
    int_cvg = HTSeq.GenomicArray("auto", stranded=strandedness, typecode=typecode)
    fract_cvg = HTSeq.GenomicArray("auto", stranded=strandedness, typecode='d')

    start_time = time.time()

    ## -------
    ## Options
    ## -------

    i = 0

    ## simply iterate through all the reads
    #for aln in itertools.islice(bam_file, 10):
    for aln in bam:

        if i > 0 and i % 100000 == 0:
            sys.stderr.write("{} alignment records processed. {} s\n".format(i,  time.time() - start_time))
            sys.stderr.flush()
        i += 1

        ## ----------------------
        ## Inspect read alignment
        ## ----------------------

        ## _mapped or _unmapped (our BAM files only contain _mapped)
        counts["_total"] += 1

        if not aln.aligned:

            counts["_unmapped"] += 1
            ## skips to next iteration
            continue

        else:

            ## invert strand - due to sequencing is strand is reversed!
            if strandedness:
                iv_seq = (invert_strand(co.ref_iv) for co in aln.cigar if (co.type in com and co.size > 0))
            else:
                iv_seq = (co.ref_iv for co in aln.cigar if (co.type in com and co.size > 0))

            for iv in iv_seq:

                ## add the value 1 at the interval to which each read was aligned to
                int_cvg[ iv ] += 1

                ## add the corresponding (fractional) value at the interval to which each read was aligned to
                fract_cvg[ iv ] += 1 / aln.optional_field('NH')


    print('Elapsed Time (Counting reads):', time.time() - start_time)

    return int_cvg, fract_cvg, counts
    

# - Compute coverage vector for a specific `set of genomic intervals`:
#    * Simplified version only used for ChIP!
#    => e.g. might have issues with spliced reads
def coverage_genomic_intervals(bam_file, genomic_intervals, count_type="int", stranded=False, fragment_size=None):

    print("Loading BAM file ...\n{}".format(bam_file))
    
    # Load `BAM` files: instantiate BAM_Reader Object
    bam = HTSeq.BAM_Reader(bam_file)
    #bam = pysam.AlignmentFile(bam_file, 'rb')
    
    print("Done.\n")

    # Instantiate a `GenomicArray` object for the coverage (`cvg`)
    if count_type == "int":
        #cvg = HTSeq.GenomicArray("auto", stranded=stranded, typecode='i') 
        cvg = HTSeq.GenomicArray("auto", stranded=stranded, typecode='d') ## needed for INPUT subtraction

    elif count_type == "frac":
        cvg = HTSeq.GenomicArray("auto", stranded=stranded, typecode='d')

    else:
        raise ValueError("Illegal strand!")
    
    padding = 0 if isinstance(fragment_size, type(None)) else fragment_size
    print('Add padding:', padding, '\n')

    # function expects - list of `genomic intervals` 
    # => dict with "chrom", "start", "end" key-value pairs
    assert isinstance(genomic_intervals, type([]))
    regions = [ "{}:{}-{}".format(g_i["chrom"], 
                                  g_i["start"] - padding, # left-boundary condition
                                  g_i["end"] + padding) if g_i["start"] >= padding else
                "{}:{}-{}".format(g_i["chrom"], 
                                  1,
                                  g_i["end"] + padding) for g_i in genomic_intervals]
                            
    print('Genomic Regions:')
    n_max_regions = 10
    # loop over - genomic intervals
    for ii, region in enumerate(regions):
        
        g_i = genomic_intervals[ii]
        
        # only print first regions
        if ii + 1 < n_max_regions:
            print("{}- {} ({:.2f} kb)".format('\t', region, (g_i["end"] - g_i["start"]) / 10**3 ))
        elif ii + 1 == n_max_regions:
            print("{}- {}".format('\t', "(...)"))
        
        # loop over - read alignments
        for aln in bam.fetch(region = region):
            
            # TODO: possible improvement - change read interval extending to fragment size
            if not isinstance(fragment_size, type(None)):
                
                if aln.iv.start >= fragment_size + 1:
                    aln.iv.length = fragment_size
                else:
                    aln.iv.length = aln.iv.length + aln.iv.start - 1
                    
            # Define count weight schemes: "int" or "frac"
            if count_type == "int":
                # add the value 1 at the interval to which each read was aligned to
                #cvg[ aln.iv ] += 1 
                count_weight = 1
                
            else:
                # add the corresponding (fractional) value at the interval to which each read was aligned to
                #cvg[ aln.iv ] += (1 / aln.optional_field('NH'))
                count_weight = (1 / aln.optional_field('NH'))
            
            # Spliced - CIGAR & Strandedness of the Data
            if stranded:
                # invert strand - due to sequencing is strand is reversed!
                iv_seq = (invert_strand(co.ref_iv) for co in aln.cigar if (co.type in com and co.size > 0))
            else:
                iv_seq = (co.ref_iv for co in aln.cigar if (co.type in com and co.size > 0))
            
            # add read 'count_weight' to coverage vector
            for iv in iv_seq:
                # add the corresponding `count_weight` value at the interval to which each read was aligned to
                cvg[ iv ] += count_weight

    return cvg


## -----------------------------
## Visualizing: Genomic Features
## -----------------------------

## - Get **genomic features** in regions of interest
def get_genomic_features_dfs(genomic_intervals, in_gtf, feature_types=['gene']):
    
    # Select gtf file: instantiate GFF_Reader Object
    gtf = HTSeq.GFF_Reader(param_in_gtf)
    
    # Init dictionary containing info about features - to later create Data Frames
    init_feature_dict = {'region_id':[], 'feature_id':[], 'feature_name':[], 'type':[], 'chrom':[], 'start':[], 'end':[], 'strand':[]}
    genomic_features = {gg['region_id']: copy.deepcopy(init_feature_dict) for gg in genomic_intervals}

    ## Loop over features in gtf:
    for feature in gtf:
        
        if feature.type in feature_types:
                
            ## Loop over genomic regions of interest: check if `feature` falls within the genomic region
            for gi in genomic_intervals:
                
                if (feature.iv.chrom == gi['chrom'] and feature.iv.start + 1 >= gi['start'] and feature.iv.end + 1 <=  gi['end']):
                    
                    genomic_features[gi['region_id']]['region_id'].append(gi['region_id'])
                    genomic_features[gi['region_id']]['feature_id'].append(feature.attr["ID"])
                    try:
                        genomic_features[gi['region_id']]['feature_name'].append(feature.attr["Name"])
                    except:
                        genomic_features[gi['region_id']]['feature_name'].append(np.nan)
                    genomic_features[gi['region_id']]['type'].append(feature.type)
                    genomic_features[gi['region_id']]['chrom'].append(feature.iv.chrom)
                    genomic_features[gi['region_id']]['start'].append(feature.iv.start)
                    genomic_features[gi['region_id']]['end'].append(feature.iv.end)
                    genomic_features[gi['region_id']]['strand'].append(feature.iv.strand)
                    
    ## Convert to Data Frames
    genomic_features_dfs = {}
    
    for region_id, gg_features in genomic_features.items():
        genomic_features_dfs[region_id] = pd.DataFrame.from_dict(gg_features)
        
    return genomic_features_dfs

def plot_non_overlapping_features(genomic_features_df, ax = None, annotate_features = True, genomic_range = None):
    
    # Create our figure
    if isinstance(ax, type(None)):
        fig, ax = plt.subplots(nrows = 1, ncols=1)
    
    # Create list for all the rectangle patches (genomic features)
    rectangle_patches = []

    ## Tracks last rectangle in the Y-axis, in order to avoid overlap between features:
    y_tracks = [0, 0, 0]

    genomic_features_df = genomic_features_df.sort_values(['chrom', 'start', 'end'])

    # Loop over features in genomic_region; create rectangle defined by each feature
    for row in genomic_features_df.itertuples():
        
        ## Left bottom x-position
        x = row.start

        track_idx = 0
        padding = 10
        ## Find y-track to avoid overlaps
        while (y_tracks[track_idx] != 0) and (y_tracks[track_idx] + padding >= row.start):
            track_idx += 1
        y_tracks[track_idx] = row.end
        y = track_idx + 0.2

        ## Rectangle - width
        width = row.end - row.start + 1
        height = 0.8
        
        ## Color - stranded
        color = 'r' if row.strand == '+' else 'b'

        ## Use patch collection
        #rect = mpatches.Rectangle(xy = (x, y), width = width, height = height)
        #rectangle_patches.append(rect)

        # Add individual patches to the Axes
        rect = mpatches.Rectangle(xy = (x, y),
                                  width = width, height = height, 
                                  facecolor=color, alpha=0.3, edgecolor='k', linewidth=1)
        ax.add_patch(rect)

        # Add feature_id's with text:
        if annotate_features:
            
            if isinstance(annotate_features, dict):
                
                label = annotate_features[row.feature_id] if row.feature_id in annotate_features else ""

            else:
                
                trim_name = 8
                label = row.feature_id if pd.isnull(row.feature_name) else row.feature_name
                #label = label if len(label) <= trim_name else label[:trim_name]
                label = label if len(label) <= trim_name else ""

            ax.text(x + width/2, y + height/2, label, fontsize=4, ha="center", va="center")

    # Create patch collection with specified colour/alpha
    #pc = PatchCollection(rectangle_patches, facecolor='r', alpha=0.3, edgecolor='k', linewidth=1)
    
    # Add collection to axes
    #ax.add_collection(pc)
    
    ## Set x-Axis limit
    if isinstance(genomic_range, type(None)):
        x_min = genomic_features_df['start'].min()
        x_max = genomic_features_df['end'].max()
        
    else:
        assert len(genomic_range) == 2
        x_min = genomic_range[0]
        x_max = genomic_range[1]
    
    ## Add plot limits:
    ax.set_xlim(xmin = x_min, xmax = x_max)

    #ax.axvline(x_min, linewidth = 2, alpha = 0.4, color = 'k')
    #ax.text(x_min, 0, 'x_min = {}'.format(x_min))

    #ax.axvline(x_max, linewidth = 2, alpha = 0.4, color = 'k')
    #ax.text(x_max, 0, 'x_max = {}'.format(x_max))
    
    ## Add padding to plot limits:
    #axis_padding = (x_max - x_min) * 0.05
    #ax.set_xlim(xmin = x_min - axis_padding, 
    #            xmax = x_max + axis_padding)
    
    ## Depends on the number of tracks - inverted
    ax.set_ylim(bottom = len(y_tracks), top = 0)

    ## Configure Axis
    #ax.xaxis.set_tick_params(labelsize=5)
    #ax.set_xlabel('genomic coordinates', fontsize = 5)

    #a = 20
    #ax.set_aspect((x_max - x_min) / a)
    
    #ax.yaxis.set_visible(False)
    #ax.xaxis.set_visible(False)
    ax.axis('off')
    #ax.grid(False)

    return ax


## -----------------------------
## Visualizing: Coverage Vectors
## -----------------------------

# - Visualize a genomic_interval coverage profile (optionally an input factor, l, can be passed for annotation of plot)
def plot_coverage_interval(genomic_interval, cvg, plot_coverage=True, strand = ".", l=None, color = 'k', figsize=None, ax=None, y_scale="linear", invert=False, show_title=True):
    """
    TODO: log plot does not look so good. Either use the:
        - ax.semilogy
        - ax.set_yscale(y_scale)
    """
    
    if isinstance(cvg, HTSeq.GenomicArray):
        window = HTSeq.GenomicInterval(genomic_interval["chrom"], 
                                       genomic_interval["start"], 
                                       genomic_interval["end"],
                                       strand)
        window_cvg = np.fromiter(cvg[ window ], dtype='d')
        
    else:
        window_cvg = cvg
    
    
    if plot_coverage:
        
        genomic_iv_size = genomic_interval["end"] - genomic_interval["start"] + 1
        
        font_size = 5
        
        # create our figure
        if isinstance(ax, type(None)):
            #fig, ax = plt.subplots()
            fig, ax = plt.subplots(figsize=figsize)
            font_size = font_size * figsize[0]/8
        
        # add 'main plot': coverage profile
        if y_scale == 'log':
            y_shift = 1
            ax.semilogy(np.arange(genomic_interval["start"], genomic_interval["end"]), y_shift + window_cvg, linewidth=0.8, alpha=0.8, color=color)
        else:
            y_shift = 0
            ax.plot(np.arange(genomic_interval["start"], genomic_interval["end"]), window_cvg, linewidth=0.8, alpha=0.8, color=color)
    
        # add title and axis names
        if show_title:
            ax.set_title("{}: {}-{} ({:.1f} kb)".format(
                genomic_interval["region_id"].replace('_', " ") if "region_id" in genomic_interval else "Chromosome " + genomic_interval["chrom"], 
                genomic_interval["start"],
                genomic_interval["end"],
                genomic_iv_size / (10**3)),
                fontsize = 1.5 * font_size
                )
                                                
        # add plot limits:
        #ax.set_ylim(bottom = 0)
        
        ## Configure Axis
        ax.xaxis.set_visible(True)
        ax.xaxis.set_tick_params(labelsize=font_size)
        ax.yaxis.set_tick_params(labelsize=font_size)
        
        ## Invert Y-Axis
        if invert:
            ax.invert_yaxis()
    
        ax.set_xlabel('genomic coordinates', fontsize = font_size)
        ax.set_ylabel('coverage', fontsize = font_size)
        #ax.set_yscale(y_scale)
        
        ## Horizontal line at 0
        #ax.axhline(y = y_shift, linewidth=0.5, color='r', alpha=0.8, linestyle='--')
        ax.axhline(y = y_shift, linewidth=0.5, color='k', alpha=0.8, linestyle='--')
        
        if not isinstance(l, type(None)):
            
            from matplotlib.offsetbox import AnchoredText
            at = AnchoredText("$\lambda$ = {:.2f}".format(l),
                              prop=dict(size=8), frameon=True,
                              loc='upper left',
                             )
            at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
            ax.add_artist(at)
        
    return window_cvg

def plot_coverage_genomic_regions(genomic_intervals, cvg, genomic_features_dfs, stranded = False, annotate_features=True, y_scale="linear", prefix='', out_dir='.'):
    
    ## define gridspec gs_kw
    if not stranded:
        nrows = 2
        ncols = 3
        
        widths = [15, 5, 5]
        heights = [5, 1]
        
        color = 'k'
        strand = "."
        
    else:
        nrows = 3
        ncols = 3
        
        widths = [15, 5, 5]
        heights = [5, 5, 1]
        
        color = 'r'
        strand = "+"
        
    ## https://matplotlib.org/3.1.0/tutorials/intermediate/gridspec.html
    gs_kw = dict(width_ratios=widths, height_ratios=heights)
    
    #fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize=(sum(widths), sum(heights)), sharex='col', gridspec_kw=gs_kw)
    #fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize=(sum(widths), sum(heights)), sharex='col', sharey='col', gridspec_kw=gs_kw)
    fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize=(sum(widths), sum(heights)), sharex='col', sharey='all', gridspec_kw=gs_kw)
    
    ## The problem is that the ticker is still the same for all axes. So one will need to:
    ##  - Remove the axes from the grouper
    ##  - Set a new Ticker with the respective new locator and formatter
    ## https://stackoverflow.com/questions/54915124/how-to-unset-sharex-or-sharey-from-two-axes-in-matplotlib
    for xx in ax[-1]:
        
        ## Remove sharing y-Axis with between coverage tracks and feature tracks
        xx.get_shared_y_axes().remove(xx)

        # Create and assign new ticker
        yticker = matplotlib.axis.Ticker()
        xx.yaxis.major = yticker

        # The new ticker needs new locator and formatters
        yloc = matplotlib.ticker.AutoLocator()
        yfmt = matplotlib.ticker.ScalarFormatter()

        xx.yaxis.set_major_locator(yloc)
        xx.yaxis.set_major_formatter(yfmt)

    # - Visualize coverage (`cvg`) in the genomic intervals of interest
    for ii, gg_i in enumerate(genomic_intervals):
        plot_coverage_interval(genomic_interval = gg_i, cvg = cvg, strand=strand, color=color, ax = ax[0][ii], y_scale=y_scale)
    
    # - Visualize coverage (`cvg`) in the genomic intervals of interest
    if stranded:
        for ii, gg_i in enumerate(genomic_intervals):
            #plot_coverage_interval(genomic_interval = gg_i, cvg = cvg, strand="-", color='b', ax = ax[1][ii], y_scale=y_scale, invert=True, show_title=False)
            plot_coverage_interval(genomic_interval = gg_i, cvg = cvg, strand="-", color='b', ax = ax[1][ii], y_scale=y_scale, show_title=False)

    # - Add feature annotation in the genomic intervals of interest
    for ii, gg_i in enumerate(genomic_intervals):
        
        genomic_features_df = genomic_features_dfs[gg_i['region_id']]
        genomic_range = (gg_i['start'], gg_i['end'])
        plot_non_overlapping_features(genomic_features_df, ax = ax[-1][ii], annotate_features = annotate_features, genomic_range = genomic_range)
    
    fig.tight_layout(pad=1.0)
    plt.savefig(os.path.join(out_dir, prefix + 'coverage_genomic_regions.pdf'))
    
    return ax


## -----------------------------
## Visualizing: Coverage Vectors (OLD)
## -----------------------------

# - Visualize each feature **individually**
def plot_feature_coverage(cvg, feature_name, features_dict, out_dir=''):
    
    feature = features_dict[feature_name][0]
    
    fig = plt.figure()
    
    #plt.plot(np.arange(feature.iv.start, feature.iv.end), list(cvg[ feature.iv ]))
    plt.plot(list(cvg[ feature.iv ]))

    plt.title('Coverage Vector: {}'.format(feature_name), fontsize=15)

    plt.xlabel('Position', fontsize=10)
    plt.ylabel('Coverage', fontsize=10)

    plt.xticks(fontsize=8, rotation=45)

    #plt.savefig(os.path.join('/home/pmonteagudo/workspace/RNAdeg/htseq/plots/coverage', feature_name + '_coverage.pdf'))
    plt.savefig(os.path.join(out_dir, feature_name + '_coverage.pdf'))

# - Visualize all features as **subplots**
def subplots_feature_coverage(cvg, features_dict, n_rows=2, n_cols=3, fig_size=(100, 50), prefix='', out_dir=''):
    
    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=fig_size)
    
    ii=0 ## rows
    jj=0 ## columns
    
    feature_names = list(features_dict)
    
    for row in ax:
        
        for col in row:
            
            feature_name = feature_names[ii * n_cols + jj]
            feature = features_dict[feature_name][0]
            
            #print(ii, jj, feature)
        
            col.plot(np.arange(feature.iv.start, feature.iv.end), list(cvg[ feature.iv ]))

            col.set_title('Coverage Vector: {}'.format(feature_name), fontsize=30)

            col.set_xlabel('Position', fontsize=20)
            col.set_ylabel('Coverage', fontsize=20)

            #col.set_xticks(fontsize=8, rotation=45)
            col.tick_params(axis="x", size=15, rotation=45) 

            jj += 1
            
        jj = 0
        ii =+1
        
        plt.savefig(os.path.join(out_dir, prefix + 'subplots_feature_coverage.pdf'))
        #plt.savefig(os.path.join('/home/pmonteagudo/workspace/RNAdeg/htseq/plots/coverage', feature_name + '_coverage.pdf'))


if __name__ == "__main__":
    
    ## -----------------
    ## Script Parameters
    ## -----------------

    params = vars(get_args())

    ## initialize args/params passed to script
    param_in_bam = params['in_bam'] ## path to alignment file (.bam or .sam)
    
    ## (Optional): INPUT subtraction
    param_in_input = params['in_chip_input'] ## path to ChIP INPUT alignment file (.bam or .sam)
    param_norm_factor = params['norm_factor']
    
    param_in_gtf = params['in_gtf'] ## Gene information table (gdf)
    param_seq_type = params['seq_type']

    param_verbose = params['verbose']
    param_out_dir = params['out_dir']  ## Save results in (default: current working directory).'

    out_dir = '.'
    if param_out_dir:
        out_dir = param_out_dir
        #out_dir = '/home/pmonteagudo/workspace/RNAdeg/htseq/plots/coverage'

    ## get sample_file
    sample_file = os.path.basename(param_in_bam)
    sample_id = sample_file.split(".")[0]

    out_dir = os.path.join(out_dir, sample_id)
    coverage_plot_dir = os.path.join(out_dir, "coverage_plots")
    if not os.path.exists(coverage_plot_dir):
        os.makedirs(coverage_plot_dir)

    # -------------------------------------------------------------
    # -------    Whole Genome Coverage & Visualization     --------
    # -------------------------------------------------------------

    # --------------------------
    # Workflow:
    # - Load `sample` BAM file.
    # - Create coverage for whole genome.
    # - Using fragment_size: None
    #     - https://htseq.readthedocs.io/en/master/tss.html#using-the-full-coverage
    # - (Optional) load corresponding `INPUT` BAM files and subtract INPUT from ChIP coverage
    # - Visualize coverage for specficic genomic regions: htc genes
    #     - `chr I   - 3,765,776 - 3,776,547`
    #     - `chr I  - 1,619,269 - 1,629,017`
    #     - `chr III - 1,093,002 - 1,105,987`
    # --------------------------
    
    # -------------------------------------------------------------
    # --  1. Compute Coverage: whole genome (might take a while) --
    # -------------------------------------------------------------
    
    if ( param_seq_type == 'ChIP' ):
        stranded = False
    else:
        stranded = True
    
    ## Compute coverage sample:
    int_cvg, fract_cvg, counts = coverage_vectors(param_in_bam, seq_type=param_seq_type, typecode='d')

    # However, a proper genome browser gives a better impression of the data.
    # The following commands write two `BedGraph` (Wiggle) files, one for the **plus** and one for the **minus strands**:
    if ( param_seq_type == 'ChIP' ) and ( isinstance(param_in_input, type(None)) ):

        int_cvg.write_bedgraph_file(os.path.join(out_dir, "int_" + sample_id + "_coverage.wig"), ".")
        fract_cvg.write_bedgraph_file(os.path.join(out_dir, "frac_" + sample_id + "_coverage.wig"), ".")

    elif ( param_seq_type == 'RNA' ):

        int_cvg.write_bedgraph_file(os.path.join(out_dir, "int_plus_" + sample_id + "_coverage.wig"), "+")
        fract_cvg.write_bedgraph_file(os.path.join(out_dir, "frac_plus_" + sample_id + "_coverage.wig"), "+")

        int_cvg.write_bedgraph_file(os.path.join(out_dir, "int_minus_" + sample_id + "_coverage.wig"), "-")
        fract_cvg.write_bedgraph_file(os.path.join(out_dir, "frac_minus_" + sample_id + "_coverage.wig"), "-")

    # These two files can then be viewed in a **genome browser** (e.g. [`IGB`](http://igb.bioviz.org/) or
    # [`IGV`](http://www.broadinstitute.org/igv/)), alongside the annotation from a GFF file.

    # int_cvg_file = os.path.join(out_dir, "int_cvg_" + sample_id + "_coverage.obj")
    # fract_cvg_file = os.path.join(out_dir, "frac_cvg_" + sample_id + "_coverage.obj")

    # # Store Objects:
    # pickle.dump(int_cvg, open(int_cvg_file, 'wb') )
    # pickle.dump(fract_cvg, open(fract_cvg_file, 'wb') )
    
    # Import Objects:
    # int_cvg = pickle.load(open(int_cvg_file, 'rb'))
    # fract_cvg = pickle.load(open(fract_cvg_file, 'rb'))
    
    ## ---------------------------------------------
    ## (Optional) Subtract INPUT from ChIP coverage
    ## ---------------------------------------------
    
    if not isinstance(param_in_input, type(None)):
        
        assert not isinstance(param_norm_factor, type(None))

        ## Compute coverage of corresponding INPUT ChIP sample
        input_int_cvg, input_fract_cvg, input_counts = coverage_vectors(param_in_input, seq_type=param_seq_type, typecode='d')
        #print(param_norm_factor)

        norm_factor = float(param_norm_factor)

        ## Integer: Subtract INPUT from ChIP coverage
        for _iv, _value in input_int_cvg.steps():
            int_cvg[ _iv ] += - norm_factor * _value

        ## Fractional: Subtract INPUT from ChIP coverage
        for _iv, _value in input_fract_cvg.steps():
            fract_cvg[ _iv ] += - norm_factor * _value

        ## Store as bedfiles
        int_cvg.write_bedgraph_file(os.path.join(out_dir, "int_subtracted_INPUT_" + sample_id + "_coverage.wig"), ".")
        fract_cvg.write_bedgraph_file(os.path.join(out_dir, "frac_subtracted_INPUT_" + sample_id + "_coverage.wig"), ".")

    # -------------------------------------------------------------
    # ---  2 . Visualize Coverage (at specific genomic regions) ---
    # -------------------------------------------------------------

    # - Define characteristics for `coverage interval`: 
    #     - Subtelomeric & telomeric genes tlh, SPAC212.10, etc...
    #     - Centromeric repeats dg, dh
    #     - Mating type region (MTR) MAT, other regions, etc...
    
    ## instead we can follow:
    # - https://htseq.readthedocs.io/en/master/tss.html#using-indexed-bam-files
    i_1 = {'region_id':'Subtelomeric_and_Telomeric_region', 'chrom':'I', 'start': 1, 'end': 55000}
    i_2 = {'region_id':'Centromeric_repeats', 'chrom':'I', 'start': 3748000, 'end': 3766000}
    #i_3 = {'region_id':'Mating_type_region', 'chrom':'II', 'start': 2113000, 'end': 2138000}
    i_3 = {'region_id':'Mating_type_region_chromosome', 'chrom':'mating_type_region', 'start': 1, 'end': 20128} ## whole chromosome
    #genomic_intervals = [i_1, i_2]
    genomic_intervals = [i_1, i_2, i_3]
    
    ## - **Fragment size** smoothes coverage profile:
    #fragment_size = 500
    #fragment_size = 200
    fragment_size = None
    
    ## - Get Data Frames **genomic features** in regions of interest:
    genomic_features_dfs = get_genomic_features_dfs(genomic_intervals, param_in_gtf, feature_types=['gene'])
                     
    ## ----------------------
    ## A. Regions of Interest: Coverage Vectors
    ## ----------------------
    
    fract_cvg = coverage_genomic_intervals(param_in_bam, genomic_intervals, count_type="frac", stranded=stranded, fragment_size=fragment_size)
    int_cvg = coverage_genomic_intervals(param_in_bam, genomic_intervals, count_type="int", stranded=stranded, fragment_size=fragment_size)

    ## (Optional) Subtract INPUT from ChIP coverage
    if not isinstance(param_in_input, type(None)):

        ## Compute coverage of corresponding INPUT ChIP sample
        input_fract_cvg = coverage_genomic_intervals(param_in_input, genomic_intervals, count_type="frac", stranded=stranded, fragment_size=fragment_size)
        input_int_cvg = coverage_genomic_intervals(param_in_input, genomic_intervals, count_type="int", stranded=stranded, fragment_size=fragment_size)

        norm_factor = float(param_norm_factor)

        ## Integer: Subtract INPUT from ChIP coverage
        for _iv, _value in input_int_cvg.steps():
            int_cvg[ _iv ] += - norm_factor * _value

        ## Fractional: Subtract INPUT from ChIP coverage
        for _iv, _value in input_fract_cvg.steps():
            fract_cvg[ _iv ] += - norm_factor * _value
        
        coverage_plot_dir = os.path.join(out_dir, "INPUT_subtracted_coverage_plots")
        if not os.path.exists(coverage_plot_dir):
            os.makedirs(coverage_plot_dir)

    ## -------------
    ## Visualization
    ## -------------

    #plot_coverage_genomic_regions(genomic_intervals, fract_cvg, genomic_features_dfs, stranded = stranded, annotate_features = visualize_gene_names, y_scale="log", prefix='frac_', out_dir=coverage_plot_dir)
    plot_coverage_genomic_regions(genomic_intervals, fract_cvg, genomic_features_dfs, stranded = stranded, annotate_features = visualize_gene_names, y_scale="linear", prefix='frac_', out_dir=coverage_plot_dir)
    
    #plot_coverage_genomic_regions(genomic_intervals, int_cvg, genomic_features_dfs, stranded = stranded, annotate_features = visualize_gene_names, y_scale="log", prefix='int_', out_dir=coverage_plot_dir)
    plot_coverage_genomic_regions(genomic_intervals, int_cvg, genomic_features_dfs, stranded = stranded, annotate_features = visualize_gene_names, y_scale="linear", prefix='int_', out_dir=coverage_plot_dir)

    # ## ------------------------
    # ## B. Heterochormatic Genes
    # ## ------------------------
    # 
    # htc_genes_dict = {}
    # 
    # # select gtf file: instantiate GFF_Reader Object
    # gtf = HTSeq.GFF_Reader(param_in_gtf)
    # 
    # ## loop over all features in gtf file - Select Heterochromatic Genes
    # for feature in gtf:
    # 
    #     ## Store all exons in our `GenomicArrayOfSets`
    #     #if feature.type == "exon":
    #     if feature.type == "gene":
    # 
    #         ## identify each `exon` feature by parent transcript/gene
    #         #gene_id = feature.attr["Parent"].split(':')[1][:-2]  ## get rid of '.1'
    #         gene_id = feature.attr["ID"] ## get rid of '.1'
    # 
    #         if gene_id in htc_genes:
    # 
    #             # Is this the first time we see this gene?
    #             if gene_id not in htc_genes_dict:
    #                 # If so, add to the 'genes_dict' an empty list
    #                 htc_genes_dict[gene_id] = list()
    # 
    #             # add the feature to the gene list
    #             htc_genes_dict[gene_id].append(feature)
    # 
    # # We can plot an excerpt of this with:
    # #feature_name = 'MAT3'
    # #plot_feature_coverage(feature_name, htc_genes_dict, out_dir=out_dir)
    # 
    # #import pdb; pdb.set_trace()
    # 
    # # - Loop over each feature **individually**
    # #for ff in htc_genes_dict.keys():
    # #    plot_feature_coverage(ff, htc_genes_dict, out_dir=out_dir)
    # 
    # #n_rows = 2
    # n_rows = 2
    # 
    # #n_cols = 3
    # n_cols = 4
    # 
    # ## Fractional Counts: (weighted by 1/NH)
    # subplots_feature_coverage(fract_cvg, htc_genes_dict, n_rows=n_rows, n_cols=n_cols, fig_size=(100, 50), prefix='fract_', out_dir=out_dir)
    # ## Integer Counts:
    # subplots_feature_coverage(int_cvg, htc_genes_dict, n_rows=n_rows, n_cols=n_cols, fig_size=(100, 50), prefix='int_', out_dir=out_dir)

