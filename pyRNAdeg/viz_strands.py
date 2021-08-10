#!/usr/bin/env python

import os
import logging

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from sklearn import linear_model

import pandas as pd
#import seaborn as sns

import pdb

logger = logging.getLogger(__name__)

# Define colors
myblue = (0.0 / 255, 153.0 / 255, 255.0 / 255)
myviol = (204.0 / 255, 0.0 / 255, 204.0 / 255)
myviol_dark = (124.0 / 255, 0.0 / 255, 124.0 / 255)
mysalmon = 'xkcd:salmon'

gene_id_col = 'gene_id'
#gene_id_col = 'gene-id'

#gene_name_col = 'gene_name'
#gene_name_col = 'gene-name'
gene_name_col = 'gene_id' ## for gff_v2 use 'gene_id' column instead


## ----------
## Parameters
## ----------

include_mat_locus_vis = False
#include_mat_locus_vis = True

########################################################
####   Divide Heterocrhomatic genes in 4 groups     ####
########################################################

# -------
# Group 1 - dg, dh: Centromeric Genes `deg1`
# -------

#deg1 = ['dh1', 'dg1'] ## unstranded
#deg1 = ['dh1+',  'dh1-', 'dg1+', 'dg1-']
deg1 = ['dg1a', 'dg1b', 'dh1a', 'dh1b'] ## gff_v2
deg1_names = {'dg1a': 'dg1+', 'dg1b': 'dg1-',
              'dh1a': 'dh1+', 'dh1b': 'dh1-'}

# --------
# Group 2 - tlh, SPAC212.10: Subtelomeric Genes `deg2`
# --------

#deg2 = ['tlh1', 'SPAC212.10'] ## unstranded
#deg2 = ['tlh1+', 'tlh1-', 'SPAC212.10']
deg2 = ['SPAC212.11', 'SPAC212.11b', 'SPAC212.10', 'SPAC212.10b'] ## gff_v2
deg2_names = {'SPAC212.11': 'tlh1-', 'SPAC212.11b': 'tlh1+',
              'SPAC212.10': 'SPAC212.10-', 'SPAC212.10b': 'SPAC212.10+'}

# --------
# Group 3: Mating type region (MTR) Genes `deg3`
# --------

# deg3 = ['MAT2', 'MAT3', 'MAT1']
# deg3_names = {'MAT2': 'MAT2',
#               'MAT3': 'MAT3',
#               'MAT1': 'MAT1'}
              
deg3 = ["SPMTR.01", "SPMTR.02",
        "FP565355_region_1..2120", "FP565355_region_3203..3259", "FP565355_region_3260..3394",
        "FP565355_region_3395..4497", "FP565355_region_4498..4556", "FP565355_region_9170..13408",
        "FP565355_region_15417..15473", "FP565355_region_15474..15608", "FP565355_region_15609..16735",
        "FP565355_region_16736..16794", "FP565355_region_18009..20128"]
deg3_names = {"SPMTR.01": "mat2-Pc",
              "SPMTR.02": "mat2-Pi",
              #"FP565355_region_1..2120": "region_1", "FP565355_region_3203..3259": "region_2", "FP565355_region_3260..3394": "region_3",
              "FP565355_region_1..2120": "1", "FP565355_region_3203..3259": "2", "FP565355_region_3260..3394": "3",
              #"FP565355_region_3395..4497": "region_4",  "FP565355_region_4498..4556": "region_5", "FP565355_region_9170..13408": "region_6", 
              "FP565355_region_3395..4497": "4",  "FP565355_region_4498..4556": "5", "FP565355_region_9170..13408": "6", 
              #"FP565355_region_15417..15473": "region_7", "FP565355_region_15474..15608": "region_8", "FP565355_region_15609..16735": "region_9", 
              "FP565355_region_15417..15473": "7", "FP565355_region_15474..15608": "8", "FP565355_region_15609..16735": "9", 
              #"FP565355_region_16736..16794": "region_10", "FP565355_region_18009..20128": "region_11"}
              "FP565355_region_16736..16794": "10", "FP565355_region_18009..20128": "11"}
              
# --------
# Group 4: rest Heterochromatic Genes (`non_degraded`)
# --------

non_degraded = ['SPAC212.09c', 'SPNCRNA.70', 'SPAC212.08c', 'SPAC212.07c', 'SPAC212.12', 'SPAC212.06c',
                'SPAC212.04c', 'SPAC212.03', 'SPAC212.02', 'SPAC212.01c', 'SPAC977.01', 'SPAC977.18',
                'SPAC977.02', 'SPAC977.03', 'SPAC977.04', 'SPAC212.05c',
                ]


def pos_maker(start, groups, samples, dist=5, pad=.9):

    positions = []
    n_samples = len(samples)

    # for each group - as many 'samples'
    for ii in range(groups):

        start_seq = start + (pad * ii)
        stop_seq = start_seq + (dist * n_samples)

        #pos_seq = np.arange(start_seq, stop_seq, step=dist)
        #print('start:', start_seq, 'stop:', stop_seq, 'step:', dist)

        pos_seq = np.linspace(start_seq, stop_seq, num=n_samples)
        #print('start:', start_seq, 'stop:', stop_seq, 'num:', 6)
        #print(pos_seq, '\n')

        assert len(pos_seq) == n_samples
        positions.append(pos_seq)

    return positions


def prepare_4cat_data(df, samples):
    
    # -----------
    # Gene counts: Protein Coding genes
    # -----------

    data1 = [] # 0
    for sample in samples:
        ## Create list with all expressions for Protein Coding genes
        #data1.append(list(df[df['category'] == 'gene'][sample].dropna()))
        data1.append(list(df[(df['bio_type'] == 'mRNA') & (df['category'] == 'gene')][sample].dropna()))

    # -------------
    # dg, dh counts: Centromeric Genes `deg1`
    # -------------

    data = [] # 1
    names = []
    for sample in samples:
        ## Create list with all expressions for Centromeric Genes
        data.append(list(df[df[gene_name_col].isin(deg1)][sample]))
        ## Add names for annotation
        gene_names = df[df[gene_name_col].isin(deg1)][gene_name_col].tolist()
        #names.append(gene_names)
        names.append([deg1_names[nn] for nn in gene_names])
    data2 = (data, names)

    # ----------------------
    # tlh, SPAC212.10 counts: Subtelomeric Genes `deg2`
    # ----------------------

    data = [] # 2
    names = []
    for sample in samples:
        ## Create list with all expressions for Centromeric Genes
        data.append(list(df[df[gene_name_col].isin(deg2)][sample]))
        ## Add names for annotation
        gene_names = df[df[gene_name_col].isin(deg2)][gene_name_col].tolist()
        #names.append(gene_names)
        names.append([deg2_names[nn] for nn in gene_names])
    data3 = (data, names)

    # --------------------------------------
    # rest Heterochromatic + mat gene counts: `non_degraded`
    # --------------------------------------

    data4 = [] # 3
    for sample in samples:
        ## Create list with all expressions for rest Heterochromatic + mat gene counts: `non_degraded`
        data4.append(list(df[df[gene_name_col].isin(non_degraded)][sample].dropna()))

    return data1, data2, data3, data4


def prepare_5cat_data(df, samples):
    
    # -----------
    # Gene counts: Protein Coding genes
    # -----------

    data1 = [] # 0
    for sample in samples:
        ## Create list with all expressions for: 'Protein Coding genes'
        #data1.append(list(df[df['category'] == 'gene'][sample].dropna()))
        data1.append(list(df[(df['bio_type'] == 'mRNA') & (df['category'] == 'gene')][sample].dropna()))

    # -------------
    # dg, dh counts: Centromeric Genes `deg1`
    # -------------

    data = []
    names = [] # 1
    for sample in samples:
        ## Create list with all expressions for: 'Centromeric genes'
        data.append(list(df[df[gene_name_col].isin(deg1)][sample]))
        ## Add names for annotation
        gene_names = df[df[gene_name_col].isin(deg1)][gene_name_col].tolist()
        #names.append(gene_names)
        names.append([deg1_names[nn] for nn in gene_names])
    data2 = (data, names)

    # ----------------------
    # tlh, SPAC212.10 counts: Subtelomeric Genes `deg2`
    # ----------------------

    data = []
    names = [] # 2
    for sample in samples:
        ## Create list with all expressions for: 'Sub-telomeric genes'
        data.append(list(df[df[gene_name_col].isin(deg2)][sample]))
        ## Add names for annotation
        gene_names = df[df[gene_name_col].isin(deg2)][gene_name_col].tolist()
        #names.append(gene_names)
        names.append([deg2_names[nn] for nn in gene_names])
    data3 = (data, names)

    # ------------------------
    # MAT1, MAT2, MAT3 counts: Mating type region (MTR) Genes `deg3`
    # ------------------------

    # data = []
    # names = []
    # for sample in samples:
    #     ## Create list with all expressions for: 'Mating Type Region genes'
    #     data.append(list(df[df[gene_name_col].isin(deg3)][sample]))
    #     ## Add names for annotation
    #     gene_names = df[df[gene_name_col].isin(deg3)][gene_name_col].tolist()
    #     #names.append(gene_names)
    #     names.append([deg3_names[nn] for nn in gene_names])
    # data4 = (data, names)

    data4 = [] # 3
    for sample in samples:
        ## Create list with all expressions for: 'rest Heterochromatic genes' (`non_degraded`)
        data4.append(list(df[df[gene_name_col].isin(deg3)][sample].dropna()))

    # ---------------------
    # rest Heterochromatic: `non_degraded`
    # ---------------------

    data5 = [] # 4
    for sample in samples:
        ## Create list with all expressions for: 'rest Heterochromatic genes' (`non_degraded`)
        data5.append(list(df[df[gene_name_col].isin(non_degraded)][sample].dropna()))

    #import pdb; pdb.set_trace()

    return data1, data2, data3, data4, data5


def single_box_plot(data, position, color, widths=None):

    #pdb.set_trace()
    plt.boxplot(
        data,
        positions=position,
        patch_artist=True, 
        widths=widths,
        #meanline=True,
        boxprops=dict(facecolor='white', color=color, linewidth=7),
        flierprops=dict(color=color, markeredgecolor=color, markersize=20),
        capprops=dict(color=color, linewidth=9),
        whiskerprops=dict(color=color, linewidth=5),
        medianprops=dict(color=color, linewidth=5),
        #medianprops=dict(color='red', linewidth=10),
        )

def single_violin_plot(data, position, color, widths=None):

    sns.violinplot(data,
                   positions=position,
                   patch_artist=True,
                   widths=widths,
                   boxprops=dict(facecolor='white', color=color, linewidth=7),
                   flierprops=dict(color=color, markeredgecolor=color, markersize=20),
                   capprops=dict(color=color, linewidth=9),
                   whiskerprops=dict(color=color, linewidth=5),
                   medianprops=dict(color=color, linewidth=5),
                   )


def annotated_scatter_plot(axes, data, names, pos, samples, color, annotate=True):
    """
    used to plot individual points together with box plots.

    :param axes:
    :param data:
    :param names:
    :param pos:
    :param samples:
    :param color:
    :param annotate:
    :return:
    """

    ## iterate over samples
    for i in range(len(samples)):

        idx = i - 1
        #idx = i

        y = data[idx]
        ## x-position fixed, repeat value for each point y
        x = [pos[idx]] * len(y)

        name = names[idx]

        ## iterate over points (genes) in scatter plot
        for xx, yy, nn in zip(x, y, name):

            if annotate:
                #axes.annotate(nn, (xx, yy), fontsize=30)
                axes.annotate(nn, (xx, yy), fontsize=25, alpha=0.8)

            if ('dg' in nn) or ('tlh' in nn):
                marker_style = "o"
            else:
                marker_style = "D"

            ## change for if n.str.contains('a') or 'b' / 'minus'or 'plus'
            if '-' in nn:
                plt.scatter(xx, yy, color=color, marker=marker_style, s=500, alpha=0.3, linewidth=3, linestyle='dashed', edgecolors='k')
                #axes.scatter(aa, bb, color=color, marker="x", s=800)

            elif '+' in nn:
                plt.scatter(xx, yy, color=color, marker=marker_style, s=500, alpha=0.8, linewidth=2, edgecolors='k')
                #axes.scatter(x, y, color=color, marker='+', s=800)

            else:
                plt.scatter(xx, yy, color=color, marker=marker_style, s=500, alpha=0.75, edgecolors='k')
                #axes.scatter(x, y, color=color, marker='o', s=800)


def multi_4cat_box_plot(data, samples, labels, outpath, figsize=(45, 25), dist=6, title=None, title_font_size=90,
                        y_label=None, y_lim=(-6, 6), format='png', annotate=False, include_mat_locus_vis=False,
                        widths=None, hlines=None, xlable_size=80, ylable_size=70):
    """
    4 categories

    :param data:
    :param samples:
    :param labels:
    :param outpath:
    :param figsize:
    :param dist:
    :param title:
    :param title_font_size:
    :param y_label:
    :param y_lim:
    :param format:
    :param annotate:
    :param widths:
    :param hlines:
    :param xlable_size:
    :param ylable_size:
    :return:
    """

    # Assign box positions:
    pos1, pos2, pos3, pos4 = pos_maker(1, 4, samples, dist=dist)

    # Assign colors
    c1 = 'dimgray' ## protein coding genes
    c2 = myviol ## centromeric genes
    c3 = myblue ## rest of heterochromatic genes
    c4 = 'g' ## subtelomeric genes

    ## ---------
    ## Box Plots
    ## ---------

    fig, ax = plt.subplots(figsize=figsize)

    ax.spines['top'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)

    ax.set_ylabel(y_label, fontsize=ylable_size, labelpad=20)
    ax.set_xticklabels(labels, fontstyle='italic', fontname='sans-serif')
    ax.tick_params(axis='x', which='major', pad=30)
    plt.xticks(fontsize=xlable_size)
    plt.yticks(fontsize=ylable_size, fontname='sans-serif')

    if title:
        plt.suptitle(title, fontsize=title_font_size, y=.95)

    # Plot gene box groups: Protein Coding genes
    single_box_plot(data[0], pos1, c1, widths=widths)

    # Plot the actual points for bp2 samples: Centromeric genes
    annotated_scatter_plot(ax, data[1][0], data[1][1], pos2, samples, c2, annotate=annotate)

    # Plot heterochromatine gene box groups: rest of Heterochromatic genes
    single_box_plot(data[3], pos3, c3, widths=widths)

    # Plot the actual points for bp4 samples: Subtelomeric genes
    annotated_scatter_plot(ax, data[2][0], data[2][1], pos4, samples, c4, annotate=annotate)

    # Plot horizontal guide-lines
    if hlines:
        start = min(pos1) - 2
        end = max(pos4) + 2
        for x, color in hlines:
            plt.hlines(x, start, end, colors=color, linestyles='-')

    # Ticks and limits
    plt.ylim(y_lim)
    plt.xlim(min(pos1) - 2, max(pos4) + 2)
    ax.set_xticks(pos2)
    ax.set_xticklabels(labels)

    # Grid
    plt.rc('grid', linestyle="--", color='silver', linewidth=3)
    ax.yaxis.grid(True)

    pdf_path = outpath.split('.')[0] + '.pdf'
    plt.savefig(outpath, format=format, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.show()


def my_multi_4cat_box_plot(df, samples, labels, outpath, figsize=(45, 25), dist=6, title=None, title_font_size=90,
                        y_label=None, y_lim=(-6, 6), format='png', annotate=False, include_mat_locus_vis=False,
                        widths=None, hlines=None, xlable_size=80, ylable_size=70):

    if include_mat_locus_vis:
        data = prepare_5cat_data(df, samples)

        # Assign box positions:
        pos1, pos2, pos3, pos4, pos5 = pos_maker(1, 5, samples, dist=dist)

    else:
        data = prepare_4cat_data(df, samples)

        # Assign box positions:
        pos1, pos3, pos4, pos5 = pos_maker(1, 4, samples, dist=dist)

    # Assign colors
    c1 = 'dimgray' ## protein coding genes
    c2 = myviol ## centromeric genes
    c3 = myblue ## rest of heterochromatic genes
    c4 = 'g' ## subtelomeric genes
    c5 = mysalmon ## mating type region

    ## ---------
    ## Box Plots
    ## ---------

    fig, ax = plt.subplots(figsize=figsize)

    ax.spines['top'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)

    ax.set_ylabel(y_label, fontsize=ylable_size, labelpad=20)
    ax.set_xticklabels(labels, fontstyle='italic', fontname='sans-serif')
    ax.tick_params(axis='x', which='major', pad=30)
    plt.xticks(fontsize=xlable_size)
    plt.yticks(fontsize=ylable_size, fontname='sans-serif')

    if title:
        plt.suptitle(title, fontsize=title_font_size, y=.95)

    # ------------
    # 1st Box Plot: Protein Coding genes
    # ------------

    single_box_plot(data[0], pos1, c1, widths=widths)

    # ------------------------
    # (Optional) 2nd Box Plot: Mating Type regions genes
    # ------------------------

    if include_mat_locus_vis:
        single_box_plot(data[3], pos2, c5, widths=widths)

    # Plot the actual points for bp2 samples: Centromeric genes
    annotated_scatter_plot(ax, data[1][0], data[1][1], pos3, samples, c2, annotate=annotate)

    # -------------
    # 3rd Box Plot: rest of Heterochromatic genes
    # -------------

    if include_mat_locus_vis:
        single_box_plot(data[4], pos4, c3, widths=widths)

    else:
        single_box_plot(data[3], pos4, c3, widths=widths)

    # Plot the actual points for bp4 samples: Sub-telomeric genes
    annotated_scatter_plot(ax, data[2][0], data[2][1], pos5, samples, c4, annotate=annotate)

    # Plot horizontal guide-lines
    if hlines:
        start = min(pos1) - 2
        end = max(pos5) + 2
        for x, color in hlines:
            plt.hlines(x, start, end, colors=color, linestyles='-')

    # Ticks and limits
    plt.ylim(y_lim)
    plt.xlim(min(pos1) - 2, max(pos5) + 2)
    ax.set_xticks(pos3)
    ax.set_xticklabels(labels)

    # Grid
    plt.rc('grid', linestyle="--", color='silver', linewidth=3)
    ax.yaxis.grid(True)

    pdf_path = outpath.split('.')[0] + '.pdf'
    plt.savefig(outpath, format=format, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.show()


def regress_pair(df, sample1, sample2, xlim=12):

    X = np.array(list(df[sample1])).reshape(-1, 1)
    Y = np.array(list(df[sample2])).reshape(-1, 1)
    assert len(X) == len(Y)

    regr = linear_model.LinearRegression()
    regr.fit(X, Y)

    coef = regr.coef_[0][0]
    intercept = regr.intercept_[0]
    reg_data = [[float(i) for i in X], [float(i * coef + intercept) for i in X]]
    reg_data[0].append(xlim)
    reg_data[1].append(float(xlim * coef + intercept))

    return reg_data, coef, intercept

def linear_regress(df, sample1, sample2, fit_intercept=True, xlim=12):
    
    # one shape dimension can be -1, the value is inferred from the length of the array
    X = np.array(df[sample1].to_list()).reshape(-1, 1)
    y = np.array(df[sample2].to_list()).reshape(-1, 1)
    assert len(X) == len(y)

    reg = linear_model.LinearRegression(fit_intercept=fit_intercept)
    reg.fit(X, y)

    lin_model = {}
    # estimated coefficients for the linear regression problem
    lin_model['a'] = reg.coef_[0][0]
    # independent term in the linear model. Set to 0.0 if fit_intercept = False.
    lin_model['b'] = reg.intercept_[0] if fit_intercept else 0
    # the coefficient of determination R^2 of the prediction.
    lin_model['r_score'] = reg.score(X, y)
    
    if isinstance(xlim, type(None)):
        xlim = max(X)[0]
        
    #lin_model['x'] = [0, xlim]
    lin_model['x'] = [2, 12]
    y = reg.predict(np.array(lin_model['x']).reshape(-1, 1)).flatten()
    lin_model['y'] = list(y)
    
    return lin_model
    
def average_slope(df, sample1, sample2, fit_intercept=True, xlim=12):
    
    # remove zeros - before dividing
    df = df[~(df[[sample1, sample2]] == 0).any(axis=1)]
    # compute slope for each point
    df['a'] = df[sample2] / df[sample1]
    
    lin_model = {}
    # estimated coefficients for the linear regression problem
    lin_model['a'] = df['a'].mean()
    # independent term in the linear model.
    lin_model['b'] = 0
    # the coefficient of determination R^2 of the prediction.
    lin_model['r_score'] = -np.inf
    
    if isinstance(xlim, type(None)):
        xlim = df[sample1].max()
    
    lin_model['x'] = [0, xlim]
    y = lin_model['a'] * np.array(lin_model['x'])
    lin_model['y'] = list(y)
    
    return lin_model
    
def scatter_plot_gene_params(nn):

    param_dict = {}

    if ('dg' in nn) or any(x in nn for x in ['tlh', 'SPAC212.11']):
        param_dict['marker_style'] = "o"
    else:
        param_dict['marker_style'] = "D"

    # default params
    param_dict['linewidth'] = 3
    param_dict['edgecolors'] = 'k'
    param_dict['linestyle'] = 'solid'

    #import pdb; pdb.set_trace()
    # check strandness
    if ('+' in nn):
        param_dict['alpha'] = 0.8
        param_dict['linewidth'] = 2

    elif ('-' in nn):
        param_dict['alpha'] = 0.3
        param_dict['linestyle'] = 'dashed'

    else:
        raise Exception('Unknown strand!')

    return param_dict

def scatter_plot(df, sample1, sample2, out_dir, fname, format=None, regressor=False, fit_intercept=True,
                 annotate=True, include_mat_locus_vis=False, include_rest_heterochromatic_repeats=False,
                 include_other_genes=None,
                 xlabel=None, ylabel=None, ribo_color='gray', xlim=12, ylim=12):

    # Plots genes + ribosomal genes + centromeric, and subtelomeric genes
    # Plot genes:
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.spines['top'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    
    # Plot regression line
    if regressor:
        #rgr_data, _, _ = regress_pair(df, sample1, sample2, xlim=xlim)
        #ax.plot(rgr_data[0], rgr_data[1], color=myviol_dark)  # , label='Average degradation line')
        lin_model = linear_regress(df, sample1, sample2, fit_intercept=fit_intercept, xlim=xlim)
        #lin_model = average_slope(df, sample1, sample2, fit_intercept=fit_intercept, xlim=xlim)
        #ax.plot(lin_model['x'], lin_model['y'], color=myviol_dark)
        
        if (lin_model['b'] > 0):
            label_leg = r'y={0:.1f}x+{1:.1f} ($R^2$={2:.2f})'.format(lin_model['a'], lin_model['b'], lin_model['r_score'])
        else:
            label_leg = r'y={0:.1f}x-{1:.1f} ($R^2$={2:.2f})'.format(lin_model['a'], -lin_model['b'], lin_model['r_score'])
        ax.plot(lin_model['x'], lin_model['y'], color=myviol_dark, label=label_leg)
        plt.legend(loc='upper left', fontsize=35)
        #import pdb; pdb.set_trace()
    
    # -------------------------
    # Plot Protein Coding Genes
    # -------------------------

    #gene_df = df[(df['type'] == 'gene') & (df['category'] != 'ribosomal')]
    gene_df = df[(df['bio_type'] == 'mRNA') & (df['category'] == 'gene')]
    ax.plot(gene_df[sample1], gene_df[sample2],
            marker='o', markerfacecolor='white', markeredgecolor='gray',
            linestyle='', ms=20, label='genes',
            zorder=-1)

    # -------------
    # Plot RP Genes
    # -------------

    #gene_df = df[(df['type'] == 'gene') & (df['category'] == 'ribosomal')]
    gene_df = df[df['category'] == 'ribosomal']
    ax.plot(gene_df[sample1], gene_df[sample2],
            #marker='o', markerfacecolor=ribo_color, markeredgecolor='black',
            marker='o', markerfacecolor='white', markeredgecolor='gray',
            linestyle='', ms=20, label='ribosomal proteins',
            zorder=-1)

    # ----------------------
    # Plot centromeric genes: `deg1`
    # ----------------------

    c_df = df[df[gene_name_col].isin(deg1)]
    if annotate:
        for index, row in c_df.iterrows():
            #ax.annotate(row[gene_name_col], (row[sample1], row[sample2]),
            ax.annotate(deg1_names[row[gene_name_col]], (row[sample1], row[sample2]),
                        size=25, alpha=0.75, rotation=45,)
                        #arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))  # rotation='vertical'

    #for xx, yy, nn in zip(c_df[sample1], c_df[sample2], c_df['gene_name']):
    for index, row in c_df.iterrows():
        param_dict = scatter_plot_gene_params(deg1_names[row[gene_name_col]])
        ax.scatter(
            #xx, yy,
            row[sample1], row[sample2],
            color=myviol,
            marker=param_dict['marker_style'],
            alpha=param_dict['alpha'],
            linewidth=param_dict['linewidth'], linestyle=param_dict['linestyle'],
            edgecolors=param_dict['edgecolors'],
            zorder=1,
            s=500,
            label='${dg, dh}$'
        )

    # ------------------------
    # Plot subtelomeric genes: `deg2`
    # ------------------------

    s_df = df[df[gene_name_col].isin(deg2)]
    if annotate:
        for index, row in s_df.iterrows():
            #ax.annotate(row[gene_name_col], (row[sample1], row[sample2]),
            ax.annotate(deg2_names[row[gene_name_col]], (row[sample1], row[sample2]),
                        size=25, alpha=0.75, rotation=45,)
                        # arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))  # rotation='vertical'
    #import pdb; pdb.set_trace()

    #for xx, yy, nn in zip(s_df[sample1], s_df[sample2], s_df['gene_name']):
    for index, row in s_df.iterrows():
        param_dict = scatter_plot_gene_params(deg2_names[row[gene_name_col]])
        ax.scatter(
            #xx, yy,
            row[sample1], row[sample2],
            color='g',
            marker=param_dict['marker_style'],
            alpha=param_dict['alpha'],
            linewidth=param_dict['linewidth'], linestyle=param_dict['linestyle'],
            edgecolors=param_dict['edgecolors'],
            zorder=1,
            s=500,
            label='${tlh, SPAC212.10}$'
        )


    # ------------------------------------
    # Plot Mating type region (MTR) genes: `deg3`
    # ------------------------------------

    if include_mat_locus_vis:
        m_df = df[df[gene_name_col].isin(deg3)]
        if annotate:
            for index, row in m_df.iterrows():
                #ax.annotate(row[gene_name_col], (row[sample1], row[sample2]),
                ax.annotate(deg3_names[row[gene_name_col]], (row[sample1], row[sample2]),
                            size=25, alpha=0.75, rotation=45,)
                            # arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))  # rotation='vertical'

        ax.plot(m_df[sample1], m_df[sample2],
                marker='o', markerfacecolor=mysalmon, markeredgecolor=mysalmon,
                linestyle='', ms=25, label='${MAT1, MAT2, MAT3}$')

    # ---------------------
    # rest Heterochromatic: `non_degraded`
    # ---------------------
    
    if include_rest_heterochromatic_repeats:
        rest_df = df[df[gene_name_col].isin(non_degraded)]
        ax.plot(rest_df[sample1], rest_df[sample2],
                marker='o', markerfacecolor=myblue, markeredgecolor=myblue,
                linestyle='', ms=25, label='${heterochromatic genes}$')
    
    if not isinstance(include_other_genes, type(None)):
        assert isinstance(include_other_genes, type(pd.DataFrame()))
        other_df = include_other_genes
        
        if annotate:
            for index, row in other_df.iterrows():
                import pdb; pdb.set_trace()
                ax.annotate(
                    row['gene_name'] if not row['gene_name'] else row[gene_name_col],
                    (row[sample1], row[sample2]),
                    size=25, alpha=0.75, rotation=45,)
                    # arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))  # rotation='vertical'

        ax.plot(other_df[sample1], other_df[sample2],
                marker='o', markerfacecolor='r', markeredgecolor='k',
                linestyle='', ms=25, label='${Regions of interest}$')

    plt.xticks(fontsize=70)
    plt.yticks(fontsize=70)

    plt.xlim(0, xlim)
    plt.ylim(0, ylim)

    plt_fname = fname
    if xlabel:
        plt.xlabel(xlabel, fontsize=80)
    if ylabel:
        plt.ylabel(ylabel, fontsize=80)

    '''
    if format:
        plt.savefig(os.path.join(out_dir, plt_fname), format=format, bbox_inches='tight')
    else:
        plt.savefig(os.path.join(out_dir, plt_fname), bbox_inches='tight')
    '''

    outpath= os.path.join(out_dir, plt_fname)
    pdf_path = outpath.split('.')[0] + '.pdf'
    plt.savefig(outpath, format=format, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.show()



if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('RNAdeg visualization module.')
    logger.info('-' * 40)





