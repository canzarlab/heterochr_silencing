#!/usr/bin/env python

import os
import logging

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from sklearn import linear_model

logger = logging.getLogger(__name__)

# Define colors
myblue = (0.0 / 255, 153.0 / 255, 255.0 / 255)
myviol = (204.0 / 255, 0.0 / 255, 204.0 / 255)
myviol_dark = (124.0 / 255, 0.0 / 255, 124.0 / 255)

# Define gene groups
deg1 = ['dh1', 'dg1']
deg2 = ['tlh1', 'SPAC212.10']
non_degraded = ['SPAC212.09c', 'SPNCRNA.70', 'SPAC212.08c', 'SPAC212.07c', 'SPAC212.12', 'SPAC212.06c',
                'SPAC212.04c', 'SPAC212.03', 'SPAC212.02', 'SPAC212.01c', 'SPAC977.01', 'SPAC977.18',
                'SPAC977.02', 'SPAC977.03', 'SPAC977.04', 'SPAC212.05c',
                'MAT2', 'MAT3', 'MAT1']


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


def pos_maker(start, groups, samples, dist=5, pad=.9):

    positions = []

    for i in range(groups):
        positions.append(np.arange(start + pad * i, len(samples) * dist + (start + pad * i), dist))

    return positions


def prepare_4cat_data(df, samples):

    # Gene counts
    data1 = []
    for sample in samples:

        data1.append(list(df[df['category'] == 'gene'][sample].dropna()))

    # dg, dh counts
    data = []
    names = []
    for sample in samples:
        data.append(list(df[df['gene-name'].isin(deg1)][sample]))
        names.append(list(df[df['gene-name'].isin(deg1)]['gene-id']))
    data2 = (data, names)

    # tlh, SPAC212.10 counts
    data = []
    names = []
    for sample in samples:
        data.append(list(df[df['gene-name'].isin(deg2)][sample]))
        names.append(list(df[df['gene-name'].isin(deg2)]['gene-id']))
    data3 = (data, names)

    # Heterochromatine + mat gene counts
    data4 = []
    for sample in samples:
        data4.append(list(df[df['gene-name'].isin(non_degraded)][sample].dropna()))

    return data1, data2, data3, data4


def single_box_plot(data, position, color, widths=None):

    plt.boxplot(data, positions=position, patch_artist=True, widths=widths,
                boxprops=dict(facecolor='white', color=color, linewidth=7),
                flierprops=dict(color=color, markeredgecolor=color, markersize=20),
                capprops=dict(color=color, linewidth=9),
                whiskerprops=dict(color=color, linewidth=5),
                medianprops=dict(color=color, linewidth=5),
                )


def annotated_scatter_plot(axes, data, names, pos, samples, color, annotate=True):

    for i in range(len(samples)):
        idx = i - 1
        y = data[idx]
        x = [pos[idx]] * len(y)
        name = names[idx]
        for a, b, n in zip(x, y, name):
            if annotate:
                axes.annotate(n, (a, b))

            plt.scatter(x, y, color=color, marker='o', s=800)


def multi_4cat_box_plot(data, samples, labels, outpath, figsize=(45, 25), dist=6, title=None, title_font_size=90,
                        y_label=None, y_lim=(-6, 6), format='png', annotate=False, widths=None,
                        hlines=None, xlable_size=80, ylable_size=70):

    # Assign box positions
    pos1, pos2, pos3, pos4 = pos_maker(1, 4, samples, dist=dist)
    # Assign colors
    c1 = 'dimgray'; c2 = myviol; c3 = myblue; c4 = 'g'

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

    # Plot gene box groups.
    single_box_plot(data[0], pos1, c1, widths=widths)
    # Plot the actual points for bp2 samples.
    annotated_scatter_plot(ax, data[1][0], data[1][1], pos2, samples, c2, annotate=annotate)
    # Plot heterochromatine gene box groups.
    single_box_plot(data[3], pos3, c3, widths=widths)
    # Plot the actual points for bp4 samples.
    annotated_scatter_plot(ax, data[2][0], data[2][1], pos4, samples, c4, annotate=annotate)

    # Plot horizontal guidlines
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


def scatter_plot(df, sample1, sample2, out_dir, fname, format=None, regressor=False, annotate=True,
                 xlabel=None, ylabel=None, ribo_color='gray', xlim=12, ylim=12):

    # Plots genes + ribosomal genes + centromeric, and subtelomeric genes
    # Plot genes:
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.spines['top'].set_linewidth(5)
    ax.spines['right'].set_linewidth(5)
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)

    # Plot genes
    gene_df = df[(df['type'] == 'gene') & (df['category']!='ribosomal')]
    ax.plot(gene_df[sample1], gene_df[sample2], marker='o', markerfacecolor='white', markeredgecolor='gray',
            linestyle='', ms=20,
            label='genes')

    # Plot RP genes
    gene_df = df[(df['type'] == 'gene') & (df['category']=='ribosomal')]
    # ax.plot(gene_df[sample1], gene_df[sample2], marker='o', markerfacecolor=ribo_color, markeredgecolor='black',
    ax.plot(gene_df[sample1], gene_df[sample2], marker='o', markerfacecolor='white', markeredgecolor='gray',
            linestyle='', ms=20, label='ribosomal proteins')

    # Plot regression line
    if regressor:
        rgr_data, _, _ = regress_pair(df, sample1, sample2, xlim=xlim)
        ax.plot(rgr_data[0], rgr_data[1], color=myviol_dark)  # , label='Average degradation line')

    # Plot centromeric genes:
    c_df = df[df['gene-name'].isin(deg1)]
    for index, row in c_df.iterrows():
        if annotate:
            ax.annotate(row['gene-name'], (row[sample1], row[sample2]), size=40)  # rotation='vertical')
    ax.plot(c_df[sample1], c_df[sample2], marker='o', markerfacecolor=myviol,
            markeredgecolor=myviol, linestyle='', ms=25, label='${dg, dh}$')

    # Plot subtelomeric genes:
    s_df = df[df['gene-name'].isin(deg2)]
    for index, row in s_df.iterrows():
        if annotate:
            ax.annotate(row['gene-name'], (row[sample1], row[sample2]), size=40)
    ax.plot(s_df[sample1], s_df[sample2], marker='o', markerfacecolor='g', markeredgecolor='g', linestyle='',
            ms=25, label='${tlh, SPAC212.10}$')

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





