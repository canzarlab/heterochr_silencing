#!/usr/bin/env python

import os
import logging

import pandas as pd
import pysam

## what is this?
from expression_data import gene_table_from_annot

pd.options.mode.chained_assignment = None

logger = logging.getLogger(__name__)


def coverage_dict_from_bed(in_bed):

    cov_dict = {}

    with open(in_bed) as itc:
        for line in itc:
            ls = line.split('\t')
            ch = ls[0].strip(' ')
            pos = ls[1].strip(' ')
            key = ch + ':' + pos
            cov_dict[key] = 0

    return cov_dict


def update_coverage_dict(coverage_dict, ch, start, end, inc):

    for i in range(start, end):
        key = ch + ':' + str(i)
        coverage_dict[key] += inc


def write_coverage_file(coverage_dict, out_file):

    with open(out_file, 'w+') as out:
        for coord, item in coverage_dict.items():
            ch = coord.split(':')[0]
            pos = coord.split(':')[1]
            line = ch + '\t' + str(pos) + '\t' + str(round(item, 2)) + '\n'
            out.write(line)


def bed_list_from_coverage_dict(cov_dict):

    bed_list = []

    for key in cov_dict.keys():
        x, y = key.split(':')
        bed_list.append((x, y, cov_dict[key]))

    return bed_list


def rrna_iter_from_gene_table(in_gdf):

    rrna_iter = in_gdf[in_gdf['bio_type'] == 'rRNA'][['chr', 'start', 'end']].values

    return rrna_iter


def rrna_count(in_file, iterator, reverse=True):

    rrna_counter = 0
    alignment_file = pysam.AlignmentFile(open(in_file, 'rb'))

    for ch, start, end in iterator:

        for read in alignment_file.fetch(region=ch, start=start, end=end):

            if reverse == read.is_reverse:
                nh = read.get_tag('NH')
                inc = 1 / nh
                rrna_counter += inc

    return round(rrna_counter, 2)


def normalize_coverage(cov_dict, rrna_counts=0):

    read_len = 50

    bed_list = bed_list_from_coverage_dict(cov_dict)
    df = pd.DataFrame(bed_list, columns=range(3))

    total = sum(df[2]) / read_len
    total -= rrna_counts
    scaling_factor = (total / 1000000)
    df[2] = round(df[2] / scaling_factor, 2)

    return df


def bam_coverage_dict(in_file, in_bed, strand=None):

    forward = True if not strand else  ('forward' == strand)
    reverse = True if not strand else ('reverse' == strand)

    coverage_dict = coverage_dict_from_bed(in_bed)
    st = pysam.AlignmentFile(open(in_file, 'rb'))

    for read in st.fetch(until_eof=True):

        r = read.is_reverse
        f = not r

        if (r and reverse) or (f and forward):

            nh = read.get_tag('NH')
            inc = 1.0 / nh
            ch = read.reference_name
            blocks = read.get_blocks()

            for start, end in blocks:
                update_coverage_dict(coverage_dict, ch, start + 1, end + 1, inc)

    return coverage_dict


def bam_file_coverage(sample, seq_type, in_file, in_gdf, in_bed, out_dir):

    if 'ChIP' in seq_type:

        strand = None
        cov_dict = bam_coverage_dict(in_file, in_bed, strand=strand)
        norm_df = normalize_coverage(cov_dict)
        write_coverage_file(cov_dict, os.path.join(out_dir, sample + '.coverage.txt'))
        norm_df.to_csv(os.path.join(out_dir, sample + '.norm.coverage.txt'), sep='\t', index=None, header=None)

    else:

        rrna_iter = rrna_iter_from_gene_table(in_gdf)

        strand = 'forward'
        cov_dict = bam_coverage_dict(in_file, in_bed, strand=strand)
        rrna_counts = rrna_count(in_file, rrna_iter, reverse=False)
        norm_df = normalize_coverage(cov_dict, rrna_counts=rrna_counts)

        write_coverage_file(cov_dict, os.path.join(out_dir, sample + '.forward.coverage.txt'))
        norm_df.to_csv(os.path.join(out_dir, sample + '.forward.norm.coverage.txt'),
                       sep='\t', index=None, header=None)

        strand = 'reverse'
        cov_dict = bam_coverage_dict(in_file, in_bed, strand=strand)
        rrna_counts = rrna_count(in_file, rrna_iter, reverse=True)
        norm_df = normalize_coverage(cov_dict, rrna_counts=rrna_counts)
        write_coverage_file(cov_dict, os.path.join(out_dir, sample + '.reverse.coverage.txt'))
        norm_df.to_csv(os.path.join(out_dir, sample + '.reverse.norm.coverage.txt'),
                       sep='\t', index=None, header=None)

    return


def compute_coverage(out_base, info_file, samples, in_gtf, in_bed):

    gdf = gene_table_from_annot(in_gtf)

    samples_df = pd.read_csv(info_file, sep='\t', names=['sample', 'seq', 'mut'])

    bam_dir = os.path.join(out_base, 'bams/')
    cov_dir = os.path.join(out_base, 'coverage/')

    for sample in samples:

        seq_type = samples_df[samples_df['sample'] == sample]['seq'].values[0]
        in_file = os.path.join(bam_dir, sample + '.bam')

        bam_file_coverage(sample, seq_type, in_file, gdf, in_bed, cov_dir)

    return


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to coverage module')

