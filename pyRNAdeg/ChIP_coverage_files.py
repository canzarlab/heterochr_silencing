#!/usr/bin/env python

import os
import logging
import argparse

import pandas as pd
import pysam


logger = logging.getLogger(__name__)


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--in_dir', '-d', type=str, help='Alignments folder (.bam)', required=True)
    #parser.add_argument(
    #    '--in_gdf', '-g', type=str, help='Gene annotation table (tab-delimited file).', required=True)
    parser.add_argument(
        '--in_bed', '-b', type=str, help='Bed scaffold file.', required=True)
    parser.add_argument(
        '--out_dir', '-o', type=str, default='.', help='Outp path (default: .)')
    parser.add_argument(
        '-s', action='store_true', help='Save report to file.')

    args = parser.parse_args()

    return args


def bed_coverage_dict(in_bed):

    cov_dict = {}

    with open(in_bed) as itc:
        for line in itc:
            ls = line.split('\t')
            ch = ls[0].strip(' ')
            pos = ls[1].strip(' ')
            key = ch + ':' + pos
            cov_dict[key] = 0

    return cov_dict


def generate_coverage_dict(alignmentFile, in_bed):

    coverage_dict = bed_coverage_dict(in_bed)

    for read in alignmentFile.fetch(until_eof=True):

        nh = read.get_tag('NH')
        inc = 1.0 / nh
        ch = read.reference_name
        blocks = read.get_blocks()

        for start, end in blocks:
            update_coverage_dict(coverage_dict, ch, start + 1, end + 1, inc)

    return coverage_dict


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


if __name__ == "__main__":

    params = vars(get_args())

    in_dir = params['in_dir']
    #in_gdf = params['in_gdf']
    in_bed = params['in_bed']
    out_dir = params['out_dir']
    save_report = params['s']

    out_directory = '.'
    if out_dir:
        out_directory = out_dir

    if save_report:
        report_file = os.path.join(out_directory, 'ChIP_coverage_files.log')
        logging.basicConfig(filename=report_file, level=logging.INFO)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.getLogger().addHandler(logging.StreamHandler())

    logger.info('Call to ChIP_coverage_files module.')
    logger.info('')
    logger.info('-' * 40)

    #  calculate coverages for bam files
    for file in os.listdir(in_dir):

        if file.endswith('.bam'):

            base_name = file.split('.')[0]
            in_path = os.path.join(in_dir, file)
            out_path = os.path.join(out_dir, base_name + '.coverage.txt')
            out_norm_path = os.path.join(out_dir, base_name + '.norm.txt')

            logger.info('bam file : %s ...' % in_path)

            # Load a bam file and create coverage information for mapped regions

            st = pysam.AlignmentFile(open(in_path, 'rb'))
            coverage_dict = generate_coverage_dict(st, in_bed)

            # Merge coverage information into a new file
            write_coverage_file(coverage_dict, out_path)
            logger.info('Coverage file saved in %s.' % out_path)

            # Load coverage file and normalize
            logger.info('Normalizing coverage file...')
            df = pd.read_csv(out_path, sep='\t', names=range(3))
            total = sum(df[2]) / 50
            # logger.info(total)
            scaling_factor = (total / 1000000)
            # logger.info(scaling_factor)
            df[2] = round(df[2] / scaling_factor, 2)
            df.to_csv(out_norm_path, sep='\t', index=None, header=None)
            logger.info('Normalized coverage file saved in: %s' %out_norm_path)

    logger.info('Finished successfully!')


