#!/usr/bin/env python

import os
import logging
import argparse

import pandas as pd
import pysam

import sys
import time

logger = logging.getLogger(__name__)


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--in_dir', '-d', type=str, help='Directory for alignment files(.bam or .sam).', required=True)
    parser.add_argument(
        '--in_gdf', '-g', type=str, help='Gene annotation table (tab-delimited file).', required=True)
    parser.add_argument(
        '--verbose', '-v', type=bool, default=True, help='Verbose (default: True).')
    parser.add_argument(
        '--out_dir', '-o', type=str, default='.', help='Save results in (default: current working directory).')
    parser.add_argument(
        '--store_report', '-s', action='store_true', help='Save report to file.')
    parser.add_argument(
        '--sample_id', '-f', type=str, default=None, help='sample_id for single sample mode execution (eg. just analyze selected samples instead of everythin in the directory)')
    parser.add_argument(
        '--ignore_files', '-i', type=str, nargs='+', default='.', help='Ignore files present in the directory.')

    args = parser.parse_args()

    return args


def out_file_name(in_file, out_dir=None, ext='', prefix=''):

    base_name = os.path.splitext(os.path.basename(in_file))[0]

    ## default is to use the same directory
    if out_dir is None:
        out_dir = os.path.dirname(in_file)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_file = os.path.join(out_dir, prefix + base_name + ext)

    return out_file


def bam_short_report(in_file, verbose=False):
    """
    loop over bam file alignments to get the total number of reads and the number of unmapped reads.
    Not to confuse with the number of alignments, each read can be aligned to several locations.

    :param in_file:
    :return:
    """

    start_time = time.time()

    # read a file in BAM format, create an AlignmentFile object
    st = pysam.AlignmentFile(in_file, 'rb')

    ## init variables
    total_alignments = []
    unmapped_reads = []

    i=0
    for r in st.fetch(until_eof=True):

        if verbose:
            if i > 0 and i % 100000 == 0:
                sys.stderr.write("{} alignment records processed. {} s\n".format(i, time.time() - start_time))
                sys.stderr.flush()
        i += 1

        if r.is_unmapped:
            unmapped_reads.append(r.query_name)
        total_alignments.append(r.query_name)

    st.close()

    if verbose:
        print('Elapsed Time (bam_short_report):', time.time() - start_time)

    ## by making the set we get the number of reads (instead of read alignments)
    return len(set(total_alignments)), len(set(unmapped_reads)), len(total_alignments)


def split_rna_iters(gdf):
    """
    iterate over features in gene information table (gdf), splitting entries in two (exclusive) categories:
        - 'not rRNA' (mRNA)
        - 'rRNA' (rRNA)
    distinguising between regions associated with ribosomal RNA (rRNA) and everything else.

    generate and return iterables for each respective category.

    Note: 'not RNA' (mRNA) includes:
    ['gene', 'pseudogene', 'ncRNA_gene', 'snoRNA_gene', 'tRNA_gene', 'snRNA_gene']
    """

    ## filter out entries containing Ribosomal RNA (rRNA), to get so called 'mRNA'
    mrna_df = gdf[~gdf['type'].str.contains('rRNA')]
    ## keep entries containing Ribosomal RNA (rRNA)
    rrna_df = gdf[gdf['type'].str.contains('rRNA')]

    assert len(gdf) == len(rrna_df) + len(mrna_df)

    ## generate `mrna` iterable - zip columns: ['chr', 'start', 'end', 'gene-id']
    mrna_iter = list(zip(mrna_df['chr'].tolist(), mrna_df['start'].tolist(), mrna_df['end'].tolist(), mrna_df['gene-id'].tolist()))
    ## generate `rrna` iterable - zip columns: ['chr', 'start', 'end']
    rrna_iter = list(zip(rrna_df['chr'].tolist(), rrna_df['start'].tolist(), rrna_df['end'].tolist()))

    return mrna_iter, rrna_iter


def count_rrna(in_file, rrna_iter, verbose=False):
    """
    loop over ribosomal RNA (rRMA) features present in gtf, to recover reads mapping to those regions.
    return the total number of reads mapped to rRNA regions.

    :param in_file:
    :param rna_iter:
    :return:
    """

    start_time = time.time()

    st = pysam.AlignmentFile(in_file, 'rb')

    ## init variables
    rrna = []

    i=0
    ## loop over ribosomal RNA (rRNA) features present in gtf
    for chro, start, end in rrna_iter:

        if verbose:
            if i > 0 and i % 10 == 0:
                sys.stderr.write("{} feature records processed. {} s\n".format(i, time.time() - start_time))
                sys.stderr.flush()
        i += 1

        st.reset()
        for r in st.fetch(chro, start, end):
            rrna.append(r.query_name)

    st.close()

    if verbose:
        print('Elapsed Time (count_rrna):', time.time() - start_time)

    ## by making the set we get the number of reads (instead of read alignments)
    return len(set(rrna))


def mrna_tagged_file(in_file, rna_iter, rrna_iter, verbose=False, out_dir=''):
    """
    Function including, ammong other operations, the main pre-processing steps:
        - A. Removing ribosomal-RNA (rRNA)
        - B. Ignoring Spliced reads
        - C. Tagging each read alignment with the associated gene-id present in the gdf for that region.

    Given a single bam file and two iterables over the gene information table (gdf), one containing regions
    associated with ribosomal RNA (rRNA) and the other everything else.

    The function will loop over the bam file, in order to:
        - get the total number of reads and the number of unmapped reads, present in bam file.
        Not to confuse with the number of read alignments present in bam file.

    then will iterate over the rrna_iter, containing the regions associated to rRNA:
        - get the total number of reads mapped to rRNA regions.

    and finally will loop over the rna_iter, containing the regions associated to RNA (excluding rRNA):
        - return a new tagged bam file, that:
            - Ignored reads aligned to rRNA
            - Ignored Spliced reads
            - Tagged each read alignment with the associated `gene-id` present in the `gdf` for that region.


    what about htseq-count?
    https://htseq.readthedocs.io/en/release_0.11.1/count.html
    """

    # read a file in BAM format, create an AlignmentFile object
    st = pysam.AlignmentFile(in_file, 'rb')

    # total number of reads, number of unmapped read and the number of read alignments present in bam file.
    n_total_reads, n_unmapped_reads, n_read_alignments = bam_short_report(in_file, verbose=verbose)

    # total number of reads mapped to rRNA regions.
    n_rrna_reads = count_rrna(in_file, rrna_iter, verbose=verbose)

    ## name tmp .bam file where to store tagged read alignments
    tmp_bam_file = out_file_name(in_file, out_dir, '.tmp.bam')
    ## initialize tmp .bam file (from template file)
    tmp_bam = pysam.AlignmentFile(tmp_bam_file, 'wb', template=st)

    start_time = time.time()

    ## initialize variables
    spliced_reads_list = []
    tagged_reads_list = []

    i=0
    ## loop over RNA features present in gdf after excluding Ribosomal RNA (rRNA) features
    for chro, start, end, gene in rna_iter:

        if verbose:
            if i > 0 and i % 10 == 0:
                sys.stderr.write("{} feature records processed. {} s\n".format(i, time.time() - start_time))
                sys.stderr.flush()
        i += 1

        st.reset()
        for r in st.fetch(chro, start, end):
            qn = r.query_name

            ## ignore spliced read alignments - N: Skipped region from the reference
            ## for mRNA-to-genome alignment, an `N` operation represents an intron
            if 'N' in r.cigarstring:
                spliced_reads_list.append(qn)

            else:
                if r.has_tag('RG'):
                    r.set_tag('RG', None)

                ## add tag with 'gene_id'
                r.tags += [('GE', gene)]

                tmp_bam.write(r)
                tagged_reads_list.append(qn)

    tmp_bam.close()

    if verbose:
        print('Elapsed Time (mrna_tagged_file):', time.time() - start_time)

    ## name .tagged.bam file name
    tagged_bam_file = out_file_name(in_file, out_dir, ext='.tagged.bam')

    ## `sort` and `index` .tagged.bam bam file
    pysam.sort("-o", tagged_bam_file, tmp_bam_file)
    pysam.index(tagged_bam_file)

    ## remove
    os.remove(tmp_bam_file)

    ## total number of reads, number of unmapped read and the number of read alignments present new 'tagged_bam' file.
    n_tagged_reads, n_unmapped_tagged, n_tagged_alignments = bam_short_report(tagged_bam_file)

    ## number of spliced read alignments
    n_spliced_reads = len(set(spliced_reads_list))
    #new_total = len(set(tagged_reads_list))

    assert len(set(tagged_reads_list)) == n_tagged_reads

    return tagged_bam_file, n_rrna_reads, n_total_reads, n_spliced_reads, n_tagged_reads


def mrna_tagged_files(in_dir, in_gdf, sample_id=None, ignore_files=".", verbose=False, out_dir=''):
    """
    Main function in `RemrRNA.py`.

    Given an `in_dir` containing a set of bam files and a gene information table (gdf), the function will loop over
    bam files (or get a single bam file, using the `sample_id` variable) applying a set of pre-processing steps to
    read alignments, returning a new set of `tagged` bam files.

    The pre-processing steps include:
        - A. Removing ribosomal-RNA (rRNA)
        - B. Ignoring Spliced reads
        - C. Tagging each read alignment with the associated gene-id present in the gdf for that region.

    :param in_dir:
    :param in_gdf:
    :param sample_id:
    :param ignore_files:
    :param out_dir:
    :return:
    """

    ## read `gdf` (gene information table)
    gene_df = pd.read_csv(in_gdf, sep='\t')

    if not 'gene-id' in gene_df.columns:

        ## rename some columns to fit parastous definitions...
        gene_df = gene_df.rename(columns={"gene_id": "gene-id", "Name":"gene-name", "seqid":"chr", "gene_length": "length"})

        ## select only necessary columns: 'bio_type' column is missing in new GTF
        gene_df = gene_df[['gene-id', 'gene-name', 'chr', 'type', 'start', 'end', 'length', 'category']]


    ## generate iterables for features present in gdf. Create two categories, distinguishing between:
    ##      - regions associated with ribosomal RNA (rRNA)
    ##      - everything else.
    rna_iter, rrna_iter = split_rna_iters(gene_df)

    logger.info('Filtering-out ribosomal-RNA (rRNA), spliced reads and tagging reads with associated `gene_id`s.')
    logger.info('Source directory:\t%s' %in_dir)
    logger.info('-' * 50)

    ## get path to sample files
    list_files = []
    ## now samples '.bam' and '.sam' files are inside individual directories.
    for root, dirs, files in os.walk(in_dir, topdown=True):
        for name in files:
            if name.endswith('.bam') or name.endswith('.sam') and (not name.startswith(ignore_files)):
                #print(os.path.join(root, name))
                list_files.append(os.path.join(root, name))

    ## single sample mode! Filter only for 'sample_id'
    if not isinstance(sample_id, type(None)):

        #list_files = [x for x in list_files if os.path.basename(x).startswith(sample_id)]
        list_files = [x for x in list_files if (os.path.basename(x).split('.')[0] == (sample_id) and os.path.basename(x).endswith('.Aligned.sortedByCoord.out.bam'))]

        #import pdb
        #pdb.set_trace()

        assert len(list_files) == 1 ## for now don't allow empty!

    ## init list to store `tagged_bam_file`s
    tagged_files = []
    ## iterate over samples '.bam' and '.sam' files in `in_dir`
    for in_file in list_files:

        ## get sample_file
        sample_file = os.path.basename(in_file)
        logger.info('Sample name:\t%s' % sample_file)

        sample_id = sample_file.split(".")[0]
        sample_out_dir = os.path.join(out_dir, sample_id)

        tagged_file, rrna, total, spliced, tagged = mrna_tagged_file(in_file, rna_iter, rrna_iter, verbose=verbose, out_dir=sample_out_dir)
        tagged_files.append(tagged_file)

        if param_verbose:
            logger.info('-' * 25)
            logger.info('Total number of reads:\t\t\t%s' % format(total, ','))
            logger.info('Ignored ribosomal-RNA (rRNA) reads:\t\t%s' %format(rrna, ','))
            logger.info('Ignored spliced reads:\t\t%s' %format(spliced, ','))
            logger.info('Total number of read alignments to genes:\t%s' % format(tagged, ','))
            logger.info('Output file:\t%s' % tagged_file)
            logger.info('')
        logger.info('')

    return tagged_files


if __name__ == "__main__":

    params = vars(get_args())

    ## initialize args/params passed to script
    param_in_dir = params['in_dir'] ## Directory with alignment files (.bam or .sam)
    param_in_gdf = params['in_gdf'] ## Gene information table (gdf)
    param_verbose = params['verbose']
    param_out_dir = params['out_dir']  ## Save results in (default: current working directory).'
    param_save_report = params['store_report']
    param_sample_id = params['sample_id']
    param_ignore_files = tuple(params['ignore_files'])

    out_directory = '.'
    if param_out_dir:
        out_directory = param_out_dir

    ## ------------------
    ## Logging Execution
    ## ------------------

    ## create report file
    if param_save_report:
        report_file = os.path.join(out_directory, 'RemRNA.report.txt')
        logging.basicConfig(filename=report_file, level=logging.INFO)

    else:
        logging.basicConfig(level=logging.INFO)

    logging.getLogger().addHandler(logging.StreamHandler())

    #logger.info('Call to RemrRNA module.')
    logger.info('Call to PreProcessRNA module.')
    logger.info('')
    logger.info('This module removes ribosomal-RNA, spliced reads and taggs reads with associated `gene_id`s from given .bam files')
    logger.info('It produces a gene-tagged .bam file as output.')
    logger.info('-' * 40)

    ## Call main function - mrna_tagged_files()
    ## wrapper around mrna_tagged_file() function (loop over samples)
    tagged_files = mrna_tagged_files(param_in_dir, param_in_gdf, sample_id=param_sample_id, ignore_files=param_ignore_files, out_dir=param_out_dir)

    logger.info('Finished successfully!')
