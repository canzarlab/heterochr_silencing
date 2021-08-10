#!/usr/bin/env python

import os
import logging
import argparse

import pandas as pd
import numpy as np
import pysam

import Util

pd.options.mode.chained_assignment = None

header = Util.header


logger = logging.getLogger(__name__)


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--in_dir', '-d', type=str, help='Directory for alignment files(.bam or .sam)', required=True)
    parser.add_argument(
        '--in_gdf', '-g', type=str, help='Gene information table', required=True)
    parser.add_argument(
        '--out_dir', '-o', type=str, default='', help='Save results in (default: current working directory).')
    parser.add_argument(
        '--prefix', '-x', type=str, default='', help='Name prefix for count tables. (default:'')')
    parser.add_argument(
        '--sample_id', '-f', type=str, default=None, help='Sample_id for single sample mode execution (eg. just analyze selected samples instead of everythin in the directory)')
    parser.add_argument(
        '-r', action='store_true')

    args = parser.parse_args()

    return args


def generate_out_path(in_file, out_dir='', prefix=''):

    dir_name = out_dir
    out_file = os.path.join(dir_name, prefix + in_file)

    return out_file


def generate_mrna_iter(gdf):
    """
    generate iterable from gene information table (gdf) after removing entries
    containing ribosomal RNA (rRNA).
    """

    ## filter out entries containing Ribosomal RNA (rRNA).
    mrna_df = gdf[~gdf['type'].str.contains('rRNA')]
    ## generate iterable - zip columns: ['chr', 'start', 'end', 'gene-id']
    mrna_iter = list(zip(mrna_df['chr'].tolist(), mrna_df['start'].tolist(), mrna_df['end'].tolist(), mrna_df['gene-id'].tolist()))

    return mrna_iter


def report_spliced(alignment_file):
    """
    iterate over the `alignment_file` to categorize every read as either: 'spliced'
    or 'not_spliced'. return the nunber of 'spliced' and 'not_spliced' reads.
    """

    ## initialize variables
    spliced = 0
    not_spliced = 0
    total = 0

    ## re-initialize iterator
    alignment_file.reset()

    ## iterate over BAM (Binary Alignment/Mapping File)
    for r in alignment_file.fetch(until_eof=True):

        if 'N' in r.cigarstring:
            spliced += 1
        else:
            not_spliced += 1
        total += 1

    assert total == spliced + not_spliced

    return spliced, not_spliced


def sample_count_column(in_file, mrna_iter, h_df, repeat_ids):
    """
    main function for counting reads.

    BAM file is parsed (fetched) in intervals defined by the `features` found in the gene information table (gdf) using
    the `mrna_iter` iterable.

    Different approach from htseq, where each individual read location is used to interrogate (features[ iv ].steps())
    the reference annotation file (gtf/gff).
    """

    gene_count = []

    # get file's prefix - without path and extensions
    sample = os.path.basename(in_file).split('.')[0]
    # read a file in BAM format, create an AlignmentFile object
    st = pysam.AlignmentFile(in_file, 'rb')

    ## return the number of 'spliced' and 'not_spliced' reads
    spliced, not_spliced = report_spliced(st)

    logger.info('Total number of alignments : %d' % (spliced + not_spliced))
    logger.info('Spliced alignments : %d' % spliced)
    logger.info('Non spliced : %d' % not_spliced)

    # for each `feature` f in the `mrna_iter` (derived from gdf by filtering rRNA)
    for chro, start, end, gene in mrna_iter:

        # reset the AlignmentFile object
        st.reset()
        # list to store `read` alignments to `feature`
        qnames = []

        ## iterate looking for `read`s r matching the `feature` f
        for r in st.fetch(chro, start, end):
            qn = r.query_name

            ## ignore spliced-reads
            if 'N' not in r.cigarstring:
                nh = r.get_tag('NH')
                
                if (gene in repeat_ids and nh <= 15) or (gene not in repeat_ids and nh == 1):
                    # Pre-processing step
                    gene_count.append((qn, gene))

        total = len(qnames)
        gene_count.append((gene, total))

    count_column = []

    df = pd.DataFrame(gene_count, columns=['QN', 'gene-id'])
    df = df.drop_duplicates()

    for gene, group in df.groupby(['gene-id']):
        count_column.append((gene, len(group)))

    column_df = pd.DataFrame(count_column, columns=['gene-id', sample])
    column_df = h_df.merge(column_df, on='gene-id', how='left')
    column = list(column_df[sample])

    return sample, column


def domain_count_table(in_dir, in_gdf, sample_id=None):

    ## read `gdf` (gene information table)
    gene_df = pd.read_csv(in_gdf, sep='\t')

    if not 'gene-id' in gene_df.columns:

        ## rename some columns to fit parastous definitions...
        gene_df = gene_df.rename(columns={"gene_id": "gene-id", "Name":"gene-name", "seqid":"chr", "gene_length": "length"})

        ## select only necessary columns: 'bio_type' column is missing in new GTF
        gene_df = gene_df[['gene-id', 'gene-name', 'chr', 'type', 'start', 'end', 'length', 'category']]


    ## header? - DataFrame with `gene-id`s
    h_df = gene_df[['gene-id']]
    repeat_ids = list(set(gene_df[gene_df['category'] == 'repeat']['gene-id']))
    mrna_iter = generate_mrna_iter(gene_df)

    ## init list of sample names
    new_cols = []
    ## init `count matrix` with an empty Column
    count_mtx = np.empty([len(h_df), 1])

    ## get sample files
    list_files = []
    ## now samples '.bam' and '.sam' files are inside individual directories.
    for root, dirs, files in os.walk(in_dir, topdown=True):
        for name in files:
            #if name.endswith('.bam') or name.endswith('.sam'):
            if name.endswith('.Aligned.sortedByCoord.out.bam') or name.endswith('.Aligned.sortedByCoord.out.tagged.bam') or name.endswith('.sam'):
                #print(os.path.join(root, name))
                list_files.append(os.path.join(root, name))

    ## single sample mode! Filter only for 'sample_id'
    if not isinstance(sample_id, type(None)):

        #list_files = [x for x in list_files if os.path.basename(x).startswith(sample_id)]
        #list_files = [x for x in list_files if (os.path.basename(x).split('.')[0] == (sample_id) and os.path.basename(x).endswith('.Aligned.sortedByCoord.out.bam'))]
        list_files = [x for x in list_files if os.path.basename(x).split('.')[0] == (sample_id)]

        #import pdb
        #pdb.set_trace()

        assert len(list_files) == 1 ## for now don't allow empty!

    ## iterate over samples '.bam' and '.sam' files in `in_dir`
    ## each sample will produce a column that will be appended to the `count_mtx`
    
    for in_file in list_files:

        # select one '.bam' and '.sam' file
        #in_path = os.path.join(in_dir, in_file)
        logger.info('Input bam: %s' % in_file)

        ## return `sample` and respective `count column`
        sample, column = sample_count_column(in_file, mrna_iter, h_df, repeat_ids)
        new_cols.append(sample)
        count_mtx = np.c_[count_mtx, np.array(column)]

    assert np.delete(count_mtx, 0, 1).shape == (len(h_df), len(new_cols))

    ## create dataframe from count matrix - delete empty column used to init the count matrix
    gene_count_df = pd.DataFrame(np.delete(count_mtx, 0, 1), columns=new_cols)
    count_df = pd.concat([h_df, gene_count_df], axis=1)
    count_df = pd.merge(gene_df[header], count_df, on='gene-id')

    ## substitute NA's by 0 - NA's appear from sample's expression not containing a particular gene?
    ## although for even one sample we already see problems!
    count_df = count_df.fillna(0)
    ## where do duplicates come from - ?
    count_df = count_df.drop_duplicates()

    ## count entries corresponding to ribosomal RNA (rRNA)
    rrna_total = len(count_df[count_df['type'].str.contains('rRNA')])
    ## filter out - rRNA entries
    count_df = count_df[~count_df['type'].str.contains('rRNA')]

    logger.info('Removed %d rRNA genes.' % rrna_total)
    logger.info('Calculated raw counts for %d genes.' % len(set(count_df['gene-id'])))

    return count_df


def domain_tpm_table(count_df):
    """
    get TPM-normed Gene Counts Table (tpm_gxt) from raw Gene Counts.

    ----
    Def. Transcripts Per Kilobase Million (TPM)
    ----
    https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

    1. Divide the read counts by the length of each gene in kilobases.
    This gives you reads per kilobase (RPK)

    2. Count up all the RPK values in a sample and divide this number
    by 1,000,000. ('per million' scaling factor).

    3. Divide the RPK values by the 'per million' scaling factor.
    This gives you transcripts per kilobase million (TPM).
    """

    ## make copy
    tpm_df = count_df.copy()

    ## transform gene length to kilobases
    tpm_df['length'] /= 1000

    # ignore header columns:
    # ['gene-id', 'gene-name', 'length', 'type', 'category', 'bio_type']
    # only keep columns containing samples `counts`
    columns = [item for item in tpm_df.columns if item not in header]
    gene_lengths = list(tpm_df['length'])

    ## each column is a `sample`
    for column in columns:

        ## 1. Reads per kilobase (RPK)
        ## Divide the `read counts` by the `length` of each gene in kilobases.
        rpk = [i / float(j) for i, j in list(zip(list(tpm_df[column]), gene_lengths))]

        ## 2. 'per million' scaling factor
        ## Count up all the RPK values in a sample and divide this number by 1,000,000. ().
        per_million = sum(rpk) / 1000000

        ## 3. Transcripts per Million (TPM)
        ## Divide the RPK values by the 'per million' scaling factor.
        tpm = [i / per_million for i in rpk]
        tpm_df[column] = pd.Series(tpm).values

    tpm_df = tpm_df.drop(['length'], axis=1)

    logger.info('Calculated TPM for %d genes.' % len(set(tpm_df['gene-id'])))

    return tpm_df


def generate_raw_tpm_csv(in_dir, in_gdf, out_dir='.', prefix='', sample_id=None):

    ## -------------
    ## Counts Table:
    ## -------------

    count_df = domain_count_table(in_dir, in_gdf, sample_id=sample_id)
    #count_file = generate_out_path(prefix + 'gene_count_table.csv', out_dir)
    if isinstance(sample_id, type(None)):
        count_file = os.path.join(out_dir, prefix + 'pombe_gene_count_matrix.csv')
    else:
        count_file = os.path.join(out_dir, sample_id, prefix + sample_id + '_pombe_gene_count_matrix.csv')
    # store Counts Table
    count_df.to_csv(count_file, sep='\t', index=None)

    ## ----------
    ## TPM Table: (using count_df)
    ## ----------

    tpm_df = domain_tpm_table(count_df)
    #tpm_file = generate_out_path(prefix + 'tpm_table.csv', out_dir)
    if isinstance(sample_id, type(None)):
        tpm_file =  os.path.join(out_dir, prefix + 'pombe_tpm_matrix.csv')
    else:
        tpm_file =  os.path.join(out_dir, sample_id, prefix + sample_id + '_pombe_tpm_matrix.csv')
    ## store TPM table
    tpm_df.to_csv(tpm_file, sep='\t', index=None)

    return count_file, tpm_file


if __name__ == "__main__":

    params = vars(get_args())

    ## initialize args/params passed to script
    param_in_dir = params['in_dir'] ## Directory with alignment files (.bam or .sam)
    param_in_gdf = params['in_gdf'] ## Gene information table (gdf)
    param_out_dir = params['out_dir'] ## Save results in (default: current working directory).
    param_prefix = params['prefix']  ## Name prefix (eg. `chip_' or 'rna_') for count tables. (default:'')
    param_sample_id = params['sample_id']
    param_report = params['r']

    if param_out_dir:
        out_dir = param_out_dir
    else:
        out_dir = param_in_dir

    ## ------------------
    ## Logging Execution
    ## ------------------

    #print(os.getcwd())

    # logging.basicConfig(level=logging.INFO)
    report_file = os.path.join(out_dir, 'GeneExpressionTable.report.txt')
    logging.basicConfig(filename=report_file, level=logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler())

    logger.info('Call to GeneExpressionTable module....')
    logger.info('')
    logger.info('This module calculates gene expression in given alignment file(s).')
    logger.info('Input: folder containing .bam file(s)')
    logger.info('Output: raw and tpm gene counts data (.csv) files')
    logger.info('-' * 40)
    
    ## Call main function - generate_raw_tpm_csv()
    ## wrapper around domain_count_table() and domain_tpm_table() functions in_dir, in_gdf, '.', None
    count_file, tpm_file = generate_raw_tpm_csv(param_in_dir, param_in_gdf, out_dir=out_dir, prefix=param_prefix, sample_id=param_sample_id)

    logger.info('Gene count matrix saved in: %s' % count_file)
    logger.info('TPM expression matrix saved in: %s' % tpm_file)

    logger.info('Finished successfully!')
