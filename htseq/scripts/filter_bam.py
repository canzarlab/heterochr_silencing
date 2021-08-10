#!/usr/bin/env python

import os
import logging
import argparse

import pandas as pd

import pysam
import HTSeq

import collections

import sys
import time

import pdb

logger = logging.getLogger(__name__)

def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--in_bam', '-f', type=str,
        default=None,
        help='Input bam file for which coverage is to be computed.', required=True)

    # parser.add_argument(
    #     '--in_gtf', '-g', type=str,
    #     help='Reference Annotation File (GTF/GFF).', required=True)

    parser.add_argument(
        '--in_gdf', '-d', type=str,
        help='Gene information table (gdf), containing gene_lengths, etc.', required=True)

    parser.add_argument(
        '--verbose', '-v', type=bool,
        default=True,
        help='Verbose (default: True).')

    parser.add_argument(
        '--out_dir', '-o', type=str,
        default='.',
        help='Save results in (default: current working directory).')

    parser.add_argument(
        '--store_report', '-r',
        action='store_true',
        help='Save report to file.')
        
    parser.add_argument(
        '--pipeline_type', '-p', type=str,
        default='RNA',
        help='Pipeline to be used depending on Sequencing technology, for now `ChIP` or `RNA` (e.g. pA-RNA, RIP, total-RNA ...), this has to do with the `strandedness` of the data.')

    args = parser.parse_args()

    return args


## -------------------
## Auxiliary Functions
## -------------------

def invert_strand(iv):
    iv2 = iv.copy()

    if iv2.strand == "+":
        iv2.strand = "-"

    elif iv2.strand == "-":
        iv2.strand = "+"

    else:
        raise ValueError("Illegal strand")

    return iv2

def out_file_name(in_file, out_dir=None, ext='', prefix=''):
    """
    Auxiliar function to return 'suffixed' and 'prefixed' version
    of an 'in_file' name.
    """

    base_name = os.path.splitext(os.path.basename(in_file))[0]

    ## default is to use the same directory
    if out_dir is None:
        out_dir = os.path.dirname(in_file)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_file = os.path.join(out_dir, prefix + base_name + ext)

    return out_file

## ---------
## Functions
## ---------

def bam_short_report(in_file, verbose=False):
    """
    Loop over `bam` file containing 'read alignments' to get:
     - the total number of 'reads'
     - the number of 'mapped reads'
     - the number of 'unmapped reads'
     - the number of 'read alignments'

    Do NOT confuse the number of 'reads' with the number of 'read alignments' present in bam file.
    (reads can be aligned more than once to multiple locations)

    :param in_file:
    :return:
    """

    start_time = time.time()

    # read a file in BAM format, create an `AlignmentFile` object
    st = pysam.AlignmentFile(in_file, 'rb')

    # init variables
    mapped_reads = []
    unmapped_reads = []
    
    i = 0
    # classify reads as: "mapped" or 
    for r_align in st.fetch(until_eof=True):

        if verbose:
            if i > 0 and i % 100000 == 0:
                sys.stderr.write("{} alignment records processed. {} s\n".format(i, time.time() - start_time))
                sys.stderr.flush()
        i += 1

        if r_align.is_unmapped:
            unmapped_reads.append(r_align.query_name)
        else:
            mapped_reads.append(r_align.query_name)

    st.close()

    if verbose:
        print('Elapsed Time (bam_short_report):', time.time() - start_time, '\n')
    
    # number of 'mapped_reads'
    n_mapped_reads = len(set(mapped_reads))
    # number of 'unmapped_reads'
    n_unmapped_reads = len(set(unmapped_reads))
    assert n_unmapped_reads == len(unmapped_reads) # unmapped reads should be present 1 or 0 times
    # total number of 'read alignments'
    n_read_alignments = len(mapped_reads)
    
    # by making the `set()` we get the "number of reads" (instead of "read alignments")
    return n_mapped_reads + n_unmapped_reads, n_mapped_reads, n_unmapped_reads, n_read_alignments


def count_rrna(in_file, rrna_df, stranded=True, verbose=False):
    """
    Loop over all ribosomal RNA (rRNA) features present in the gdf and recover read alignments mapping to those regions.

    More specifically, create (non-overlaping) disjunct genomic regions from gdf using htseq `GenomicArrayOfSets` and use pysam to fetch()
    read alignments associated to that specific genomic region. 
    => Attention! This step is crutial to avoid double counting of read alignments!
        
    :param in_file:
    :param rrna_df:
    :return:
    """
    
    start_time = time.time()
    
    # init counter
    rrna_counts = collections.Counter()
    seen_alignments = collections.Counter()

    ## ---------------------
    ## Genomic Array of Sets: define non-overlapping genomic features to iterate over
    ## ---------------------
    
    ## Attention! 
    ## => Defining non-overlapping features is crutial, othwerwise we induce over-counting.
    
    # - instantiate `GenomicArrayOfSets` for the `exons`/`gene` features: (`RNA` **stranded** / `DNA` **unstranded**)
    rrna_features = HTSeq.GenomicArrayOfSets("auto", stranded=stranded)

    # loop over 'genomic features' present in 'rrna_df' - only including Ribosomal RNA (rRNA) features
    for ff_row in rrna_df.itertuples():
        
        # each feature has an associated 'gene_id' (by construction of the gdf)
        gene_id = ff_row.gene_id
        
        ## -----------------
        ## Genomic Interval: 
        ## -----------------
        ## => https://htseq.readthedocs.io/en/release_0.11.1/genomic.html
        
        ## (from htseq documentation)
        ## - The start of the interval. Note that all positions should be given and interpreted as 0-based value!
        ## => Our case we need to substract 1, since the annotation starts at 1.
        
        ## - The end of the interval. Following Python convention for ranges, this in one
        ## more than the coordinate of the last base that is considered part of the sequence.
        ## => Our case since we should substract 1, we don't do anything and that should be it.
        feature_iv = HTSeq.GenomicInterval(ff_row.seqid, ff_row.start - 1, ff_row.end, ff_row.strand)
        
        # define `GenomicFeature`
        feature = HTSeq.GenomicFeature(gene_id, ff_row.type, feature_iv)
        # add feature to `GenomicArrayOfSets`
        rrna_features[ feature_iv ] += gene_id

    print('Total number of `genes` in `rrna_features` GenomicArrayOfSets Obj: {}'.format(len(rrna_df)))
    
    ## ---------
    ## BAM file: read alignments
    ## ---------
    
    # read a file in BAM format, create an `AlignmentFile` object
    st = pysam.AlignmentFile(in_file, 'rb')
    
    ii = 0
    # loop over non-overlaping 'genomic features' present in `GenomicArrayOfSets` - from ribosomal RNA (rRNA) features present in gdf
    for iv, gene_ids in rrna_features.steps():
        
        gene_ids = sorted(gene_ids)
        
        # if len(gene_ids) > 1:
        #     import pdb; pdb.set_trace()
        
        # region is not empty - contains at least 1 gene
        if len(gene_ids) > 0:
        
            if verbose:
                if ii > 0 and ii % 10 == 0:
                    sys.stderr.write("{} feature records processed. {} s\n".format(ii, time.time() - start_time))
                    sys.stderr.flush()
            ii += 1
            
            ## -----------
            ## Fetch reads
            ## -----------
            
            feature_strand = iv.strand

            st.reset() # reset pysam.AlignmentFile
            
            ## Attention! fetch returns:
            ##    - reads from both strands
            ##    - reads that 'start' or 'end' in interval: [iv.start,  iv.end]
            ##    => e.g. issues in boundaries (double-counting)
            ## Note: this is NOT important in the rrna case
            ##    => we do not introduce duplicated 'read_alignments' in any bam_file, since we are just counting!
            
            # loop over 'read alignments' present in selected 'genomic region'
            for r_align in st.fetch(iv.chrom, iv.start, iv.end):
                
                # invert strand - due to sequencing the strand is reversed!
                read_strand = '+' if r_align.is_reverse else '-'
                
                # idenfify 'read alignment' -  with 'read_name' + 'start_location'
                r_align_id = ':'.join([r_align.query_name, str(r_align.reference_start)])
                
                # TODO: Is `seen_alignments` necessary?
                # => Yes, it avoids artifacts in boundaries (double-counting) due to splitting regions with non-overlaping 'genomic features'.
                # => Using the `seen_alignments` corrects for that.
                if (read_strand == feature_strand or not stranded):
                    
                    if not seen_alignments[r_align_id]:
                        # init counter that will contain 'gene_id's read alignment overlaps
                        seen_alignments[r_align_id] = collections.Counter()
                        unseen_genes = gene_ids
                        
                    else: # seems like this never happens in rrna case
                        #import pdb; pdb.set_trace()
                        unseen_genes = [kk for kk in gene_ids if kk not in seen_alignments[r_align_id]]  
                        if len(unseen_gene) == 0:
                            continue
                    
                    # only count 'genes' that have not been visited before
                    for gene_id in unseen_genes:
                        # count-it!
                        rrna_counts[ gene_id ] += 1

                
    # close connection pysam.AlignmentFile to .bam file
    st.close()

    ## -----------------
    ## Summarize Counts: rRNA counts DataFrame
    ## -----------------
    
    # - Convert `counter` to DataFrame:
    if len(rrna_counts):
        rrna_counts_df = pd.DataFrame.from_dict(rrna_counts, orient='index').reset_index()
        rrna_counts_df = rrna_counts_df.rename(columns={'index':'gene_id', 0:'count'})
    else:
        # empty
        rrna_counts_df = pd.DataFrame(columns=['gene_id', 'count'])
    #rrna_counts_df.head(20)

    if verbose:
        print('Elapsed Time (count_rrna):', time.time() - start_time, '\n')

    return rrna_counts_df, seen_alignments


def filter_bam(in_bam, mrna_df, stranded=True, read_length=50, verbose=False, out_dir=''):
    """
    Loop over all features present in the gdf that are NOT ribosomal RNA (rRNA) and recover read alignments mapping to those regions.

    More specifically, create (non-overlaping) disjunct genomic regions from gdf using htseq `GenomicArrayOfSets` and use pysam to fetch()
    read alignments associated to that specific genomic region. 
    => Attention! This step is crutial to avoid double counting of read alignments!
    
    The function creates a new 'tagged_bam' file where:
        - read alignments that map to rRNA features have been filtered-out
        => remove rrna regions
        - read alignments that do NOT map to a known genomic feature (present in the `gdf`) have been filtered-out
        => remove intergenic regions
        - individual read alignments are tagged with compatible genomic features

    :param in_file:
    :param mrna_df:
    :return:
    """

    start_time = time.time()
    
    # init counter
    gene_counts = collections.Counter()
    seen_alignments = collections.Counter()
    
    ## ---------------------
    ## Genomic Array of Sets: define non-overlapping genomic features to iterate over
    ## ---------------------
    
    ## Attention! 
    ## => Defining non-overlapping features is crutial, othwerwise we induce over-counting.
    
    # - instantiate `GenomicArrayOfSets` for the `exons`/`gene` features: (`RNA` **stranded** / `DNA` **unstranded**)
    features = HTSeq.GenomicArrayOfSets("auto", stranded=stranded)

    # loop over 'genomic features' present in 'mrna_df' - after excluding Ribosomal RNA (rRNA) features
    for ff_row in mrna_df.itertuples():
        
        # each feature has an associated 'gene_id' (by construction of the gdf)
        gene_id = ff_row.gene_id
        
        ## -----------------
        ## Genomic Interval:
        ## -----------------
        ## => https://htseq.readthedocs.io/en/release_0.11.1/genomic.html
        
        ## (from htseq documentation)
        ## - The start of the interval. Note that all positions should be given and interpreted as 0-based value!
        ## => Our case we need to substract 1, since the annotation starts at 1.
        
        ## - The end of the interval. Following Python convention for ranges, this in one
        ## more than the coordinate of the last base that is considered part of the sequence.
        ## => Our case since we should substract 1, we don't do anything and that should be it.
        feature_iv = HTSeq.GenomicInterval(ff_row.seqid, ff_row.start - 1, ff_row.end, ff_row.strand)
        
        # define `GenomicFeature`
        feature = HTSeq.GenomicFeature(gene_id, ff_row.type, feature_iv)
        ## add feature to `GenomicArrayOfSets`
        features[ feature_iv ] += gene_id

    print('Total number of `genes` in `genomic features` GenomicArrayOfSets Obj: {}'.format(len(mrna_df)))
    
    ## ---------
    ## BAM file: read alignments
    ## ---------
    
    # read a file in BAM format, create an `AlignmentFile` object
    st = pysam.AlignmentFile(in_bam, 'rb')
    
    # name 'tmp' .bam file where to store tagged read alignments
    tmp_bam_file = out_file_name(in_bam, out_dir, '.tmp.bam')
    # initialize `.tmp.bam` file (from template file)
    tmp_bam = pysam.AlignmentFile(tmp_bam_file, 'wb', template=st)
    
    ii = 0
    # loop over non-overlaping 'genomic features' present in `GenomicArrayOfSets` - from 'gdf' after excluding Ribosomal RNA (rRNA) features
    for iv, gene_ids in features.steps():
        
        gene_ids = sorted(gene_ids)
        
        # if len(gene_ids) > 1:
        #     import pdb; pdb.set_trace()
        
        # region is not empty - contains at least 1 gene
        if len(gene_ids) > 0:
        
            if verbose:
                if ii > 0 and ii % 1000 == 0:
                    sys.stderr.write("{} feature records processed. {} s\n".format(ii, time.time() - start_time))
                    sys.stderr.flush()
            ii += 1
            
            ## -----------
            ## Fetch reads
            ## -----------
            
            feature_strand = iv.strand
            
            st.reset() # reset pysam.AlignmentFile
            
            #flag_align = False
            
            ## Attention! fetch returns:
            ##    - reads from both strands
            ##    - reads that 'start' or 'end' in interval: [iv.start,  iv.end]
            ##    => e.g. issues in boundaries (double-counting)
            ## Note: we take care of this issues with 'seen_alignments' dictionary
            ##    => make sure NO repeated 'read_alignment' are introduced in bam_file.

            # start_region = iv.start + read_length
            # start_region = min(start_region, iv.end)
            # end_region = max(start_region, iv.end)

            # loop over 'read alignments' present in selected 'genomic region'
            for r_align in st.fetch(iv.chrom, iv.start, iv.end):
            #for r_align in st.fetch(iv.chrom, start_region, end_region):
                
                # invert strand - due to sequencing the strand is reversed!
                read_strand = '+' if r_align.is_reverse else '-'
                
                # idenfify 'read alignment' -  with 'read_name' + 'start_location'
                r_align_id = ':'.join([r_align.query_name, str(r_align.reference_start)])
                
                # TODO: Is `seen_alignments` necessary?
                # => Yes, it avoids artifacts in boundaries (double-counting) due to splitting regions with non-overlaping 'genomic features'.
                # => Using the `seen_alignments` corrects for that.
                if (read_strand == feature_strand or not stranded):
                    
                    if not seen_alignments[r_align_id]:
                        
                        #flag_align = True
                        
                        # init counter that will contain 'gene_id's read alignment overlaps
                        seen_alignments[r_align_id] = collections.Counter()
                        unseen_genes = gene_ids
                    
                        ## - SAMtags: RG - Read group
                        ## https://samtools.github.io/hts-specs/SAMtags.pdf
                        if r_align.has_tag('RG'):
                            r_align.set_tag('RG', None)
            
                        ## TODO: not perfect, I guess one option would be to do this at the end by looping over `seen_alignments`
                        ## - SAMtags: GE (not defined) -
                        ## add tag with 'gene_id'
                        r_align.tags += [('GE', '-'.join(gene_ids))]
            
                        # write tagged read alignment in .tmp.bam file
                        tmp_bam.write(r_align)
                        
                    else: # seems like this also never happens in mrna case
                        #import pdb; pdb.set_trace()
                        unseen_genes = [kk for kk in gene_ids if kk not in seen_alignments[r_align_id]]  
                        if len(unseen_gene) == 0:
                            continue
                        
                    # only count 'genes' that have not been visited before
                    for gene_id in unseen_genes:
                        # count-it!
                        gene_counts[ gene_id ] += 1
                
            # if (flag_align):
            #     r_aligns = [r_align for r_align in st.fetch(iv.chrom, iv.start, iv.end)]
            #     #r_align = r_aligns[-1]
            #     r_align = r_aligns[0]
            #     print(r_align.to_dict())
            #     import pdb; pdb.set_trace()

    # close connection pysam.AlignmentFile to .bam file
    st.close()
    # close connection pysam.AlignmentFile to .tmp.bam file
    tmp_bam.close()

    ## ---------------
    ## Sort and Index: tagged.bam bam file
    ## ---------------

    # name .tagged.bam file name
    tagged_bam_file = out_file_name(in_bam, out_dir, ext='.tagged.bam')

    # `sort` and `index` .tagged.bam bam file
    pysam.sort("-o", tagged_bam_file, tmp_bam_file)
    pysam.index(tagged_bam_file)

    # remove `.tmp.bam` file
    os.remove(tmp_bam_file)

    ## -----------------
    ## Summarize Counts: mRNA counts DataFrame
    ## -----------------
    
    # - Convert `counter` to DataFrame:
    if len(gene_counts):
        gene_counts_df = pd.DataFrame.from_dict(gene_counts, orient='index').reset_index()
        gene_counts_df = gene_counts_df.rename(columns={'index': 'gene_id', 0: 'count'})
    else:
        # empty
        gene_counts_df = pd.DataFrame(columns=['gene_id', 'count'])
    #gene_counts_df.head(20)

    if verbose:
        print('Elapsed Time (filter_bam):', time.time() - start_time, '\n')

    return gene_counts_df, seen_alignments


def mrna_tagged_file(in_bam, in_gdf, stranded=True, verbose=False, out_dir=''):
    """
    Main function in `filter_bam.py`, new version of the `RemrRNA.py`.

    Given an "input" `bam` file and a gene information table (gdf), the function will loop over read alignments,
    applying a set of pre-processing steps and returning a new "tagged" `bam` file.

    The pre-processing steps include:
        - A. Removing ribosomal-RNA (rRNA)
        - B. Tagging each read alignment with the associated gene-id present in the gdf for that region.

    Given a single bam file, the function will split the gene information table (gdf) into two distinct DataFrames:
        - 'rrna_df': contains regions associated with ribosomal RNA (rRNA)
        - 'mrna_df': contains everything else (poor naming)

    The function will loop over the `bam` file, in order to get:
        - the total number of 'reads'
        - the number of 'unmapped reads' present in the bam file
    Not to be confused with the number of 'read alignments' present in bam file.
    (reads can be aligned more than once to different locations).

    After that, it will iterate over the `rrna_iter`, containing the regions associated to ribosomal-RNA (rRNA):
        - get the total number of reads mapped to rRNA regions.

    Finally, it will loop over the 'rna_iter' containing the regions of interest associated to RNA (excluding rRNA):
        - return a new tagged bam file, that:
            - Ignored reads aligned to rRNA
            - Tagged each read alignment with the associated `gene-id` present in the `gdf` for that region.

    :param in_dir:
    :param in_gdf:
    :param sample_id:
    :param ignore_files:
    :param out_dir:
    :return:
    """

    ## get sample_file
    sample_file = os.path.basename(in_bam)
    logger.info('Sample name:\t%s' % sample_file)

    sample_id = sample_file.split(".")[0]
    sample_out_dir = os.path.join(out_dir, sample_id)

    # -------------------------------------------------------------
    # --------    Import Gene Annotation Data: in_gdf     ---------
    # -------------------------------------------------------------

    ## read `gdf` (gene information table)
    gdf = pd.read_csv(in_gdf, sep='\t')

    if not 'gene_id' in gdf.columns:
        ## rename some columns to fit parastous definitions...
        gdf = gdf.rename(columns={"gene-id": "gene_id", "gene-name": "gene_name", "chr": "seqid", "length": "gene_length"})

    elif not 'gene-name' in gdf.columns:
        ## rename some columns to fit parastous definitions...
        #gene_df = gene_df.rename(columns={"gene_id": "gene-id", "Name": "gene-name", "seqid": "chr", "gene_length": "length"})
        gdf = gdf.rename(columns={"Name": "gene_name"})

    ## select only necessary columns: 'bio_type' column is missing in new GTF
    #gene_df = gene_df[['gene-id', 'gene-name', 'chr', 'type', 'start', 'end', 'length', 'category']]
    #gdf = gdf[['gene_id', 'gene_name', 'seqid', 'type', 'start', 'end', 'strand', 'gene_length', 'category']]
    gdf = gdf[['gene_id', 'gene_name', 'seqid', 'type', 'start', 'end', 'strand']]

    ## --------------------------------------------------
    ## - Create iterators from Gene Annotation Data (gdf)
    ## --------------------------------------------------

    ## Generate iterables for features in gene information table (gdf), split entries in two (disjoint) categories:
    ##     - 'not rRNA' (mRNA)
    ##     - 'rRNA' (rRNA)
    ## distinguish between regions associated with ribosomal RNA (rRNA) and everything else.

    ## Note that, 'not RNA' (mRNA) category also includes:
    ## ['gene', 'pseudogene', 'ncRNA_gene', 'snoRNA_gene', 'tRNA_gene', 'snRNA_gene']

    logger.info('Filtering-out ribosomal-RNA (rRNA) and tagging reads with associated `gene_id`s.')
    logger.info('-' * 50)
    #rna_iter, rrna_iter = split_rna_iters(gene_df)

    # filter out entries containing Ribosomal RNA (rRNA), to get so called 'mRNA'
    #mrna_df = gdf[~gdf['type'].str.contains('rRNA')]
    mrna_df = gdf[~gdf['gene_id'].str.contains('SPRRNA', na=False)]
    
    # keep entries containing Ribosomal RNA (rRNA)
    #rrna_df = gdf[gdf['type'].str.contains('rRNA')]
    rrna_df = gdf[gdf['gene_id'].str.contains('SPRRNA', na=False)]
    
    assert len(gdf) == len(rrna_df) + len(mrna_df)

    # -------------------------------------------------------------
    # ---------    Import read Alignment file: in_bam   -----------
    # -------------------------------------------------------------

    # total number of 'reads', number of 'mapped/unmapped reads' and the number of 'read alignments' present in bam file.
    n_total_reads, n_mapped_reads, n_unmapped_reads, n_read_alignments = bam_short_report(in_bam, verbose=verbose)

    # ------------------------------
    # A. Analyze (rRNA) in bam file:
    # ------------------------------

    # total number of reads mapped to rRNA regions.
    #rrna_counts_df = count_rrna(in_bam, rrna_df, stranded=True, verbose=verbose)
    rrna_counts_df, rrna_alignments_counter = count_rrna(in_bam, rrna_df, stranded=stranded, verbose=verbose)

    # -------------------
    # B. Filter bam file: (everything but rRNA)
    # -------------------

    #gene_counts_df = filter_bam(in_bam, mrna_df, stranded=True, verbose=verbose, out_dir=sample_out_dir)
    gene_counts_df, gene_alignments_counter = filter_bam(in_bam, mrna_df, stranded=stranded, verbose=verbose, out_dir=sample_out_dir)
    
    # -------------------------------------------------------------
    # ------------     Summary: rRNA & other RNA      -------------
    # -------------------------------------------------------------

    logger.info('-' * 25)
    logger.info('Total number of reads:\t\t\t%s' % format(n_total_reads, ','))
    logger.info('Number of read alignments:\t\t\t%s' % format(n_read_alignments, ','))
    logger.info('Number of mapped reads:\t\t\t%s' % format(n_mapped_reads, ','))
    logger.info('Number of unmapped reads:\t\t\t%s' % format(n_unmapped_reads, ','))
    
    logger.info('Ignored ribosomal-RNA (rRNA) reads:\t\t%s' % format(rrna_counts_df['count'].sum(), ','))
    #logger.info('Ignored ribosomal-RNA (rRNA) reads:\t\t%s' % format(rrna_counts_df['read_count'].sum(), ','))

    logger.info('Number of read alignments ( to feature: `genes`):\t%s' % format(gene_counts_df['count'].sum(), ','))
    #logger.info('Number of read alignments ( to feature: `genes`):\t%s' % format(gene_counts_df['read_count'].sum(), ','))

    #logger.info('Output file:\t%s' % tagged_file)

    logger.info('')
    logger.info('')

    # -------------------------------------------------------------
    # --------       Store counts: rRNA & other RNA       ---------
    # -------------------------------------------------------------

    rrna_counts_df.to_csv(os.path.join(sample_out_dir, 'rrna_counts.csv'), sep='\t', index=False)
    gene_counts_df.to_csv(os.path.join(sample_out_dir, 'gene_counts.csv'), sep='\t', index=False)

    return rrna_counts_df, gene_counts_df


if __name__ == "__main__":

    ## -----------
    ## Parameters
    ## -----------

    params = vars(get_args())

    ## initialize args/params passed to script
    param_in_bam = params['in_bam'] ## path to alignment file (.bam or .sam)
    param_in_gdf = params['in_gdf'] ## Gene information table (gdf), containing gene_lengths, etc.
    param_verbose = params['verbose']
    param_out_dir = params['out_dir']  ## Save results in (default: current working directory).'
    param_save_report = params['store_report']
    
    ## Sequencing technology, ChIP, pA-RNA, RIP, total-RNA will determine:
    ##  - How to we deal with the `strandedness` of the data.
    param_pipeline_type =  params['pipeline_type']
    if param_pipeline_type == 'RNA':
        ## Default: True
        stranded = True
        
    elif param_pipeline_type == 'ChIP':
        ## Default: False
        stranded = False
        
    else:
        raise ValueError("Unknown `seq_type`: {} specified.".format(seq_type))

    out_directory = '.'
    if param_out_dir:
        out_directory = param_out_dir

    ## ------------------
    ## Logging Execution
    ## ------------------

    ## create report file
    if param_save_report:
        report_file = os.path.join(out_directory, 'filter_bam.report.txt')
        logging.basicConfig(filename=report_file, level=logging.INFO)

    else:
        logging.basicConfig(level=logging.INFO)

    logging.getLogger().addHandler(logging.StreamHandler())

    # logger.info('Call to RemrRNA module.')
    logger.info('Call to PreProcess_RNA (filter_bam) module.')
    logger.info('')
    logger.info('This module removes ribosomal-RNA, (optionally) spliced-reads and tags reads with associated `gene_id`s from a given .bam file')
    logger.info('The script returns a gene-tagged .bam file as output.')
    logger.info('-' * 40)

    ## Call main function - mrna_tagged_files()
    ## wrapper around mrna_tagged_file() function (loop over samples)
    rrna_df, mrna_df = mrna_tagged_file(param_in_bam, param_in_gdf, stranded=stranded, verbose=param_verbose, out_dir=param_out_dir)
    
    logger.info('Finished successfully!')
