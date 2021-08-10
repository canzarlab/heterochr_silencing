#!/usr/bin/env python

import itertools
import time
import os
import sys
import argparse

import HTSeq
import pandas as pd

import collections


## ----------
## Parameters
## ----------

# CIGAR match characters (including alignment match, sequence match, and sequence mismatch)
com = ('M', '=', 'X')


## -------------------
## Auxiliary Functions
## -------------------

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--in_bam', '-f', type=str,
        default=None,
        #help='Input bam file for which gene counts are to be computed.', required=True)
        help='Input bam file for which gene counts are to be computed.') # can pass counts directly!
    
    ## (Optional): INPUT subtraction
    parser.add_argument(
        '--in_chip_input', '-i',
        type=str, default=None,
        help='Input ChIP INPUT bam file which counts scaled by a norm factor needs to be subtracted from ChIP counts.')
        
    ## (Optional): INPUT subtraction - Input factor
    parser.add_argument(
        '--norm_factor', '-n',
        type=str, default=None,
        help='Norm factor used to scale the ChIP INPUT coverage before subtraction from coverage.')
        
    ## (Optional): INPUT subtraction - ChIP counts
    # => use counts if have already been computed
    parser.add_argument(
        '--in_count',
        type=str, default=None,
        help='Count file for ChIP data during INPUT subtraction, instead of re-computing.') 
        
    ## (Optional): INPUT subtraction - INPUT counts
    # => use counts if have already been computed
    parser.add_argument(
        '--in_chip_input_count',
        type=str, default=None,
        help='Count file for ChIP INPUT data, instead of re-computing. INPUT counts scaled by a norm factor need to be subtracted from ChIP counts.')

    parser.add_argument(
        '--in_gtf', '-g', type=str,
        #help='Reference Annotation File (GTF/GFF).', required=True) ## not anymore, with gdf might be enough!
        help='Reference Annotation File (GTF/GFF).', required=False)

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
        '--prefix', '-x', type=str,
        default='',
        help='Name prefix for count tables (default:'').')

    parser.add_argument(
        '--seq_type', '-s', type=str,
        default='ChIP',
        help='Sequencing technology, for now `ChIP` or `RNA` (e.g. pA-RNA, RIP, total-RNA ...), this has to do with the `strandedness` of the data.')

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

def invert_strand(iv):
    iv2 = iv.copy()

    if iv2.strand == "+":
        iv2.strand = "-"

    elif iv2.strand == "-":
        iv2.strand = "+"

    else:
        raise ValueError("Illegal strand")

    return iv2

def seq_param_config(seq_type, stranded=None, feature_type="default", length='length'):
    """
    There are default configurations for each type of Sequencing technology: ChIP, pA-RNA, RIP, total-RNA
        - How to we deal with the `strandedness` of the data, need to distinguish between:
            * `ChIP` - DNA
            * `RNA`
        - Which features to use for counting:
            * `genes`: ChIP, (RIP)
            * `exons`: pA-RNA, total-RNA
        - Which length use to Normalize:
            * `gene_length`: ChIP, RIP
            * `transcript_length`: pA-RNA, total-RNA
    """
    config = {}
    
    ## Used for selecting `length` col
    if isinstance(feature_type, type(list())):
        feature_name = feature_type[0]
    
    if seq_type in ['RNA', 'pA-RNA', 'total-RNA', 'simulated-data']:
        ## Default: True
        config['stranded'] = True if isinstance(stranded, type(None)) else stranded
        ## Default: `exon`
        config['feature_type'] = "gene" if feature_type == "default" else feature_type
        #config['feature_type'] = "exon" if feature_type == "default" else feature_type
        #config['feature_type'] = ["CDS", "UTR"] if feature_type == "default" else feature_type
        ## Default: `transcript_length`
        config['length'] = 'transcript_length' if feature_type == "default" else feature_name + '_length'
        
    elif seq_type in ['RIP', 'S2-RIP', 'S5-RIP']:
        ## Default: True
        config['stranded'] = True if isinstance(stranded, type(None)) else stranded
        ## Default: `gene`
        config['feature_type'] = "gene" if feature_type == "default" else feature_type
        ## Default: `gene_length`
        config['length'] = 'gene_length' if feature_type == "default" else feature_name + '_length'
        
    elif seq_type in ['ChIP', 'S2-ChIP', 'S5-ChIP', 'S2-ChIP-INPUT', 'S2-ChIP-OIN', 'H3K9me2']:
        ## Default: False
        config['stranded'] = False if isinstance(stranded, type(None)) else stranded
        ## Default: `gene`
        config['feature_type'] = "gene" if feature_type == "default" else feature_type
        ## Default: `gene_length`
        config['length'] = 'gene_length' if feature_type == "default" else feature_name + '_length'
        
    else:
        raise ValueError("Unknown `seq_type`: {} specified.".format(seq_type))
        
    return config
    

## ---------------------
## Reference Annotation:  GFF - feature Parsers
## ---------------------

def parse_gff_feature_cds(features, gff_file, genes_dict=None, include_utrs=False):

    if isinstance(genes_dict, type(None)):
        genes_dict = {}

    ## loop over all features in gtf file
    for feature in gff_file:

        ## store all exons in our `GenomicArrayOfSets`
        if feature.type == "CDS":

            ## identify each `CDS` feature by parent transcript/gene
            #gene_id = feature.attr["Parent"].split(':')[1][:-2]  ## get rid of '.1'
            gene_id = feature.attr["Parent"][:-2]  ## get rid of '.1'
            features[ feature.iv ] += gene_id

            # Is this the first time we see this gene?
            if gene_id not in genes_dict:
                # If so, add to the 'genes_dict' an empty list
                genes_dict[gene_id] = list()

            # add the feature to the gene list
            genes_dict[gene_id].append(feature)

        elif (feature.type == 'three_prime_UTR' or feature.type == 'five_prime_UTR') and include_utrs:

            ## identify each `*_UTR` feature by parent transcript/gene
            #gene_id = feature.attr["Parent"].split(':')[1][:-2]  ## get rid of '.1'
            gene_id = feature.attr["Parent"][:-2]  ## get rid of '.1'
            features[ feature.iv ] += gene_id

            # Is this the first time we see this gene?
            if gene_id not in genes_dict:
                # If so, add to the 'genes_dict' an empty list
                genes_dict[gene_id] = list()

            # add the feature to the gene list
            genes_dict[gene_id].append(feature)

    return features, genes_dict

def parse_gff_feature_intron(features, gff_file, genes_dict=None):

    if isinstance(genes_dict, type(None)):
        genes_dict = {}

    ## loop over all features in gtf file
    for feature in gff_file:

        ## store all exons in our `GenomicArrayOfSets`
        if feature.type == "intron":

            ## identify each `intron` feature by parent transcript/gene
            gene_id = feature.attr["Parent"][:-2]  ## get rid of '.1'
            features[ feature.iv ] += gene_id

            # Is this the first time we see this gene?
            if gene_id not in genes_dict:
                # If so, add to the 'genes_dict' an empty list
                genes_dict[gene_id] = list()

            # add the feature to the gene list
            genes_dict[gene_id].append(feature)

    return features, genes_dict

def parse_gff_feature_exon(features, gff_file, genes_dict=None):

    if isinstance(genes_dict, type(None)):
        genes_dict = {}

    ## loop over all features in gtf file
    for feature in gff_file:

        ## store all exons in our `GenomicArrayOfSets`
        if feature.type == "exon":

            ## identify each `exon` feature by parent transcript/gene
            gene_id = feature.attr["Parent"].split(':')[1][:-2]  ## get rid of '.1'
            features[feature.iv] += gene_id

            # Is this the first time we see this gene?
            if gene_id not in genes_dict:
                # If so, add to the 'genes_dict' an empty list
                genes_dict[gene_id] = list()

            # add the feature to the gene list
            genes_dict[gene_id].append(feature)

    return features, genes_dict

def parse_gff_feature_gene(features, gff_file, genes_dict=None):

    if isinstance(genes_dict, type(None)):
        genes_dict = {}

    ## loop over all features in gtf file
    for feature in gff_file:

        ## parse features contained in `features_of_interest`
        # if feature.type in features_of_interest:
        if 'gene' == feature.type:

            # get `gene` feature id
            try:
                ## identify each `gene` feature by `gene_id` attribute: transcript/gene
                #gene_id = feature.attr["gene_id"]
                gene_id = feature.attr["ID"]

            except:

                ## sub-set of pseudogenes that behave as transcripts
                assert feature.type == 'pseudogene'
                gene_id = feature.attr["Parent"].split(':')[1]

            ## add `gene` feature to `GenomicArrayOfSets`
            features[feature.iv] += gene_id

            # Is this the first time we see this gene?
            if gene_id not in genes_dict:
                # If so, add to the 'genes_dict' an empty list
                genes_dict[gene_id] = list()

            # add the feature to the gene list
            genes_dict[gene_id].append(feature)

    return features, genes_dict

def parse_gff_feature(features, gff_file, feature_type = ['region'], features_dict=None):

    if isinstance(features_dict, type(None)):
        features_dict = {}

    if isinstance(feature_type, type('chr')):
        feature_type = [feature_type]

    ## loop over all features in gtf file
    for feature in gff_file:

        ## parse features contained in `feature_type`
        if feature.type in feature_type:
        #if feature_type == feature.type:

            ## identify each `feature` by `ID` attribute:
            try:
                feature_id = feature.attr["ID"]

            except:

                raise ValueError("Unknown `feature_type`: {} specified.".format(feature_type))

            ## add feature to `GenomicArrayOfSets`
            features[feature.iv] += feature_id

            # Is this the first time we see this feature?
            if feature_id not in features_dict:
                # If so, add to the 'feature_dict' an empty list
                features_dict[feature_id] = list()

            # add the feature to the gene list
            features_dict[feature_id].append(feature)

    return features, features_dict

def parse_gff(features, gff_file, feature_type = 'gene'):
    """
    Note that multiple 'feature_type's can be passed.

    It is the job of the user to make sure that the are no undesired interactions
    between features.
    """

    genes_dict = {}

    ## multiple features can be parsed!
    if isinstance(feature_type, type('chr')):
        feature_type = [feature_type]

    # - Parse the `gff `file, adding features of interest to the `GenomicArrayOfSets
    known_features = ['exon', 'cds + utrs', 'cds', 'intron', 'gene']

    if 'cds + utrs' in feature_type:
        features, genes_dict = parse_gff_feature_cds(features, gff_file, include_utrs = True, genes_dict = genes_dict)

    if 'cds' in feature_type:
        features, genes_dict = parse_gff_feature_cds(features, gff_file, genes_dict = genes_dict)

    if 'intron' in feature_type:
        features, genes_dict = parse_gff_feature_intron(features, gff_file, genes_dict = genes_dict)

    if 'exon' in feature_type:
        features, genes_dict = parse_gff_feature_exon(features, gff_file, genes_dict = genes_dict)

    if 'gene' in feature_type:
        features, genes_dict = parse_gff_feature_gene(features, gff_file, genes_dict = genes_dict)

    # General purpose parsing for features might not work!
    unknown_features = list(set(feature_type).difference(known_features))
    if len(unknown_features) > 0:
        print("Warning! Computing `feature_type`: {} specified using general purpose parse_gff_feature() function. \nMight behave unexpectedly!\n".format(unknown_features))
        features, genes_dict =  parse_gff_feature(features, gff_file, feature_type = unknown_features, features_dict = genes_dict)
        #raise ValueError("Unknown `feature_type`: {} specified.".format(feature_type))

    return features, genes_dict


## ---------------------
## Reference Annotation:  `GenomicArrayOfSets` &  `genes_dict`
## ---------------------

def get_features_from_gff(in_gff, stranded=True, feature_type='exon'):

    # - Select gtf file: instantiate GFF_Reader Object
    gff_file = HTSeq.GFF_Reader(in_gff)

    # - Obtain `GenomicArrayOfSets` for the features in the `GFF` File (for`RNA` stranded/ `DNA` unstranded)
    # instantiate `GenomicArrayOfSets` for the `exons`/`gene` features: 
    features = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    
    # - Parse the `gff `file, adding features of interest to the `GenomicArrayOfSets`
    features, genes_dict = parse_gff(features, gff_file, feature_type = feature_type)

    print('Total number of `genes` in `features` GenomicArrayOfSets Obj: {}'.format(len(genes_dict)))

    return features, genes_dict

def get_features_from_gdf(gdf, stranded=True):
    """
    Most features we are looking at (e.g. gene's, exons, introns, etc...) are what I called `gene features`.
    These features are part of a hierarchical framework where each feature can be associated uniquely to a `gene_id`.
        
    Other features that we might be interested in like `region`, `repeats`, etc...
    are NOT `gene features` and therefore do NOT conform to the defintion above.
    
    - Important!!!
    => For those features (non `gene features`), when creating the `gdf` I take care of creating a `gene_id` 
    column using their by `ID` attribute. 
    
    """

    # - instantiate `GenomicArrayOfSets` for the `exons`/`gene` features: (`RNA` **stranded** / `DNA` **unstranded**)
    features = HTSeq.GenomicArrayOfSets("auto", stranded=stranded)

    # - Iterate over the `gdf` adding features of interest to the `GenomicArrayOfSets`
    genes_dict = {}

    ## loop over all features in gtf file
    for ff_row in gdf.itertuples():
        
        ## each feature has an associated 'gene_id' (by construction of the gdf)
        gene_id = ff_row.gene_id
        
        ## feature id
        feature_id = ff_row.ID
        
        ## ------------------------
        ## feature genomic interval: https://htseq.readthedocs.io/en/release_0.11.1/genomic.html
        ## ------------------------
        ## - The start of the interval. Note that all positions should be given and  interpreted as 0-based value!
        ## => Our case we need to substract 1, since the annotation starts at 1.
        ## - The end of the interval. Following Python convention for ranges, this in one
        ## more than the coordinate of the last base that is considered part of the sequence.
        ## => Our case since we should substract 1, we don't do anything and that should be it.
        feature_iv = HTSeq.GenomicInterval(ff_row.seqid, ff_row.start - 1, ff_row.end, ff_row.strand)
        if (feature_iv.end == feature_iv.start):
            import pdb; pdb.set_trace()
        feature = HTSeq.GenomicFeature(feature_id, ff_row.type, feature_iv)
        
        ## add feature to `GenomicArrayOfSets
        features[feature_iv] += gene_id

        # Is this the first time we see this feature?
        if gene_id not in genes_dict:
            # If so, add to the 'feature_dict' an empty list
            genes_dict[gene_id] = list()

        # add the feature to the gene list
        genes_dict[gene_id].append(feature)

    print('Total number of `genes` in `features` GenomicArrayOfSets Obj: {}'.format(len(genes_dict)))

    return features, genes_dict
    

## --------------------------
## Read Alignments File - BAM
## --------------------------

def summarize_reads_dict_counter(item, feature_id, debug = False):
    
    if debug:
    
        import numpy as np
        
        col_names = ["read_name", "align_frac", "n_features", "strand"]
        
        # creating the dataframe 
        df = pd.DataFrame(data = np.vstack(item),columns = col_names)
        assert(len(df['strand'].unique()) == 1)
        
        # Parsing cols - simulated data
        df['original_strand'] = [xx.split(':')[1] for xx in df['read_name']]
        df['chrm'] = [xx.split(':')[2] for xx in df['read_name']]
        df['start'] = [xx.split(':')[3].split('-')[1] for xx in df['read_name']]
        
        # add 'gene_id'
        df['gene_id'] = feature_id
        
        df = df.sort_values(by=['chrm', 'start', 'original_strand'])
        
        # if (feature_id == 'FP565355_region_15417..15473'):
        #     import pdb; pdb.set_trace()
        import pdb; pdb.set_trace()

    # [aln.read.name, read_align_frac, n_features]
    #read_count = sum([1 / (xx[1] * xx[2]) for xx in set(item)])
    alignment_count = sum([1 / (xx[1] * xx[2]) for xx in item])
    
    #return read_count, alignment_count
    return alignment_count
    

# ![img_overlap](https://htseq.readthedocs.io/en/release_0.11.1/_images/count_modes.png)
def get_counts(in_bam, features, repeat_ids, filter_region=None, count_mode="union", ambiguous_assignment_mode='none', multimapped_mode='none', max_nh=16):
    """

    """

    ## Init Counters: What's the difference?
    counts = collections.Counter( )
    #reads_dict = collections.defaultdict(list)

    start_time = time.time()

    ## ---------
    ## BAM file:
    ## ---------

    # select bam file: instantiate BAM_Reader Object (RNA-seq)
    bam_file = HTSeq.BAM_Reader(in_bam)
    # variable to store max value of "NH" field in 'bam_file'
    bam_max_nh = 1
    
    if not isinstance(filter_region, type(None)):
        iv = filter_region[0].iv
        #bam_file.reset()
        bam_file = bam_file.fetch(iv.chrom, iv.start, iv.end)
        
    # - Parse each read alignment in `BAM` File and count/assign them to `feature`s
    i = 0
    #for aln in itertools.islice(bam_file, 10000): ## debugging
    for aln in bam_file:

        if i > 0 and i % 100000 == 0:
            sys.stderr.write("{} alignment records processed. {} s\n".format(i,  time.time() - start_time))
            sys.stderr.flush()
        i += 1

        ## --------------------------------------------------------------------
        ## ------------      Inspect read alignment         -------------------
        ## --------------------------------------------------------------------

        ## both _mapped or _unmapped (our BAM files only contain _mapped)
        counts["_total_read_alignments"] += 1

        ## our bam files only contain aligned reads
        if not aln.aligned:
            counts["_unmapped_read_alignments"] += 1
            continue  # skips to next iteration

        ## just in case, don't know if its possible
        #assert aln.optional_field("NH") != 0

        ## ---------------------------------------------------------------------
        ## ----  1. Reads mapped to multiple locations (multimapped_mode)  -----
        ## ---------------------------------------------------------------------

        ## Multimapped reads are represented as separate entries in the BAM file.

        ## We will potentially allow for multimapped-read alignments (e.g. read alignments
        ## that originate from the same read and have been mapped to multiple locations
        ## by the aligner), but only if the read aligns to less than `max_nh` locations.

        ## - How to deal with read alignments that originate from a read mapped to multiple locations?
        ##      A. `none`: ignore all multimapped reads.
        ##      B. `fractional`: allow them but assign multimapped reads a fractional count.
        ##      C. `ignore_secondary`: allow them but only count them once, use the primary_alignment.

        ##  ("none", "all", "fractional", "ignore_secondary")
        
        nh_field = aln.optional_field("NH")

        ## init `multimapped_flag`
        multimapped_flag = nh_field > 1
        ## default is to count 'all' read alignments in full, including multimaps (read_align_frac = 1)
        read_align_frac = 1

        ## find max value of NH in 'bam_file'
        if nh_field > bam_max_nh:
            bam_max_nh = nh_field

        ## Potentially allow multimapped reads - with less than max_nh (16) hits
        ## - We only allow this if the read maps to a repeat feature.
        if (nh_field > 1) and (nh_field <= max_nh):

            counts["_read_alignments_not_unique"] += 1

            ## ----------
            ## A. `none`: ignore all multimapped reads
            ## ----------
            if multimapped_mode == 'none':
                continue # skips to next iteration

            ## ----------------
            ## B. `fractional`: count multimapped reads fractionally
            ## ----------------
            ## - Assign fractional weight-count (~NH): corresponding to number of times the read was aligned (NH)
            elif multimapped_mode == 'fractional':
                read_align_frac = nh_field

            ## ----------------------
            ## C. `ignore_secondary`: only count the primary alignment of multimapped reads
            ## ----------------------
            # - C.I. `ignore_secondary`: ignore secondary_alignment reads ....
            elif (multimapped_mode == 'ignore_secondary') and (aln.not_primary_alignment):
                counts["_not_primary_alignment"] += 1
                continue # skips to next iteration

            ## C.II. ignore_secondary`: let pass `primary_alignment` reads
            ##  - From now own these read alignments behave as a unique alignments
            elif (multimapped_mode == 'ignore_secondary') and (not aln.not_primary_alignment):
                multimapped_flag = False # invert flag - treat as unique alignment from now on
                #counts["_alignment_not_unique"] -= 1

        ## - Multimapped reads - with more than max_nh (16) hits will always be ignored!
        elif nh_field > max_nh:
            counts["_read_alignments_max_NH_" + str(max_nh)] += 1
            continue # skips to next iteration

        #minaqual = 10
        #if aln.aQual < minaqual:
        #    #import pdb; pdb.set_trace()
        #    counts["_too_low_aQual"] +=  1
        #    continue


        ## ---------------------------------------------------------------------
        ## -----    2. `Read` and `Feature` overlap mode (count_mode)    -------
        ## ---------------------------------------------------------------------

        ## Given a single `read alignment` that overlaps (partially or fully) with a particular
        ## feature (or set of features), there are several ways to define consistency
        ## between the `read alignment` and the feature/features it overlaps.

        ## - How to deal with overlap of `read alignment` and `feature`?
        ## => How to decide if a `read alignment` is consistent with a `feature`?
        ##      A. `union`: the union of all the sets S(i). This mode is recommended for most use cases.
        ##      B. `intersection-strict`: the intersection of all the sets S(i).
        ##      C. `intersection-nonempty`: the intersection of all non-empty sets S(i).

        ## If or how to assign a count to the resulting features in ambiguous cases is controled by 
        ## another paramater (or mode): 
        ## => `ambiguous_assignment_mode`.

        ## These definitions are taken from htseq and can be useful to see the summary table:
        ## https://htseq.readthedocs.io/en/master/count.html#count
        
        ## The result of this step should be a `feature_set` associated with the corresponding `read alignment`.
        ## Note that: if there is no associated feature `feture_set` will be and empty set or 'None'
        
        ## -------------
        ## Invert Strand
        ## -------------
        
        ## invert strand - due to sequencing the strand is reversed!
        iv_seq = (invert_strand(co.ref_iv) for co in aln.cigar if (co.type in com and co.size > 0))
        #iv_seq = (co.ref_iv for co in aln.cigar if (co.type in com and co.size > 0))

        ## In the following step, we loop over CIGAR operations (cig_op) from the aligned read using
        ## the corresponding 'count_mode' scheme.

        ## ---------
        ## A. Union: the union of all the sets S(i). This mode is recommended for most use cases.
        ## ---------

        if count_mode == "union":

            ## feature set
            gene_ids = set()

            for iv in iv_seq:
                #if iv.chrom not in exon_features.chrom_vectors:
                #    raise UnknownChrom
                for iv2, fs in features[ iv ].steps():
                    gene_ids = gene_ids.union(fs)

        ## ----------------------
        ## B. Intersection-strict: the intersection of all the sets S(i).
        ## ----------------------

        ## -------------------------
        ## C. Intersection-nonempty: the intersection of all non-empty sets S(i).
        ## -------------------------

        elif count_mode in ("intersection-strict", "intersection-nonempty"):

            ## feature set
            gene_ids = None

            for iv in iv_seq:
                #if iv.chrom not in exon_features.chrom_vectors:
                #    raise UnknownChrom
                for iv2, fs in features[ iv ].steps():
                    if ((len(fs) > 0) or (count_mode == "intersection-strict")):
                        if gene_ids is None:
                            gene_ids = fs.copy()
                        else:
                            gene_ids = gene_ids.intersection(fs)

        ## Other: Ilegal!
        else:
            sys.exit("Illegal overlap mode.")


        ## ---------------------------------------------------------------------
        ## -----  3. Weight overlapping features (ambiguous_assignment_mode) ---
        ## ---------------------------------------------------------------------

        ## - How should we assign counts to each pair `read alignment`/`feature`?
        ##      I. Read aligned to unknown feature: (not ambigous)
        ##          - feature set is empty/None.

        ##      II. Read aligned to a unique feature: (not ambigous)
        ##          - feature set contains one unique element.

        ## - Specially how should we do this for ambigous cases?
        ##      III. Read aligned to a region with overlapping features:
        ##          - feature set contains more than one element!

        ## Case III (ambigous) is controled by yet another paramater (or mode): 
        ## => `ambiguous_assignment_mode`
        ##      A. 'none': ignore all reads mapped to a region with overlapping features.
        ##      B. 'fractional': assign a fractional weight-count to each overlaped feature.
        ##      C. 'all': assign a full count to each overlaped feature.
        
        ## The result of this step should be a second scaling variable: `n_features`

        ## ----------------------------
        ## I. Mapped to unknown feature: (the interval is empty)
        ## ----------------------------

        if (gene_ids is None) or (len(gene_ids) == 0):
            counts["_no_feature_alignment"] += 1
            #count_flag = False
            continue # skips to next iteration

        ## -----------------------------
        ## II. Mapped to a unique feature:  (contains one unique element)
        ## -----------------------------

        elif len(gene_ids) == 1:
            #count_flag = True
            n_features = 1

        ## --------------------------------------------------
        ## III. Mapped to a region with overlapping features: (contains more than one element)
        ## --------------------------------------------------

        ## See next how to deal with this ambiguous read alignments!
        else:
            #count_flag = True
            counts["_ambiguous_alignment"] += 1

            ## -----------------------------------------------------------------------------
            ## - How to count (ambigous) reads mapped to a region with overlapping features?
            ## -----------------------------------------------------------------------------

            ## A. 'none': ignore all reads mapped to a region with overlapping features.
            if ambiguous_assignment_mode == 'none':
                #count_flag = False # don't count
                continue  # skips to next iteration

            ## B. 'fractional': assign a fractional weight-count to each overlaped feature.
            elif ambiguous_assignment_mode == 'fractional':
                n_features = len(gene_ids)

            ## C. 'all': assign a full count to each overlaped feature.
            elif ambiguous_assignment_mode == 'all':
                n_features = 1

            else:
                sys.exit("Illegal ambiguous_assignment mode.")

        ## ---------------------------------------------------------------------
        ## ---------     4. Counting Features (using weights)       ------------
        ## ---------------------------------------------------------------------

        ## -------------------
        ## How to count reads?
        ## -------------------
        
        ## Here is where the actual counting takes place.
        for fsi in list(gene_ids):
            
            # only allow multimapped reads for 'repeat' features
            if (fsi in repeat_ids) or (not multimapped_flag):
                
                counts[ fsi ] += 1 / (read_align_frac * n_features)
                
                # store for each 'feature' the associated 'read alignment's
                #reads_dict[ fsi ].append((aln.read.name, read_align_frac, n_features, aln.iv.strand))

            else:
                counts["_ignore_multimapped_gene"] += 1

    print('Elapsed Time (Counting reads):', time.time() - start_time)

    #{k: v for k, v in sorted(counts.items(), key=lambda item: item[1], reverse=True)}

    ## - Add max_nh value:
    counts["_bam_max_nh"] = bam_max_nh

    ## --------------------
    ## Counts as DataFrames
    ## --------------------
    
    # - Convert `counter` to DataFrame:
    counts_df = pd.DataFrame.from_dict(counts, orient='index').reset_index()
    counts_df = counts_df.rename(columns={'index': 'gene_id', 0: 'count'})
    #counts_df.head(20)
    
    # # - Convert `counter` to DataFrame: (This is the latest version)
    # gene_counts_df = {k: summarize_reads_dict_counter(v, feature_id=k, debug=True) for (k, v) in reads_dict.items()}
    # gene_counts_df = pd.DataFrame.from_dict(gene_counts_df, orient='index').reset_index()
    # gene_counts_df = gene_counts_df.rename(columns={'index': 'gene_id', 0: 'count', 1: 'read_count'})
    
    #return counts_df, gene_counts_df
    return counts_df


## ------------------
## Gene Counts Matrix
## ------------------

def get_tpm_gene_counts(count_df, count_col = 'count', length_col= 'length', read_length = 50):
    """
    Get TPM-normed Gene Counts Table (tpm_gxt) from raw Gene Counts.

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

    tpm_df = count_df.copy()
    
    ## Transform gene lengths - into kilobases (kb)
    #tpm_df[length_col] /= 1000
    ## Transform "effective" gene lengths - into kilobases (kb)
    tpm_df["norm_length"] = (tpm_df[length_col] - (read_length - 1)) / (10 ** 3)

    ## --------------------------
    ## 1. Reads per kilobase (RPK)
    ## --------------------------
    
    # no negative counts!
    assert not (tpm_df[count_col] < 0).any()
    if (tpm_df[count_col] < 0).any():
        tpm_df[tpm_df[count_col] < 0]
        import pdb; pdb.set_trace()
    

    ## Divide the `read counts` by the `length` of each gene in kilobases.
    rpk = tpm_df[count_col] / tpm_df["norm_length"]
    
    ## -------------------------------
    ## 2. 'per million' scaling factor
    ## -------------------------------

    ## Sum up all the RPK values in a sample and divide this number by 1,000,000.
    #per_million = rpk.sum() / 1000000
    per_million = rpk.sum() / (10 ** 6)

    ## --------------------------------
    ## 3. Transcripts per Million (TPM)
    ## --------------------------------

    ## Divide the RPK values by the 'per million' scaling factor.
    tpm_df[count_col] = rpk / per_million

    print('\nCalculated TPM for {} genes.'.format(len(count_df['gene_id'].unique())))
    print('Used `length_col`: {} for normalization.'.format(length_col))

    return tpm_df



def raw_gene_counts_csv(in_bam, in_gdf,
                        in_gff=None,
                        stranded=None,
                        feature_type="gene",
                        count_mode="union", ambiguous_assignment_mode='none', multimapped_mode='none', max_nh=16):
    """
    Main function in `gene_counts.py`, new version of the `gene_expression_table.py`.

    Acts as a wrapper around the get_counts() function where most of the computation takes place.

    Given an `in_bam` file, a gene information table (gdf) and a list of `repeats`
    (special features that will be treated differently), this function will count
    the number of read alignments (from `in_bam` file) consistent with the features present in the gdf.

    Note that, the actual genome annotation (.gff file) is not needed.
    The gdf already contains all the information we need from the gff and additional information like
    feature length, necessary for proper normalization.
    
    There are several parameters that will determine the exact configuration of the execution.
    
    There are default configurations for each type of Sequencing technology: ChIP, pA-RNA, RIP, total-RNA
        - How to we deal with the `strandedness` of the data, need to distinguish between `ChIP` and `RNA`.
        - Which features to use for counting:
            * `genes`: ChIP, RIP
            * `exons`: pA-RNA, total-RNA

    Additionally, 3 parameters will define the way we want to count, using: get_counts()

        A. 'count_mode': How to deal with overlap of READ and FEATURE?
        Options: ("union", "intersection-strict", "intersection-nonempty")

        B. 'ambiguous_assignment_mode': How to deal with reads multiple-overlapping FEATURES? (according to count_mode)
        Options: ("none", "all", "fractional")

        C. 'multimapped_mode': How to deal with multi-mapped reads? including secondary alignments
        Options: ("none", "all", "fractional", "ignore_secondary")

    Once the counts are returned the function takes care of extracting:
        - summary_counts_df:  [_total, _alignment_not_unique, _no_feature, _ambiguous]
        - count_df: actual raw Gene Counts Table (gxt)

    """

    ## get sample_file
    sample_file = os.path.basename(in_bam)
    sample_id = sample_file.split(".")[0]

    # -------------------------------------------------------------
    # ------    Feature Annotation Data: in_gdf vs in_gff     -----
    # -------------------------------------------------------------

    ## Import Gene Data Table (gdf)
    gdf = pd.read_csv(in_gdf, sep='\t')

    ## - Get 'repeat_ids' for using during counting
    repeat_ids = gdf[gdf['category'] == 'repeat']['gene_id'].tolist()

    ## - Get 'features' used for counting:
    ## For sequences, (strings, lists, tuples), use the fact that empty sequences are false.
    if not feature_type:
        feature_type = gdf['type'].unique()

    print("\nGetting `raw_gene_counts` using feature: {}".format(feature_type))
    print("Using `stranded` mode: {}".format(stranded))
    print("Using Count Parameters: \n     - `count_mode`: {}\n     - `ambiguous_assignment_mode`: {}\n     - `multimapped_mode`: {}\n".format(count_mode, ambiguous_assignment_mode, multimapped_mode))
    
    ## ------------------
    ## A. Gene Data Info: (currently, seems more practical to just use this)
    ## ------------------

    ## features: as `GenomicArrayOfSets` htseq object.
    features, genes_dict = get_features_from_gdf(gdf, stranded=stranded) ## use same feature as Parastou

    ## ------------------
    ## B. GFF Annotation: (not used anymore in our current pipe-line)
    ## ------------------

    ## features: as `GenomicArrayOfSets` htseq object.
    ## use same feature as Parastou - `gene`
    #features, genes_dict = get_features_from_gff(in_gff, stranded=stranded, feature_type=feature_type) 

    # -------------------------------------------------------------
    # -----------         Get feature counts          -------------
    # -------------------------------------------------------------

    # Parse `.bam` file and count 'reads' falling into each `feature`s genomic location
    #count_df, gene_counts_df = get_counts(in_bam, features, repeat_ids,
    count_df = get_counts(
        in_bam, features, repeat_ids, 
        count_mode=count_mode,
        ambiguous_assignment_mode=ambiguous_assignment_mode,
        multimapped_mode=multimapped_mode,
        max_nh=max_nh
        )
                                            
    ## -----------------------
    ## A. Summary of Counting:
    ## -----------------------

    ## Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
    summary_counts_df = count_df[count_df['gene_id'].str.startswith('_')]
    ## remove summary entries from `count_df`
    count_df = count_df[~count_df['gene_id'].str.startswith('_')]

    ## -------------------------
    ## B. raw Gene Counts Table: (gxt)
    ## -------------------------
    
    ## raw Gene Counts Table: (gxt)
    try:
        count_df = pd.merge(count_df, gdf[['gene_id', 'gene_name', 'transcript_length', 'gene_length', 'type', 'category']], on='gene_id', how='outer')
        #gene_counts_df = pd.merge(gene_counts_df, gdf[['gene_id', 'gene_name', 'transcript_length', 'gene_length', 'type', 'category']], on='gene_id', how='outer')
    except:
        gdf_cols = ['gene_id', 'gene_name', 'type', 'category']
        gdf_cols.extend([cc for cc in gdf.columns if 'length' in cc])
        count_df = pd.merge(count_df, gdf[gdf_cols], on='gene_id', how='outer')
        #gene_counts_df = pd.merge(gene_counts_df, gdf[gdf_cols], on='gene_id', how='outer')
        
    return count_df, summary_counts_df
    #return gene_counts_df, summary_counts_df


if __name__ == "__main__":

    ## -----------------
    ## Script Parameters
    ## -----------------

    params = vars(get_args())

    # initialize args/params passed to script ------------------------------
    
    param_in_bam = params['in_bam'] # path to alignment file (.bam or .sam)
    
    ## (Optional): INPUT subtraction - from BAM
    param_in_input = params['in_chip_input'] # path to `ChIP INPUT` alignment file (.bam or .sam)

    ## (Optional): INPUT subtraction - from counts
    param_in_count = params['in_count'] # path to `ChIP` count file 
    param_in_input_count = params['in_chip_input_count'] # path to `ChIP INPUT` count file
    
    ## (Optional): INPUT subtraction factor
    param_norm_factor = params['norm_factor']
    if not isinstance(param_norm_factor, type(None)):
        # at least one 'ChIP' file
        #bam_mode = ( not isinstance(param_in_bam, type(None)) ) and ( not isinstance(param_in_input, type(None)) )
        is_chip = ( not isinstance(param_in_bam, type(None)) ) or ( not isinstance(param_in_count, type(None)) )
        # at least one 'INPUT' file
        #count_mode = ( not isinstance(param_in_count, type(None)) ) and (not isinstance(param_in_input_count, type(None)) )
        is_input = ( not isinstance(param_in_input, type(None)) ) or (not isinstance(param_in_input_count, type(None)) )
        # one of them should be used
        #assert bam_mode or count_mode
        # both should exist
        assert is_chip and is_input
        print('ChIP INPUT subtraction factor:', param_norm_factor)
        
    else:
        # in_bam is required
        assert not isinstance(param_in_bam, type(None))
        
    param_in_gtf = params['in_gtf'] # reference Annotation File (.gdf)
    param_in_gdf = params['in_gdf'] # gene information table (gdf), containing features and corresponding attributes: gene_lengths, etc.
    
    param_prefix = params['prefix'] # name prefix (eg. `chip_' or 'rna_') for count tables. (default:'')
    param_verbose = params['verbose']
    param_out_dir = params['out_dir'] # save results in (default: current working directory).'
    
    ## --------------------------------
    ## Sequencing technology Parameters
    ## --------------------------------
    
    # Feature/s to be parsed from GFF and used for counting - currently, only used for naming the 'output_file's
    param_feature_type = list(set(params['feature_type']))
        
    ## Sequencing technology, ChIP, pA-RNA, RIP, total-RNA will determine:
    ##  - How to we deal with the `strandedness` of the data.
    ##  - Which features to use for counting: `genes` or `exons`
    ## `default` can be overwritten by explicitly passing a value.
    param_seq_type = params['seq_type']
    
    seq_param = seq_param_config(seq_type=param_seq_type, feature_type=param_feature_type)
    stranded = seq_param['stranded']
    feature_type = seq_param['feature_type']
    norm_length = seq_param['length']
    
    ## ----------------
    ## Count Parameters
    ## ----------------
    
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

    # -------------------------------------------------------------
    # -------------          Gene Counts             --------------
    # -------------------------------------------------------------
    
    # --------------------------
    # Workflow:
    # - Load `sample` BAM file.
    # - Get raw gene counts for features specified in Gene information table (gdf).
    # - Get TPM-normed expression for features specified in Gene information table (gdf).
    # - (Optional) load corresponding `INPUT` BAM file and subtract INPUT from ChIP counts.
    # --------------------------
    
    ## get sample_file
    if not isinstance(param_in_bam, type(None)):
        sample_file = os.path.basename(param_in_bam)
        sample_id = sample_file.split(".")[0]
    else:
        sample_file = os.path.basename(param_in_count)
        sample_id = '_'.join(sample_file.split("_")[0:3])
    
    # -------------------------------------------------------------
    # -----     Compute Gene Counts: for features in gdf      -----
    # -------------------------------------------------------------

    ## Call main function - raw_gene_counts_csv()
    ## Once the counts are returned the function takes care of extracting:
    ##    - summary_counts_df:  [_total, _alignment_not_unique, _no_feature, _ambiguous]
    ##    - count_df: actual raw Gene Counts Table (gxt)
    if isinstance(param_in_count, type(None)):
        
        count_df, summary_counts_df = raw_gene_counts_csv(
            param_in_bam,
            param_in_gdf,
            param_in_gtf,
            feature_type=feature_type,
            stranded=stranded,
            count_mode=param_count_mode,
            ambiguous_assignment_mode=param_ambiguous_assignment_mode,
            multimapped_mode=param_multimapped_mode,
            max_nh=param_max_nh
            )
            
    else:
        count_df = pd.read_csv(param_in_count, sep='\t')
        summary_counts_df = None
    
    # -------------------------------------------------------------
    # ------   (Optional) Subtract INPUT from ChIP counts   -------
    # -------------------------------------------------------------
    
    input_prefix = ''
    if not isinstance(param_norm_factor, type(None)):
        
        print('\n{}'.format('-'*45))
        print('Subtracting INPUT from ChIP `raw_gene_counts`')
        print('-'*45)
        
        # use prefix to rename files
        input_prefix = "subtracted_INPUT_"
        
        # Compute raw counts of corresponding: 'INPUT ChIP' sample
        if isinstance(param_in_input_count, type(None)):
            
            input_count_df, _ = raw_gene_counts_csv(
                param_in_input,
                param_in_gdf,
                param_in_gtf, 
                feature_type=feature_type,
                stranded=stranded,
                count_mode=param_count_mode, 
                ambiguous_assignment_mode=param_ambiguous_assignment_mode,
                multimapped_mode=param_multimapped_mode,
                max_nh=param_max_nh
                )
             
        else:
            input_count_df = pd.read_csv(param_in_input_count, sep='\t')
                                                        
        ## Get INPUT norm factor:                                                
        norm_factor = float(param_norm_factor)
        # Store copy of used INPUT norm factor:
        norm_factor_file = os.path.join(out_dir, sample_id, 'input_norm_factor.csv')
        pd.DataFrame({'ip_norm_factor':[norm_factor]}).to_csv(norm_factor_file, sep='\t', index=None)

        # Subtract 'INPUT' from 'ChIP': raw gene counts --------------
                
        # order both DataFrames with same column 'gene_id'         
        count_df = count_df.sort_values('gene_id').reset_index(drop=True)
        input_count_df = input_count_df.sort_values('gene_id').reset_index(drop=True)
        assert count_df.shape[0] == input_count_df.shape[0]
        
        # compute subtraction
        #count_df[['count', 'read_count']] += -norm_factor * input_count_df[['count', 'read_count']].fillna(0)
        count_df[['count']] += -norm_factor * input_count_df[['count']].fillna(0)
        
        # Attention! After INPUT subtraction there are counts with negative values!
        # => Need to remove them not to cause problems with TPM-normalization
        # TODO: how to check how many negative genes! (we are turning them to zero here)
        count_df['count'] = count_df['count'].clip(lower=0)
        

    ## Furthermore, convert these raw counts to TPM-normed expression using get_tpm_gene_counts().
    ## - tpm_df:  TPM-normed Gene Expression Table (tpm_gxt)
    tpm_df = get_tpm_gene_counts(count_df, count_col = 'count', length_col = norm_length)
    #tpm_df = get_tpm_gene_counts(tpm_df, count_col = 'read_count', length_col = norm_length)

    # -------------------------------------------------------------
    # ------------        Store Count Matrices        -------------
    # -------------------------------------------------------------
    
    prefix = param_prefix
    
    ## Used for naming output_files - select first feature
    if isinstance(feature_type, type(list())):
        feature_type = feature_type[0]

    ## Naming of Files:
    summary_count_file = os.path.join(out_dir, sample_id, 'summary_{}{}_{}pombe_{}_count_matrix.csv'.format(prefix, sample_id, input_prefix, feature_type))
    count_file = os.path.join(out_dir, sample_id, '{}{}_{}pombe_{}_count_matrix.csv'.format(prefix, sample_id, input_prefix, feature_type))
    tpm_file = os.path.join(out_dir, sample_id, '{}{}_{}pombe_tpm_matrix.csv'.format(prefix, sample_id, input_prefix))

    ## A. Store - Summary of Count Data as .csv files
    if not isinstance(summary_counts_df, type(None)):
        summary_counts_df.to_csv(summary_count_file, sep='\t', index=None)

    ## B. Store - raw Gene Counts Table (gxt) as .csv files
    count_df.to_csv(count_file, sep='\t', index=None)
    #gene_counts_df.to_csv(gene_count_file, sep='\t', index=None)

    ## C. Store - TPM-normed Gene Expression Table (tpm_gxt) as .csv files
    tpm_df.to_csv(tpm_file, sep='\t', index=None)


    #counts_df.shape
    # - Show counts for genes of interest: **Heterochromatic Genes**
    #counts_df[counts_df['gene-id'].isin(htc_genes)]
    
