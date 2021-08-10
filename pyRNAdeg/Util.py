#!/usr/bin/env python

import logging
from math import log
from math import isclose
import numpy as np
import os
import pandas as pd

import pdb

# pd.options.mode.chain_assignment = None

logger = logging.getLogger(__name__)


#header = ['gene-id', 'gene-name', 'length', 'type', 'category', 'bio_type']
header = ['gene_id', 'gene_name', 'length', 'type', 'category', 'bio_type']

## used to filter-out all columns but samples - works both for Parastou and new version
long_header = [
  'gene_id', 'gene-id', 'gene_name', 'gene-name', 
  'seqid', 'chr', 'type', 'start', 'end', 'strand',
  'cds_length', 'utr_length', 'length', 'category', 'bio_type'
  ]

## used to filter-out all columns but samples - works both for Parastou and new version
# => contains `gene_length` and `intron_length`
long_header_v2 = [
  'gene_id', 'gene-id', 'gene_name', 'gene-name', 
  'seqid', 'chr', 'type', 'start', 'end', 'strand',
  'gene_length', 'intron_length', 'cds_length', 'utr_length', 'length', 'category', 'bio_type'
  ]
  
## ----------------------------
## reference `GTF/GFF` **file**
## ----------------------------

def parse_attribute_col(element):

    if not pd.isnull(element):
        ## key-value pair - k:v
        return {ii.split('=')[0]:ii.split('=')[1] for ii in element.split(';')}

    else:
        #import pdb
        #pdb.set_trace()

        return {}


def parastou_gff_to_gdf(in_gff):
    columns = ['chr', '1', 'type', 'start', 'end', '2', '3', '4', 'info']

    df = pd.read_csv(in_gff, sep='\t', comment='#', names=columns)
    df = df[['chr', 'type', 'start', 'end', 'info']]
    df = df[df['type'].str.contains('gene')]

    df['gene-id'] = df['info'].apply(lambda x: x.split('=gene:')[1].split(';')[0].strip(' ";,'))
    df['gene-name'] = df['info'].apply(lambda x: x.split(';Name=')[1].split(';')[0].strip(' ";,'))
    df['bio_type'] = df['info'].apply(lambda x: x.split(';biotype=')[1].split(';')[0].strip(' ";,'))

    df = df[['gene-id', 'gene-name', 'chr', 'type', 'start', 'end', 'bio_type']]

    return df


def gff_to_gdf(in_gff):

    # columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    columns = ['chr', '1', 'type', 'start', 'end', '2', 'strand', '4', 'info']

    df = pd.read_csv(in_gff, sep='\t', comment='#', names=columns)
    ## Select columns of interest
    df = df[['chr', 'type', 'start', 'end', 'strand', 'info']]

    ## Filter for type of feature: 'genes':
    df = df[df['type'].str.contains('gene')]

    ## Parse `attributes` colummn
    attributes_df = df['info'].apply(lambda x: parse_attribute_col(x))
    # attributes_df = attributes_df.apply(pd.Series) ## slower
    attributes_df = pd.DataFrame(attributes_df.values.tolist(), index=df.index)
    df = pd.concat([df, attributes_df], axis=1)

    ## because Parastou used the `ID` instead of `gene_id` attribute - otherwise the manual entries will contains NaNS
    df['gene-id'] = df['ID'].apply(lambda x: x.split(':')[1])

    ## rename some columns to fit parastous definitions...
    df = df.rename(columns={"Name": "gene-name", "biotype": "bio_type"})

    ## Select columns of interest
    df = df[['gene-id', 'gene-name', 'chr', 'type', 'start', 'end', 'strand', 'bio_type']]

    return df


def add_gene_length(in_gff, gdf):

    columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    df = pd.read_csv(in_gff, sep='\t', comment='#', names=columns)
    ## Select columns of interest
    df = df[['seqid', 'type', 'start', 'end', 'strand', 'attributes']]

    ## --------------------------
    ## Filter for type of feature: 'transcript'
    ## --------------------------

    transcript_df = df[df['type'].str.contains('transcript')]

    ## Parse `attributes` colummn for 'transcript_df'
    attributes_df = transcript_df['attributes'].apply(lambda x: parse_attribute_col(x))
    # attributes_df = attributes_df.apply(pd.Series) ## slower
    attributes_df = pd.DataFrame(attributes_df.values.tolist(), index=transcript_df.index)
    transcript_df = pd.concat([transcript_df[['seqid', 'type', 'start', 'end', 'strand']], attributes_df], axis=1)

    ## used the `Parent` to get the `gene-id` attribute
    transcript_df['gene-id'] = transcript_df['Parent'].apply(lambda x: x.split(':')[1])

    ## rename some columns to fit parastous definitions...
    transcript_df = transcript_df.rename(columns={"transcript_id": "transcript-id"})

    ## --------------------------
    ## Filter for type of feature: 'exons'
    ## --------------------------

    exon_df = df[df['type'].str.contains('exon')]

    ## Parse `attributes` colummn for 'exon_df'
    attributes_df = exon_df['attributes'].apply(lambda x: parse_attribute_col(x))
    # attributes_df = attributes_df.apply(pd.Series) ## slower
    attributes_df = pd.DataFrame(attributes_df.values.tolist(), index=exon_df.index)
    exon_df = pd.concat([exon_df[['seqid', 'type', 'start', 'end', 'strand']], attributes_df], axis=1)

    ## used the `Parent` to get the `transcript-id` attribute
    exon_df['transcript-id'] = exon_df['Parent'].apply(lambda x: x.split(':')[1])

    ## compute exon lengths
    exon_df['length'] = exon_df['end'] - exon_df['start']

    ## ------------
    ## Summarize by: 'transcript'
    ## ------------

    transcript_lengths_df = exon_df.groupby('transcript-id').sum()
    transcript_lengths_df = pd.merge(transcript_lengths_df, transcript_df[['gene-id', 'transcript-id']],
                                     on='transcript-id')

    return transcript_lengths_df


def get_category(gene_id, gene_name, repeat_names):

    if not pd.isnull(gene_name):

        if 'rpl' in gene_name or 'rps' in gene_name:
            return 'ribosomal'

    if gene_id in repeat_names:

        return 'repeat'

    else:
        return 'gene'



## --------------------
## Transform DataFrames
## --------------------

def to_log2_tpm(tpm_df, cols=None, gene_id_col='gene-id', shift=1):
    
    # select 'header' columns in df
    #hh = [xx for xx in tpm_df.columns if xx in long_header]
    hh = [xx for xx in tpm_df.columns if xx in long_header_v2]
    new_df = tpm_df[hh].copy()
    
    # only apply log transformation to a subset of samples
    if cols:
        sample_columns = cols
    else:
        # use all samples
        sample_columns = [ss for ss in tpm_df.columns if ss not in hh]
        
    # # keep tracks of columns containing zeros
    # zero_cols = 0
    # # loop over samples to apply log transformation
    # for sample_col in sample_columns:
    #                 
    #     # apply log transformation to column
    #     sample_df, zero_col = to_log2_col(tpm_df[[gene_id_col, sample_col]], sample_col, shift=shift)
    # 
    #     new_df = pd.merge(new_df, sample_df, on=gene_id_col, how='left')
    #     
    #     new_df = new_df.drop_duplicates()
    #     zero_cols += zero_col
    
    # TODO: how to handle NAs during log transformation?
    assert not tpm_df[sample_columns].isnull().values.any()
    
    # loop over samples to apply log transformation
    for sample_col in sample_columns:

        # apply log transformation to column
        new_df[sample_col] = tpm_df[sample_col].apply(
          lambda x: log((x + shift), 2) if not isclose(x + shift, 0., abs_tol=0.0) else -np.inf)
    
    # keep tracks of columns that contained zeros and might create issues (for shift=0)
    # => log[0] = -inf
    zero_cols = np.isinf(new_df[sample_columns]).any(0).sum()
    
    if zero_cols:
        #print('%d columns contained zero values. Their log-transformed results are NaNs' % zero_cols)
        print('%d columns contained zero values. Their log-transformed results are -inf' % zero_cols)

    return new_df


def to_log2_col(df, col, shift=1):

    zero_shift = 0
    
    # without shift - remove rows with zeros
    # => log[0] = -inf
    if shift == 0:
        
        # filter rows containing zero
        new_df = df[df[col] != 0]
        
        # keep track that at least one row was removed
        if len(new_df) < len(df):
            zero_shift = 1
            
    # with shift - make copy of df
    else:
        new_df = df.copy()
    
    # apply transform
    new_df[col] = new_df[col].apply(lambda x: log((x + shift), 2))

    return new_df, zero_shift


def zero_adjust(x, shift=0.1):
  
    if x[0] == 0.0 or x[1] == 0.0:
        x[0] += shift
        x[1] += shift
        
    return x[0], x[1]


def zero_adjust_df(df, pairs):

    new_header = [x for x in header if x in list(df.columns)]
    za_df = df[new_header]

    for pair in pairs:
        t_s = df[list(pair)].apply(zero_adjust, axis=1)
        t_d = t_s.apply(pd.Series)
        t_d.columns = list(pair)
        za_df = pd.concat([za_df, t_d], axis=1)

    return za_df


def ratio_df(df, pairs, new_cols=None):

    #new_header = [x for x in header if x in list(df.columns)]
    new_header = [x for x in header if x in df.columns]

    if new_cols:
        assert len(pairs) == len(new_cols)

    else:
        new_cols = [x + ',' + y for x, y in pairs]

    for index, pair in enumerate(pairs):

        x, y = pair
        name = new_cols[index]
        df[name] = df[x] / df[y]

    df = df.replace([np.inf, -np.inf], np.nan)
    gxp_df = df[new_header + new_cols]

    return gxp_df


def ratio_table_old(df, pairs, new_cols=None, out_dir='.', file_name=None):

    if not file_name:
        file_name = 'te_table.csv'

    out_path = os.path.join(out_dir, file_name)
    gxp_df = ratio_df(df, pairs, new_cols=new_cols)
    gxp_df.to_csv(out_path, sep='\t', index=None)

    return gxp_df

  
def ratio_table(df, pairs_dict, shift=0, out_dir='.', file_name=None):
    
    #new_header = [x for x in header if x in df.columns]
    new_header = [x for x in header if x in list(df.columns)]
    new_df = df[new_header].copy()

    for pair_id, pair in pairs_dict.items():
        
        if len(pair) == 2:
          
            # get 'seq' pair for ratio x/y
            x, y = pair
            
            try:
                new_df.loc[:, pair_id] = (shift + df[x]) / (shift + df[y])
            
            # mutant samples might have been excluded (not present)
            except:
                print('Mutant: {}'.format(pair_id))
                print(' - Sample {} in Data Frame: {}'.format(x, x in df.columns) )
                print(' - Sample {} in Data Frame: {}'.format(y, y in df.columns) )

        elif len(pair) == 1:
            new_df[pair_id] = np.nan
            
        else:
            import pdb; pdb.set_trace()
    
    assert shift >= 0 # won't work for negative shifts
    if not shift:
      # take care of (inf, -inf)        
      new_df = new_df.replace([np.inf, -np.inf], np.nan)
      # TODO: should do the same for zeros?
      new_df = new_df.replace(0, np.nan)
    
    # store copy
    if not file_name:
        file_name = 'te_table.csv'
    out_path = os.path.join(out_dir, file_name)
    new_df.to_csv(out_path, sep='\t', index=None)

    return new_df


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to RNAdeg utility module.')
    logger.info('')
    logger.info('-' * 40)





