"""
Author: P. Monteagudo
Affiliation: LMU, GC
Aim: Snakefile containing rules for ChIP-seq Data Analysis Pipe-line.

Date: November 25  2019
Run:    snakemake -s preprocess_ChIP.smk -n
        snakemake -s preprocess_ChIP.smk -p
        snakemake -s preprocess_ChIP.smk -p -j 20 -R gene_expression_table
        snakemake -s preprocess_ChIP.smk -p -j 20 -R chip_gene_expression_tables ## make sure we use input_subtracted
        
        # Touch output files (mark them up to date without really changing them) instead of running their commands.
        snakemake -s preprocess_ChIP.smk -t

        # Number of Jobs
        snakemake -s preprocess_ChIP.smk -p -j 20
        # DAG Jobs
        snakemake -s preprocess_ChIP.smk --dag | dot -Tpdf > jobs_dag.pdf
        # DAG Rules
        snakemake -s preprocess_ChIP.smk --rulegraph | dot -Tsvg > rules_dag.svg
        ## Go on with independent jobs if a job fails.
        snakemake -s preprocess_ChIP.smk -p -j 20 -k
        ## Force the re-execution or creation of the given rules or files.
        snakemake -s preprocess_ChIP.smk -p -j 20 -R construct_TI_models visualize_TI_models

        ## Cluster
        snakemake -s preprocess_ChIP.smk -p --cluster "sbatch" --jobs 4
        snakemake -s preprocess_ChIP.smk -p --cluster "sbatch --mem=2G --partition=hi_mem --nodelist=node3 --output=slurm_job_files/slurm-%j.out" --jobs 100
        snakemake -s preprocess_ChIP.smk -p --cluster "sbatch --mem=2G --partition=full_high --nodelist=master --output=slurm_job_files/slurm-%j.out" --jobs 1
        snakemake -s preprocess_ChIP.smk -p --cluster "sbatch --mem=2G --partition=full_high --output=slurm_job_files/slurm-%j.out" --jobs 1000 -k
        snakemake -s preprocess_ChIP.smk -p --cluster "sbatch --mem=2G --time=0-03:00:00 --partition=full_high --output=slurm_job_files/slurm-%j.out" --jobs 500 -k -R unzip_fastq_files

Latest modification:
  - todo
"""

import os

################
# Config Files #
################

#configfile: "../config.yaml"
configfile: "../config_algbio.yaml"

#####################
# Working Directory #
#####################

# All paths in the snakefile are interpreted relative to the directory snakemake is executed in.
# This behaviour can be overridden by specifying a workdir:
workdir: config["work_dir"]

#-------------------------------------------------------------
#----------------  Important Parameters ----------------------
#-------------------------------------------------------------

#MAX_THREADS = 32
MAX_THREADS = 1
N_BENCHMARKS = 1

######################
# Set up Singularity #
######################
#singularity:
#	"docker://continuumio/miniconda3:4.5.12"d

#################
# Set up Report #
##################
#report: "report/workflow.rst"

############################
# Include other Snakefiles #
############################
#include:
#	"rules/misc.smk"

import pandas as pd

#-------------------------------------------------------------
#------------        Auxiliary functions         -------------
#-------------------------------------------------------------

def concatenate_gene_count_tables(list_filenames, list_df = None, assert_size = None, script_mode='htseq', count_col='count', join='outer'):
    
    # pass objects as DataFrames instead of as files to be loaded
    if isinstance(list_df, type(None)):
        list_df = []
    else:
        list_df = list_df.copy()
    
    # iterate over files to be loaded, rename columns to common framework:
    for filename in list_filenames:

        sample_id = os.path.basename(os.path.dirname(filename))

        if script_mode == 'htseq':
            df = pd.read_csv(filename, sep='\t', index_col='gene_id', usecols=['gene_id', count_col])
            df = df.rename(columns={count_col: sample_id})

        else:
            df = pd.read_csv(filename, sep='\t', usecols=['gene-id', sample_id])
            df = df.rename(columns={'gene-id':'gene_id'})
            df = df.set_index('gene_id')

        #df = df.sort_values(['seqid', 'start', 'end', 'strand'])
        df = df.sort_index()

        if not isinstance(assert_size, type(None)):
            assert len(df.index) == assert_size

        list_df.append(df)
    
    # pool - concatenate all DataFrames
    return pd.concat(list_df, axis=1, sort=True, join=join)

def concatenate_dataframes(list_filenames, list_df = None, assert_size = None, script_mode='htseq', count_col='count', join='outer'):
    
    # pass objects as DataFrames instead of as files to be loaded
    if isinstance(list_df, type(None)):
        list_df = []
    else:
        list_df = list_df.copy()
    
    # iterate over files to be loaded, rename columns to common framework:
    for filename in list_filenames:

        sample_id = os.path.basename(os.path.dirname(filename))

        # if script_mode == 'htseq':
        #     df = pd.read_csv(filename, sep='\t', index_col='gene_id', usecols=['gene_id', count_col])
        #     df = df.rename(columns={count_col: sample_id})
        # 
        # else:
        #     df = pd.read_csv(filename, sep='\t', usecols=['gene-id', sample_id])
        #     df = df.rename(columns={'gene-id':'gene_id'})
        #     df = df.set_index('gene_id')
            
        df = pd.read_csv(filename, sep='\t', index_col='gene_id', usecols=['gene_id', count_col])
        df = df.rename(columns={count_col: sample_id})
            
        #df = df.sort_values(['seqid', 'start', 'end', 'strand'])
        df = df.sort_index()

        if not isinstance(assert_size, type(None)):
            assert len(df.index) == assert_size

        list_df.append(df)
    
    # pool - concatenate DataFrames
    return pd.concat(list_df, axis=1, sort=True)
    
    
#-------------------------------------------------------------
#---------    Pipeline Design - Count Parameters     ---------
#-------------------------------------------------------------

## Potentially allow multimapped-reads - with less than max_nh (default:16) hits
## - We only allow this if the read maps to a repeat feature.
max_nh = 16
#max_nh = 15 # test if inconsistency arises from this!

## - How to deal with overlap of 'READ' and 'FEATURE'?
#count_mode = 'intersection-nonempty'
count_mode = 'intersection-strict'
#count_mode = 'union'
# ("union", "intersection-strict", "intersection-nonempty")

## - How to deal with multiple-overlapping 'FEATURES'? (according to count_mode)
ambiguous_assignment_mode = 'all'
#ambiguous_assignment_mode = 'all' 
# ("none", "all", "fractional")

## - How to deal with multi-mapped 'READS'?
#multimapped_mode = 'fractional' # with NH-norm, does it matter due to ratios? => Not for ratios, but for some plots!
multimapped_mode = 'all' # with-out NH-norm
# ("none", "all", "fractional", "ignore_secondary")

## - Prefix - with Count Parameters:
prefix = multimapped_mode + '_' + count_mode + '_'


#-------------------------------------------------------------
#---------  Pipeline Design - Data Sets Information  ---------
#-------------------------------------------------------------
   
# use `simulated_data` mode: (or `H3K9me2`)
# => Attention! 
#    For simulated-data ChIP, can't apply regular pipeline because of INPUT subtraction.
#    To get corresponding 'chip_pombe_tpm_matrix.csv' remember to modify input in:
#        - group_gene_count_tables rule
#        - chip_gene_expression_tables rule
simulated_data = False
#simulated_data = True

## ---------------
## - ChIP Datasets: | S2-ChIP | S5-ChIP | INPUT |
## ---------------

# import DataFrame containing `sample` annotation and filter for samples that belong to the ChIP Pipeline
all_samples_df = pd.read_csv(config["sample_annotation_file"], sep = '\t')

if not simulated_data:
  all_samples_df = all_samples_df[all_samples_df['pipeline_type'] == 'ChIP']
else:
  # select 'simulated-data' with mutant_id 'chip-fake-reads'
  #all_samples_df = all_samples_df[ (all_samples_df['pipeline_type'] == 'simulated-data', 'H3K9me2') & (all_samples_df['mutant_id'].str.contains("chip")) ]
  # select 'H3K9me2' and 'simulated-data' with mutant_id 'chip-fake-reads' (NOT 'rna-fake-reads')
  all_samples_df = all_samples_df[ (all_samples_df['seq_category'].isin(['simulated-data', 'H3K9me2'])) & (~all_samples_df['mutant_id'].str.contains("rna")) ]
  # over-write 'pipeline_type'
  all_samples_df['pipeline_type'] = 'ChIP'

# Ignore S5-samples for now!
#all_samples_df = all_samples_df[all_samples_df['seq_type'] != 'S5-ChIP']

# rename 'sample_id' column to 'id'
all_samples_df = all_samples_df.rename(columns = {'sample_id' : 'id'})
#all_samples_df = all_samples_df.iloc[0:2, :]

## -----------------
## - Select Data Set
## -----------------
#import pdb; pdb.set_trace()

# by `dataset_id`
select_dataset = ""
#select_dataset = "1022_S2-ChIP_1"
#select_dataset = "WT_S2-ChIP-INPUT_1"
select_dataset = "WT_S2-ChIP_3"
#select_dataset = "WT_S2-ChIP-INPUT_1"
#select_dataset = "WT_H3K9me2_1"
#select_dataset = "chip-fake-reads_simulated-data_2"
#all_samples_df = all_samples_df[all_samples_df.id == select_dataset]

# by `dataset_id`s
select_datasets = ["1022_S2-ChIP_1", "1022_S2-ChIP-INPUT_1"]
#all_samples_df = all_samples_df[all_samples_df.id.isin(select_datasets)]

# by `mutant_id`
select_mutant = "WT_"
select_mutant = "WT_S2"
# select_mutant = "WT_H3K9me2"
# select_mutant = "638_H3K9me2"
# select_mutant = "H3K9me2"
#select_mutant = "chip-fake-reads"
#all_samples_df = all_samples_df[all_samples_df['id'].str.contains(select_mutant)]

## ------------------
## - Ignore Data Sets
## ------------------

# Ignore Data Sets
ignore_datasets = []
all_samples_df = all_samples_df[~all_samples_df.id.isin(ignore_datasets)]

## --------------------
## - Truncate Data Sets
## --------------------

## Error - EXITING because of FATAL ERROR in reads input: short read sequence line: 1
## Need to trim reads:
# truncate_datasets = ['63_S2ChIPp']
# #truncate_datasets.extend([])
truncate_datasets = all_samples_df[all_samples_df['trimmed']].id.values
#all_samples_df = all_samples_df[all_samples_df.id.isin(truncate_datasets)]


## --------------------------------------------------------------------------------
## Big difference in the ChIP-seq pipeline is that we need to use the INPUT samples
## to subtract noise from ChIP samples coverage/counts.
## --------------------------------------------------------------------------------

## Note: `input_types` refers to INPUT samples, but note that we also include other 'ChIP'-like samples
## that do NOT need 'INPUT substraction' during pre-processing. (e.g. 'H3K9me2', 'simulated-data' ).
#input_types = ['S2-ChIP-OIN', 'S2-ChIP-INPUT', 'simulated-data', 'H3K9me2']
input_types = ['INPUT', 'H3K9me2', 'simulated-data']
#input_types = ['INPUT', 'simulated-data']

## -----------------
## A. ChIP Datasets:
## -----------------

## Distinguish between (actual) 'ChIP' and 'INPUT' samples
## => samples that NEED INPUT subtraction
#datasets_df = all_samples_df[~all_samples_df['seq_type'].isin(input_types)]
datasets_df = all_samples_df[~all_samples_df['seq_category'].isin(input_types)]
#datasets_df = datasets_df.iloc[0:2, :]
#datasets_df = datasets_df.iloc[0:1, :]


## ------------------
## B. INPUT Datasets: (and `simulated-data`)
## ------------------

#data_input_batch = 'sequencing_new/ChIP'

## Distinguish between (actual) 'ChIP' and 'INPUT'/simulated-data samples
## => samples that DO NOT NEED INPUT subtraction
#input_datasets_df = all_samples_df[all_samples_df['seq_type'].isin(input_types)]
input_datasets_df = all_samples_df[all_samples_df['seq_category'].isin(input_types)]
#input_datasets_df = input_datasets_df.iloc[0:2, :]
 
# input_dataset = {"524_1_S2ChIP": "524_1_INPUT",
#                 "63_S2ChIPp": "63_OIN"}
#                 
# input_norm_factor = {"524_1_S2ChIP": 0.597,
#                     #"63_S2ChIPp": 0.785} ## THOR
#                     "63_S2ChIPp": 0.137} ## Centromeric Norm

#import pdb; pdb.set_trace()

#######################################
# Convienient rules to define targets #
#######################################

## Mark a rule as local, so that it is not submitted to the cluster and instead executed on the host node
localrules: all

rule all:
    input:
        # # config["chip_star_idx"],
        # ## ---------------------
        # ## 1. Unzip fastq files:
        # ## ---------------------
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/fastq/{dataset_id}.fastq"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        # ## ---------------------
        # ## 2. Align fastq files: using `STAR`
        # ## ---------------------
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.out.bam"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        # ## ----------------------------------
        # ## 3. Sort and Index read alignments: using `samtools`
        # ## ----------------------------------
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        ## -------------------------------------------------------
        ## 4. (Optional) Remove rRNA/Tag bam files with gene-id's
        ## ------------------------------------------------------
        ## => In ChIP data rRNA should not be an issue, but running into issues with normalization
        # by sample: tagged/filtered bam files
        expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        ## summarize samples by 'seq_category': (naive) ribosomal RNA and other RNA features Counts as .csv files:
        expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/ribosomal_rna_pombe_gene_count_matrix.csv"), seq_category=all_samples_df['seq_category'].unique()),
        expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/chip_pombe_gene_count_matrix.csv"), seq_category=all_samples_df['seq_category'].unique()),
        ## --------------------------------------
        ## 5. Estimate subtraction coefficients: ChIP - l * INPUT
        ## --------------------------------------
        ## ChIP & INPUT: BAM files
        expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/INPUT_factors/{dataset_id}/INPUT_factors.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        ## -------------------
        ## 5. Compute Coverage
        ## -------------------
        # ## ChIP & INPUT: Coverage (normally NOT necessary)
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/int_{dataset_id}_coverage.wig"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/frac_{dataset_id}_coverage.wig"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        # ## ChIP: Coverage subtracted INPUT
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/int_subtracted_INPUT_{dataset_id}_coverage.wig"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/frac_subtracted_INPUT_{dataset_id}_coverage.wig"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        ## -------------------------------------------------
        ## 6. Compute (raw) and TPM-normed Gene Count tables
        ## -------------------------------------------------
        ## ---------
        ## 6a. Genes
        ## ---------
        # # by sample: raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        # # by sample (INPUT subtracted): raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_subtracted_INPUT_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_subtracted_INPUT_pombe_tpm_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # summarize samples by 'seq_category'
        expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_gene_count_matrix.csv"), seq_category=all_samples_df['seq_category'].unique()),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_tpm_matrix.csv"), seq_category=all_samples_df['seq_category'].unique()),
        # summarize all `ChIP`samples
        os.path.join(config["data_dir"], "results/xp_data/ChIP/chip_pombe_gene_count_matrix.csv"),
        os.path.join(config["data_dir"], "results/xp_data/ChIP/chip_pombe_tpm_matrix.csv"),
        ## -----------
        ## 6b. Introns
        ## -----------
        # by sample: raw Intron Counts Table (gxt) and TPM-normed Intron Expression Table (tpm_gxt) as .csv files
        expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_intron_count_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        # summarize samples by 'seq_category'
        # #expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/pombe_intron_read_count_matrix.csv"), seq_category=all_samples_df['seq_category']),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/pombe_intron_count_matrix.csv"), seq_category=all_samples_df['seq_category'].unique()),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/pombe_tpm_matrix.csv"), seq_category=all_samples_df['seq_category'].unique()),

################
# Rules Proper #
################

'''
- 0. Get list of files in Directory:
'''

# not necesary for now

'''
- 1. Unzip fastq files:
'''

rule unzip_fastq_files:
    input:
        ## (zipped) raw data:
        raw = os.path.join(config["data_dir"], "seq_data/{seq_category}/raw_data/{dataset_id}.fastq.bz2"),

    output:
        ## .fastq file
        fastq = os.path.join(config["data_dir"], "seq_data/{seq_category}/fastq/{dataset_id}.fastq"),

    params:
        fastq_dir =  os.path.join(config["data_dir"], "seq_data/{seq_category}/fastq"),
        unzip_file = lambda wildcards, input: os.path.splitext(input.raw)[0],
        ## handling truncated reads
        bioawk = config["bioawk"],
        ztr_fastq = os.path.join(config["data_dir"], "seq_data/{seq_category}/fastq/{dataset_id}.ztr.fastq"),

    threads: 1

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/1_unzip_fastq_files/{dataset_id}.txt", N_BENCHMARKS)

    # log:
    #     os.path.join(config["data_dir"], "seq_data/{seq_category}/fastq/log/{dataset_id}.log")

    # shell:
    #     "{config[bzip2]} -kdv {input.raw}; mv {params.unzip_file} {output.fastq} 2> {log}"

    run:
        shell("{config[bzip2]} -kdv {input.raw}; mv {params.unzip_file} {output.fastq}")
        if wildcards['dataset_id'] in truncate_datasets:

            shell(
            """
            mv {output.fastq} {params.ztr_fastq}
            {params.bioawk} -cfastx 'length($seq) > 0 {{print "@"$name"\\n"$seq"\\n+\\n"$qual}}' {params.ztr_fastq} >> {output.fastq}
            """
            )

'''
- 2.  Align fastq files: using `STAR`
'''

# Note: Very fast, I guess because there is no GTF/GFF
rule star_index_chip:
    input:
        # reference FASTA file
        #fa = os.path.join(config["genome_dir"], 'spombe_v2/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa'),
        fa = config["fasta_file"],

    output:
        ## Contains no gtf - diffferent from RNA analysis!!
        #directory(os.path.join(config["genome_dir"], "spombe_v2/star_nogtf_idx")),
        directory(config["chip_star_idx"])

    threads: 20 # set the maximum number of available cores

    shell:
        'mkdir {output} && '
        '{config[star]} --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fa} '
        #'--genomeSAindexNbases 11 ' ## small genomes: spombe ~ 12.82 Mb, min(14, log2(GenomeLength/2) - 1) ~ 11
        ## WARNING: --genomeSAindexNbases 11 is too large for the genome size=13893632, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 10
        '--genomeSAindexNbases 10 ' ## small genomes: spombe ~ 12.82 Mb, min(14, log2(GenomeLength/2) - 1) ~ 11
        ## WARNING: --genomeSAindexNbases 14 is too large for the genome size=13893632, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 10
        #'2> {log}'

rule align_reads:
    input:
        ## STAR index
        #ref_dir = os.path.join(config["genome_dir"], "spombe_v2/star_nogtf_idx"),
        ref_dir = config["chip_star_idx"],
        ## .fastq file
        fastq = os.path.join(config["data_dir"], "seq_data/{seq_category}/fastq/{dataset_id}.fastq"),

    output:
        ## .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.out.bam"),
        #bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),

    params:
        bam_prefix = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}."),
        max_nh = max_nh,

    threads: 20

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/2_align_reads/{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/log/align_reads.log")

    shell:
        # EXITING because of FATAL ERROR: Genome version: 20201 is INCOMPATIBLE with running STAR version: 2.7.3a
        # SOLUTION: please re-generate genome from scratch with running version of STAR, or with version: 2.7.1a
        ## conda install -c bioconda star=2.7.1a-0
        '{config[star]} --runThreadN {threads} '
        '--genomeDir {input.ref_dir} '
        '--readFilesIn {input.fastq} '
        '--outFileNamePrefix {params.bam_prefix} '
        '--outSAMtype BAM Unsorted ' # sort manually with samtools, due to error
        #'--outSAMtype BAM SortedByCoordinate '
        ## max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped.
        '--outFilterMultimapNmax {params.max_nh} '
        #'--outFilterMultimapNmax 16 '
        #'--outFilterMultimapNmax 25 '
        #'--outFilterMultimapNmax 99 '
        ## this does not allow for splicing!
        '--alignIntronMax 1 '
        '--alignEndsType EndToEnd '
        '2> {log}'

'''
- 3.  Sort and Index read alignments: using `samtools`
'''

rule sort_aligned_reads:
    input:
        ## .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.out.bam"),

    output:
        ## (indexed) .bam.bai file
        bam_sort = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/3a_sort_aligned_reads/{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/log/sort_aligned_reads.log")

    shell:
        "{config[samtools]} sort -o {output.bam_sort} {input.bam} 2> {log}"

rule index_aligned_reads:
    input:
        ## .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),

    output:
        ## (indexed) .bam.bai file
        bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/3_index_aligned_reads/{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/log/index_aligned_reads.log")

    shell:
        "{config[samtools]} index {input.bam} 2> {log}"


'''
- 4. Remove rRNA/tag bam files with gene-id's
'''

rule filter_bam:
    input:
        ## .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        ## (indexed) .bam.bai file
        bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),

    output:
        ## (tagged) .bam file
        tagged_bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"),
        ## (tagged) indexed .bam.bai file
        tagged_bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam.bai"),
        ## Naive Ribosomal RNA (rRNA) and mRNA Gene Counts computed when parsing BAM file as .csv files
        rrna_counts_df = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/rrna_counts.csv"),
        gene_counts_df = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/gene_counts.csv"),

    params:
        ## This script removes ribosomal-RNA (rRNA) related reads from given .bam files
        rem_rrna_script = 'htseq/scripts/filter_bam.py', ## my script
        #rem_rrna_script = 'pyRNAdeg/RemrRNA.py', ## Parastou's script
        ## Main `read alignments` dir
        bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/bam'),
        ## Gene Information Table (GDF) file
        gdf = config["gdf_file"],
        tagged_bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/tagged_bam'),
        ## -----------
        ## Parameters:
        ## -----------
        ## - `pipeline_type`: distinguish between ChIP (Unstranded) and RNA (Stranded)
        pipeline_type = lambda wildcards: all_samples_df[all_samples_df['id'] == wildcards.dataset_id].pipeline_type.values[0],

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/4a_remove_rrna_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/log/remove_rrna.log")

    shell:
        "python {params.rem_rrna_script} "
        "--in_bam {input.bam} "
        "--in_gdf {params.gdf} "
        "--pipeline_type {params.pipeline_type} "
        "-o {params.tagged_bam_dir} 2> {log}"

# NOTE: Need to modify input files in `simulated_data` mode: (or `H3K9me2`)
rule group_gene_count_tables:
    input:
        ## Ribosomal RNA (rRNA) and mRNA Gene Counts computed when parsing BAM file as .csv files
        all_rrna_counts_df = expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/rrna_counts.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        all_gene_counts_df = expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/gene_counts.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        #all_rrna_counts_df = expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/rrna_counts.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']), # for simulated-data
        #all_gene_counts_df = expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/gene_counts.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']), # for simulated-data

    output:
        ## raw and TPM-normed Gene Counts Data (gxt)  as .csv files
        rrna_counts_df = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/ribosomal_rna_pombe_gene_count_matrix.csv"),
        gene_counts_df = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/chip_pombe_gene_count_matrix.csv"),

    params:
        ## Distinguish between: 'parastous' and 'htseq'
        script_mode = 'htseq',
        ## Gene Information Table GDF file
        gdf = config["gdf_file"],
        ## Column to select from Data Frame:
        #count_col = 'read_count',
        count_col = 'count',

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/4b_group_gene_count_tables.txt", N_BENCHMARKS)
        
    run:
        ## Init lists of individual sample DataFrames
        list_df = []
        
        ## ---------------------------
        ## Load Gene Information Table: gdf
        ## ---------------------------
        
        gdf = pd.read_csv(params['gdf'], sep='\t')
        
        ## rename columns to common framework: (parastous version does not contain so much information)
        if 'gene-id' in gdf.columns:
            gdf = gdf.rename(columns={'gene-id':'gene_id', 'gene-name':'gene_name'})
            gdf_header = ['chr', 'type', 'start', 'end', 'gene_id', 'gene_name', 'length', 'category', 'bio_type']
        
        ## my version:
        else:
            #gdf = gdf.rename(columns={'Name':'gene_name', 'gene_length':'length'})
            #gdf_header = ['gene_id', 'gene_name', 'seqid', 'type', 'start', 'end', 'strand',  'cds_length', 'utr_length', 'length', 'category', 'bio_type']
            gdf_header = ['seqid', 'type', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'cds_length', 'utr_length', 'intron_length', 'gene_length', 'category', 'bio_type']

        gdf = gdf[gdf_header]
        gdf = gdf.set_index('gene_id')
        gdf = gdf.sort_index()

        ## Add gdf to list of DataFrames
        list_df.append(gdf)

        ## -----------------------------------------
        ## A. ribosomal RNA (rRNA) Gene Counts Data:
        ## -----------------------------------------
        
        ## Concatenate DataFrames: returns a DataFrame with all features but most will be empty (only rRNA entries should contain counts)
        grouped_df = concatenate_gene_count_tables(input['all_rrna_counts_df'], list_df, script_mode=params['script_mode'], count_col= params['count_col'], join='outer')
        #grouped_df = concatenate_gene_count_tables(input['all_rrna_counts_df'], script_mode=params['script_mode'], join='outer')
        grouped_df.to_csv(output['rrna_counts_df'], sep='\t', index_label = ['gene_id'] )

        ## -------------------------
        ## B. mRNA Gene Counts Data: (everything but rRNA)
        ## -------------------------

        ## Concatenate DataFrames
        grouped_df = concatenate_gene_count_tables(input['all_gene_counts_df'], list_df, script_mode=params['script_mode'],count_col= params['count_col'], join='outer')
        #grouped_df = concatenate_gene_count_tables(input['all_gene_counts_df'], script_mode=params['script_mode'], join='outer')
        grouped_df.to_csv(output['gene_counts_df'], sep='\t', index_label = ['gene_id'])


'''
- 5. Estimate subtraction coefficients: ChIP - l * INPUT
'''

def chip_input_bam_files(wildcards):

    ## Attention! if `odf` raw files have been updated
    ## => Need to re-run:
    ##    - (1st) `sample_tables.ipynb` to obtain `chip_input_map.tsv`.
    ##    - (2nd) `sample_names.ipynb` to obtain `file_annotation.csv` needed for snakemake.
    
    ## First, try to use `INPUT_1`
    input_dataset_id = datasets_df[datasets_df['id'] == wildcards.dataset_id].INPUT_1
    
    ## If it does not exist, use `OIN`. There is at least one INPUT or OIN for each ChIP sample
    if not input_dataset_id.isnull().all(): 
      input_dataset_id = input_dataset_id.values[0]
    else:
      input_dataset_id = datasets_df[datasets_df['id'] == wildcards.dataset_id].OIN.values[0]
    
    # For INPUT subtraction is better to use regular .bam instead of (filtered) .tagged.bam
    input_bam = os.path.join(config["data_dir"], "seq_data/INPUT/bam", input_dataset_id, input_dataset_id + ".Aligned.sortedByCoord.out.bam")
    input_bam_bai = os.path.join(config["data_dir"], "seq_data/INPUT/bam", input_dataset_id, input_dataset_id + ".Aligned.sortedByCoord.out.bam.bai")
    
    return input_bam, input_bam_bai
    
checkpoint estimate_input_subtraction_factor:
    input:
        ## ChIP: .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        #bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"), # filtered bam
        ## ChIP: (indexed) .bam.bai file
        bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),
        #bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.tagged.bai"),  # filtered bam
        ## INPUT: .bam & (indexed) .bam.bai files
        chip_input_bam = chip_input_bam_files,
        
    output:
        input_factors = os.path.join(config["data_dir"], "seq_data/{seq_category}/INPUT_factors/{dataset_id}/INPUT_factors.csv"),

    params:
        ## Gene Information Table: GDF file
        gdf = config["gdf_file"],
        ## This script calculates Coverage (cvg) given an alignment file.
        input_factor_script = 'htseq/scripts/estimate_input_factor.py',
        out_dir = os.path.join(config["data_dir"], "seq_data/{seq_category}/INPUT_factors/{dataset_id}"),

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/5_estimate_input_factor/{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "seq_data/{seq_category}/INPUT_factors/log/{dataset_id}_estimate_input_factor.log")

    shell:
        "python {params.input_factor_script} "
        "--in_bam {input.bam} "
        "--in_chip_input {input.chip_input_bam[0]} "
        "--in_gdf {params.gdf} "
        "-o {params.out_dir} 2> {log}"
        

'''
- 5. Compute Coverage
'''

rule coverage_files:
    input:
        ## .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        ## (indexed) .bam.bai file
        bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),

    output:
        ## raw coverage files
        int_coverage = os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/int_{dataset_id}_coverage.wig"),
        frac_coverage = os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/frac_{dataset_id}_coverage.wig"),

    params:
        ## This script calculates Coverage (cvg) given an alignment .bam file.
        coverage_script = 'htseq/scripts/coverage.py',
        ## Reference Annotation: GTF/GFF file
        gff = config["gtf_file"],
        ## Parameters:
        ## - `seq_type`: distinguish between ChIP (Unstranded) and RNA (Stranded)
        seq_type = lambda wildcards: all_samples_df[all_samples_df['id'] == wildcards.dataset_id].pipeline_type.values[0],
        #seq_type = 'ChIP',
        out_dir = os.path.join(config["data_dir"], 'results/{seq_category}/coverage'),

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/5a_coverage_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/log/coverage.log")

    shell:
        "python {params.coverage_script} --in_bam {input.bam} --in_gtf {params.gff} --seq_type {params.seq_type} -o {params.out_dir} 2> {log}"

def get_input_factor(wildcards):
  
    #input_factor_file = os.path.join(config["data_dir"], "seq_data/{}/INPUT_factors/{}/INPUT_factors.csv".format(wildcards.seq_category, wildcards.dataset_id))
    input_factor_file =  checkpoints.estimate_input_subtraction_factor.get(**wildcards).output['input_factors']    
    input_factor_df = pd.read_csv(input_factor_file, sep ='\t', index_col = False)
    
    # distinguish between 'Pol II' and 'H3K9me2' ChIP
    if 'H3K9me2' in wildcards.dataset_id:
      
        # 95th Percentile
        def q95(x):
            return x.quantile(0.95)
        # 98th Percentile
        def q98(x):
            return x.quantile(0.98)
        
        input_factor_summary = input_factor_df.groupby('method').agg(
        # Get max of the 'input_factor' column for each group
            max_if =('input_factor', max),
            # Get min of the 'input_factor' column for each group
            min_if=('input_factor', min),
            # Get mean of the 'input_factor' column for each group
            mean_if=('input_factor', 'mean'),
            # Get median of the 'input_factor' column for each group
            median_if=('input_factor', 'median'),
            # Get standard deviation of the 'input_factor' column for each group
            std_if=('input_factor', 'std'),
            # Get 95th quantile of the 'input_factor' column for each group
            quantile_95_if=('input_factor', q95),
            # Get 98th quantile of the 'input_factor' column for each group
            quantile_98_if=('input_factor', q98),
            # Apply a lambda to date column
            #num_days=("date", lambda x: (max(x) - min(x)).days)    
            )
            
        # only use 'global lambdas' - tend to be more conservative (remove more signal)
        # => summarize with 'quantile_98_if' (extremely conservative!!)
        #input_factor = input_factor_summary.loc['global_lambda']['quantile_98_if'].round(3)  
        # => summarize with 'quantile_95_if' (slightly less but extremely conservative!!)
        #input_factor = input_factor_summary.loc['global_lambda']['quantile_95_if'].round(3)
        # => summarize with 'mean' (less conservative, but more than 'median')
        #input_factor = input_factor_summary.loc['global_lambda']['mean_if'].round(3)
        # => summarize with 'median' (even less conservative than 'mean')
        input_factor = input_factor_summary.loc['global_lambda']['median_if'].round(3)
        
        #import pdb; pdb.set_trace()
            
    else:
        
        # only use 'global lambdas' - tend to be more conservative (remove more signal)
        method_id = 'global_lambda'
        # only use 'count lambdas' 
        #method_id = 'count_lambda'
        input_factor = input_factor_df[input_factor_df['method'] == method_id]

        # summarize with 'median'
        #input_factor = 1.22 if wildcards.dataset_id == "WT_S2-ChIP_1" else 0.671
        input_factor = input_factor['input_factor'].median().round(3)
        #input_factor = input_factor['input_factor'].max().round(3)
        #input_factor = input_factor['input_factor'].max().round(3) * 3
        #input_factor = input_factor['input_factor'].max().round(3) * 3
        
    return input_factor

rule subtract_input_coverage_files:
    input:
        ## ChIP: .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        ## ChIP: (indexed) .bam.bai file
        bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),
        ## INPUT factor:
        input_factors = os.path.join(config["data_dir"], "seq_data/{seq_category}/INPUT_factors/{dataset_id}/INPUT_factors.csv"),
        ## INPUT: .bam & (indexed) .bam.bai files
        chip_input_bam = chip_input_bam_files,

    output:
        ## INPUT subtracted coverage files
        int_coverage = os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/int_subtracted_INPUT_{dataset_id}_coverage.wig"),
        frac_coverage = os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/frac_subtracted_INPUT_{dataset_id}_coverage.wig"),

    params:
        ## This script calculates Coverage (cvg) given an alignment .bam file.
        coverage_script = 'htseq/scripts/coverage.py',
        input_factor = get_input_factor,
        ## Reference Annotation: GTF/GFF file
        gff = config["gtf_file"],
        ## Parameters:
        ## - `seq_type`: distinguish between ChIP (Unstranded) and RNA (Stranded)
        seq_type = lambda wildcards: datasets_df[datasets_df['id'] == wildcards.dataset_id].pipeline_type.values[0],
        #seq_type = 'ChIP',
        out_dir = os.path.join(config["data_dir"], 'results/{seq_category}/coverage'),

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/5b_subtracted_INPUT_coverage_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/log/subtracted_INPUT_coverage.log")

    shell:
        "python {params.coverage_script} --in_bam {input.bam} --in_chip_input {input.chip_input_bam[0]} --norm_factor {params.input_factor} --in_gtf {params.gff} --seq_type {params.seq_type} -o {params.out_dir} 2> {log}"
        
        
'''
- 6. Compute (raw) and TPM-normed `Gene` & `Intron` Counts Tables
'''

## Differences with respect to RNA-pipeline
## - we use the (raw) .bam files.
## - ignore strandedness of reads alignments.
## - for ChIP samples:
##     - use whole 'gene' feature (INCLUDING intronic region).
##     - norm by 'gene_length'
## - also have to compute INPUT subtracted counts.

## --------
## A. Genes
## --------

rule gene_expression_table:
    input:
        ## .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        #bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"),  # filtered bam
        ## (indexed) .bam.bai file
        bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),
        #bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam.bai"),  # filtered bam

    output:
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"),
        ## - Summary of counting data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/summary_{dataset_id}_pombe_gene_count_matrix.csv")

    params:
        ## This script calculates Gene Expression Tables (gxt) given an alignment file(s).
        gxt_script = 'htseq/scripts/gene_counts.py', ## my script
        ## read alignments dir: (in case of ChIP not tagged)
        #bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/bam'),
        #bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/tagged_bam'),  # filtered bam
        ## Reference Annotation: GTF/GFF file
        #gff = config["gtf_file"],
        ## Gene Information Table: GDF file
        gdf = config["gdf_file"],
        ## Gene Counts dir
        gxt_dir = os.path.join(config["data_dir"], 'results/{seq_category}/xp_data'),
        ## -----------
        ## Parameters:
        ## -----------
        ## - `seq_type`: distinguish between ChIP (Unstranded, gene), RIP (Stranded, gene) and pA-RNA & total-RNA (Stranded, transcript)
        seq_type = lambda wildcards: all_samples_df[all_samples_df['id'] == wildcards.dataset_id].pipeline_type.values[0],
        ## - Feature type: (default uses features in 'gdf')
        #feature_type = expand("--feature_type {feature_type} ", feature_type = ['gene', 'region']),
        ## - Feature type: used for naming output_files
        feature_type = 'gene',
        ## ----------------
        ## Count Parameters
        ## ----------------
        ## - Potentially allow multimapped-reads - with less than max_nh (default:16) hits
        ## => We only allow this if the read maps to a repeat feature!
        max_nh = max_nh,
        ## - How to deal with overlap of READ and FEATURE?
        count_mode = count_mode, 
        #count_mode = 'union', # ("union", "intersection-strict", "intersection-nonempty")
        ## - How to deal with multiple-overlapping FEATURES? (according to count_mode)
        ambiguous_assignment_mode = ambiguous_assignment_mode,
        #ambiguous_assignment_mode = 'all', # ("none", "all", "fractional")
        ## - How to deal with multi-mapped reads?
        multimapped_mode = multimapped_mode, 
        #multimapped_mode = 'fractional', # ("none", "all", "fractional", "ignore_secondary"),
        ## - Prefix - with Count Parameters:
        #prefix = 'chip_' + prefix,
        #prefix = prefix,

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/6a_gene_expression_tables_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/log/gene_expression_tables.log")

    shell:
        ## htseq's script
        "python {params.gxt_script} "
        "--in_bam {input.bam} "
        #"--in_gtf {params.gff} "
        "--in_gdf {params.gdf} "
        "--seq_type {params.seq_type} "
        #"--prefix {params.prefix} "
        #"{params.feature_type}"
        "--feature_type {params.feature_type} "
        "--max_nh {params.max_nh} "
        "--count_mode {params.count_mode} "
        "--ambiguous_assignment_mode {params.ambiguous_assignment_mode} "
        "--multimapped_mode {params.multimapped_mode} "
        "-o {params.gxt_dir} 2> {log}"

rule group_gene_expression_tables:
    input:
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        ## Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/summary_{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
    
    output: 
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_gene_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_tpm_matrix.csv"),
        ## Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/summary_pombe_tpm_matrix.csv"),

    params:
        script_mode = 'htseq',
        #script_mode = 'parastou',
        ## Gene Information Table GDF file
        gdf = config["gdf_file"],
        ## Column to select from Data Frame:
        count_col = 'count',

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/6a_gene_expression_table.txt", N_BENCHMARKS)

    run:
        ## Init lists of individual sample DataFrames
        list_df = []
        
        ## ---------------------------
        ## Load Gene Information Table: gdf
        ## ---------------------------
        
        gdf = pd.read_csv(params['gdf'], sep='\t')
        
        gdf_header = ['seqid', 'type', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'cds_length', 'utr_length', 'intron_length', 'gene_length', 'category', 'bio_type']
        # keep only rows of interest
        gdf = gdf[gdf_header]
        gdf = gdf.set_index('gene_id')
        gdf = gdf.sort_index()

        # add gdf to list of DataFrames
        list_df.append(gdf)
        
        ## -------------------------
        ## A. Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        ## -------------------------

        ## Concatenate DataFrames:
        grouped_df = concatenate_dataframes(input['summary_tpm_gxt'], count_col = params['count_col'])
        grouped_df.to_csv(output['summary_tpm_gxt'], sep='\t')

        ## ------------------------
        ## B. raw Gene Counts Data: (gxt)
        ## ------------------------

        ## Concatenate DataFrames:
        #grouped_df = concatenate_dataframes(input['gxt'], list_df, assert_size=len(gdf.index), script_mode=params['script_mode'])
        grouped_df = concatenate_dataframes(input['gxt'], list_df, script_mode=params['script_mode'], count_col=params['count_col'])
        grouped_df.to_csv(output['gxt'], sep='\t', index_label = ['gene_id'] )

        ## ------------------------------
        ## C. TPM-normed Gene Counts Data: (tpm_gxt)
        ## -------------------------------

        ## Concatenate DataFrames
        #grouped_df = concatenate_dataframes(input['tpm_gxt'], list_df, assert_size=len(gdf.index), script_mode=params['script_mode'])
        grouped_df = concatenate_dataframes(input['tpm_gxt'], list_df, script_mode=params['script_mode'], count_col=params['count_col'])
        grouped_df.to_csv(output['tpm_gxt'], sep='\t', index_label = ['gene_id'])

def chip_input_tagged_bam_files(wildcards):
    
    ## First, try to use `INPUT_1`
    input_dataset_id = datasets_df[datasets_df['id'] == wildcards.dataset_id].INPUT_1
    
    ## If it does not exist, use `OIN`. There is at least one INPUT or OIN for each ChIP sample
    if not input_dataset_id.isnull().all(): 
      input_dataset_id = input_dataset_id.values[0]
    else:
      input_dataset_id = datasets_df[datasets_df['id'] == wildcards.dataset_id].OIN.values[0]
    
    input_bam = os.path.join(config["data_dir"], "seq_data/INPUT/tagged_bam", input_dataset_id, input_dataset_id + ".Aligned.sortedByCoord.out.tagged.bam")  # filtered bam
    input_bam_bai = os.path.join(config["data_dir"], "seq_data/INPUT/tagged_bam", input_dataset_id, input_dataset_id + ".Aligned.sortedByCoord.out.tagged.bam.bai") # filtered bam
    
    #import pdb; pdb.set_trace()

    return input_bam, input_bam_bai
    
def chip_input_count_file(wildcards):
    
    ## First, try to use `INPUT_1`
    input_dataset_id = datasets_df[datasets_df['id'] == wildcards.dataset_id].INPUT_1
    
    ## If it does not exist, use `OIN`. There is at least one `INPUT` or `OIN` for each ChIP sample
    if not input_dataset_id.isnull().all(): 
      input_dataset_id = input_dataset_id.values[0]
    else:
      input_dataset_id = datasets_df[datasets_df['id'] == wildcards.dataset_id].OIN.values[0]
    
    count_file = os.path.join(config["data_dir"], "results/INPUT/xp_data", input_dataset_id, input_dataset_id + "_pombe_gene_count_matrix.csv")

    #import pdb; pdb.set_trace()

    return count_file

# NOTE: Version with INPUT subtraction
rule subtract_input_gene_expression_table:
    input:
        ## INPUT factor:
        input_factors = os.path.join(config["data_dir"], "seq_data/{seq_category}/INPUT_factors/{dataset_id}/INPUT_factors.csv"),
        ## A. ChIP: count file
        chip_gxp = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"),
        ## B. INPUT: count file
        input_gxp = chip_input_count_file
        # ## A. ChIP:.bam & (indexed) .bam.bai files
        # # .bam file
        # bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        # #bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"),  # filtered bam
        # # (indexed) .bam.bai file
        # bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),
        # #bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam.bai"),  # filtered bam
        # ## B. INPUT: .bam & (indexed) .bam.bai files
        # chip_input_bam = chip_input_bam_files,
        # #chip_input_bam = chip_input_tagged_bam_files, # filtered bam
        
    output:
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_subtracted_INPUT_pombe_gene_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_subtracted_INPUT_pombe_tpm_matrix.csv"),

    params:
        ## This script calculates Gene Expression Tables (gxt) given an alignment file(s).
        gxt_script = 'htseq/scripts/gene_counts.py', ## my script
        ## read alignments dir: (in case of ChIP not tagged)
        #bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/bam'),
        #bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/tagged_bam'), # filtered bam
        ## Reference Annotation: GTF/GFF file
        #gff = config["gtf_file"],
        ## Gene Information Table: GDF file
        gdf = config["gdf_file"],
        ## Gene Counts dir
        gxt_dir = os.path.join(config["data_dir"], 'results/{seq_category}/xp_data'),
        ## -----------
        ## Parameters:
        ## -----------
        ## - `seq_type`: distinguish between ChIP (unstranded, gene), RIP (Stranded, gene) and pA-RNA & total-RNA (Stranded, transcript)
        seq_type = lambda wildcards: all_samples_df[all_samples_df['id'] == wildcards.dataset_id].pipeline_type.values[0],
        ## - Feature type: (default uses features in 'gdf')
        #feature_type = expand("--feature_type {feature_type} ", feature_type = ['gene', 'region']),
        ## - Feature type: used for naming output_files
        feature_type = 'gene',
        input_factor = get_input_factor,
        ## - Prefix - with Count Parameters:
        #prefix = 'chip_' + prefix,
        #prefix = prefix,
        ## ----------------
        ## Count Parameters
        ## ----------------
        ## - Potentially allow multimapped-reads - with less than max_nh (default:16) hits
        ## => We only allow this if the read maps to a repeat feature!
        max_nh = max_nh,
        ## - How to deal with overlap of READ and FEATURE?
        count_mode = count_mode, 
        #count_mode = 'union', # ("union", "intersection-strict", "intersection-nonempty")
        ## - How to deal with multiple-overlapping FEATURES? (according to count_mode)
        ambiguous_assignment_mode = ambiguous_assignment_mode,
        #ambiguous_assignment_mode = 'all', # ("none", "all", "fractional")
        ## - How to deal with multi-mapped reads?
        multimapped_mode = multimapped_mode, 
        #multimapped_mode = 'fractional', # ("none", "all", "fractional", "ignore_secondary"),

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/6b_subtracted_INPUT_gene_expression_tables_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/log/subtracted_INPUT_gene_expression_tables.log")

    shell:
        ## htseq's script
        "python {params.gxt_script} "
        #"--in_bam {input.bam} "
        #"--in_chip_input {input.chip_input_bam[0]} "
        "--in_count {input.chip_gxp} "
        "--in_chip_input_count {input.input_gxp} "
        "--norm_factor {params.input_factor} "
        #"--in_gtf {params.gff} "
        "--in_gdf {params.gdf} "
        "--seq_type {params.seq_type} "
        #"--prefix {params.prefix} "
        #"{params.feature_type}"
        "--feature_type {params.feature_type} "
        "--max_nh {params.max_nh} "
        "--count_mode {params.count_mode} "
        "--ambiguous_assignment_mode {params.ambiguous_assignment_mode} "
        "--multimapped_mode {params.multimapped_mode} "
        "-o {params.gxt_dir} 2> {log}"

# NOTE: Need to modify input files in `simulated_data` mode: (or `H3K9me2`)
# ChIP: [S2-ChIP, S5-ChIP]
rule chip_gene_expression_tables:
    input:
        ## raw Gene Counts Table (gxt) as .csv file
        #gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']), # for simulated-data
        #gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']), # without INPUT subtraction
        gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_subtracted_INPUT_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        ## TPM-normed Gene Expression Table (tpm_gxt) as .csv file
        #tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']), # for simulated-data
        #tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']), # without INPUT subtraction
        tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_subtracted_INPUT_pombe_tpm_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        ## Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        #summary_tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/summary_{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']), # for simulated-data
        #summary_tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/summary_{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']), # without INPUT subtraction
        summary_tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/summary_{dataset_id}_subtracted_INPUT_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),

    output: 
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/xp_data/ChIP/chip_pombe_gene_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/xp_data/ChIP/chip_pombe_tpm_matrix.csv"),
        ## Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = os.path.join(config["data_dir"], "results/xp_data/ChIP/chip_summary_pombe_tpm_matrix.csv"),

    params:
        script_mode = 'htseq',
        #script_mode = 'parastou',
        ## Gene Information Table GDF file
        gdf = config["gdf_file"],
        ## Column to select from Data Frame:
        count_col = 'count',

    benchmark:
        repeat("benchmarks/1_preprocess/6_chip_gene_expression_table.txt", N_BENCHMARKS)

    run:
        ## Init lists of individual sample DataFrames
        list_df = []
        
        ## ---------------------------
        ## Load Gene Information Table: gdf
        ## ---------------------------
        
        gdf = pd.read_csv(params['gdf'], sep='\t')
        
        gdf_header = ['seqid', 'type', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'cds_length', 'utr_length', 'intron_length', 'gene_length', 'category', 'bio_type']
        # keep only rows of interest
        gdf = gdf[gdf_header]
        gdf = gdf.set_index('gene_id')
        gdf = gdf.sort_index()

        # add gdf to list of DataFrames
        list_df.append(gdf)
        
        ## -------------------------
        ## A. Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        ## -------------------------
        
        ## Concatenate DataFrames:
        grouped_df = concatenate_dataframes(input['summary_tpm_gxt'], count_col=params['count_col'])
        grouped_df.to_csv(output['summary_tpm_gxt'], sep='\t')
        
        ## ------------------------
        ## B. raw Gene Counts Data: (gxt)
        ## ------------------------

        ## Concatenate DataFrames:
        #grouped_df = concatenate_dataframes(input['gxt'], list_df, assert_size=len(gdf.index), script_mode=params['script_mode'])
        grouped_df = concatenate_dataframes(input['gxt'], list_df, script_mode=params['script_mode'], count_col=params['count_col'])
        grouped_df.to_csv(output['gxt'], sep='\t', index_label=['gene_id'] )

        ## ------------------------------
        ## C. TPM-normed Gene Counts Data: (tpm_gxt)
        ## -------------------------------

        ## Concatenate DataFrames
        #grouped_df = concatenate_dataframes(input['tpm_gxt'], list_df, assert_size=len(gdf.index), script_mode=params['script_mode'])
        grouped_df = concatenate_dataframes(input['tpm_gxt'], list_df, script_mode=params['script_mode'], count_col=params['count_col'])
        grouped_df.to_csv(output['tpm_gxt'], sep='\t', index_label=['gene_id'])
        
        
## ----------
## B. Introns
## ----------

rule intron_expression_table:
    input:
        ## .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        ## (indexed) .bam.bai file
        bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),

    output:
        ## raw Intron Counts Table (gxt) and TPM-normed Intron Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_intron_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"),
        ## - Summary of counting data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/summary_{dataset_id}_pombe_intron_count_matrix.csv")

    params:
        ## This script calculates Intron Expression Tables (gxt) given an alignment file(s).
        gxt_script = 'htseq/scripts/gene_counts.py', ## my script
        ## read alignments dir: (in case of ChIP not tagged)
        bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/bam'),
        ## Intron Information Table: GDF file
        #gdf = config["gdf_file"],
        gdf = config["intron_df_file"],
        ## Intron Counts dir
        gxt_dir = os.path.join(config["data_dir"], 'results/{seq_category}/intron_xp_data'),
        ## -----------
        ## Parameters:
        ## -----------
        ## - `seq_type`: distinguish between ChIP (Unstranded, gene), RIP (Stranded, gene) and pA-RNA & total-RNA (Stranded, transcript)
        seq_type = lambda wildcards: all_samples_df[all_samples_df['id'] == wildcards.dataset_id].pipeline_type.values[0],
        ## - Feature type: used for naming output_files
        feature_type = 'intron',
        ## - Prefix - with Count Parameters:
        #prefix = 'chip_' + prefix,
        #prefix = prefix,
        ## ----------------
        ## Count Parameters
        ## ----------------
        ## - Potentially allow multimapped-reads - with less than max_nh (default:16) hits
        ## => We only allow this if the read maps to a repeat feature!
        max_nh = max_nh,
        ## - How to deal with overlap of READ and FEATURE?
        count_mode = count_mode,
        #count_mode = 'intersection-strict', # ("union", "intersection-strict", "intersection-nonempty")
        ## - How to deal with reads multiple-overlapping FEATURES? (according to count_mode)
        ambiguous_assignment_mode = ambiguous_assignment_mode, 
        #ambiguous_assignment_mode = 'all', # ("none", "all", "fractional")
        ## - How to deal with multi-mapped reads?
        multimapped_mode = multimapped_mode, 
        #multimapped_mode = 'fractional', # ("none", "all", "fractional", "ignore_secondary")


    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/6b_intron_expression_tables_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/log/intron_expression_tables.log")

    shell:
        ## htseq's script
        "python {params.gxt_script} "
        "--in_bam {input.bam} "
        "--in_gdf {params.gdf} "
        "--seq_type {params.seq_type} "
        "--feature_type {params.feature_type} "
        "--max_nh {params.max_nh} "
        "--count_mode {params.count_mode} "
        "--ambiguous_assignment_mode {params.ambiguous_assignment_mode} "
        "--multimapped_mode {params.multimapped_mode} "
        "-o {params.gxt_dir} 2> {log}"

rule group_intron_expression_tables:
    input:
        ## raw Intron Counts Table (gxt) and TPM-normed Intron Expression Table (tpm_gxt) as .csv files
        gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_intron_count_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
        ## Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/summary_{dataset_id}_pombe_intron_count_matrix.csv"), zip, seq_category=all_samples_df['seq_category'], dataset_id=all_samples_df['id']),
    
    output: 
        ## raw Intron Counts Table (gxt) and TPM-normed Intron Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/pombe_intron_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/pombe_tpm_matrix.csv"),
        ## Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/summary_pombe_tpm_matrix.csv"),

    params:
        script_mode = 'htseq',
        #script_mode = 'parastou',
        ## Intron Information Table GDF file
        #gdf = config["gdf_file"],
        gdf = config["intron_df_file"],
        ## Column to select from Data Frame:
        count_col = 'count',

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/6b_intron_expression_table.txt", N_BENCHMARKS)

    run:
        ## Init lists of individual sample DataFrames
        list_df = []
        
        ## ---------------------------
        ## Load Gene Information Table: gdf
        ## ---------------------------
        
        gdf = pd.read_csv(params['gdf'], sep='\t', index_col = 'gene_id')
        gdf = gdf.sort_index()

        ## Add gdf to list of DataFrames
        list_df.append(gdf)
        
        ## -------------------------
        ## A. Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        ## -------------------------

        ## Concatenate DataFrames:
        grouped_df = concatenate_dataframes(input['summary_tpm_gxt'], count_col = params['count_col'])
        grouped_df.to_csv(output['summary_tpm_gxt'], sep='\t')

        ## ------------------------
        ## B. raw Gene Counts Data: (gxt)
        ## ------------------------

        ## Concatenate DataFrames:
        #grouped_df = concatenate_dataframes(input['gxt'], list_df, assert_size=len(gdf.index), script_mode=params['script_mode'])
        grouped_df = concatenate_dataframes(input['gxt'], list_df, script_mode=params['script_mode'], count_col= params['count_col'])
        grouped_df.to_csv(output['gxt'], sep='\t', index_label = ['gene_id'] )

        ## ------------------------------
        ## C. TPM-normed Gene Counts Data: (tpm_gxt)
        ## -------------------------------

        ## Concatenate DataFrames
        #grouped_df = concatenate_dataframes(input['tpm_gxt'], list_df, assert_size=len(gdf.index), script_mode=params['script_mode'])
        grouped_df = concatenate_dataframes(input['tpm_gxt'], list_df, script_mode=params['script_mode'], count_col = params['count_col'])
        grouped_df.to_csv(output['tpm_gxt'], sep='\t', index_label = ['gene_id'])
   
