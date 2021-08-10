"""
Author: P. Monteagudo
Affiliation: LMU, GC
Aim:    Snakefile containing rules for pA-RNA-seq (PolyA mRNA), RIP-seq
        (nascent RNA) and regular RNA-seq (total RNA) Data Analysis Pipe-line.

Date: November 27  2019
Run:    snakemake -s preprocess_RNA.smk -n
        snakemake -s preprocess_RNA.smk -p
        snakemake -s preprocess_RNA.smk -p -j 20 -R gene_expression_table

        # Touch output files (mark them up to date without really changing them) instead of running their commands.
        snakemake -s preprocess_RNA.smk -t

        # Number of Jobs
        snakemake -s preprocess_RNA.smk -p -j 20
        # DAG Jobs
        snakemake -s preprocess_RNA.smk --dag | dot -Tpdf > jobs_dag.pdf
        # DAG Rules
        snakemake -s preprocess_RNA.smk --rulegraph | dot -Tsvg > rules_dag.svg
        ## Go on with independent jobs if a job fails.
        snakemake -s preprocess_RNA.smk -p -j 20 -k
        ## Force the re-execution or creation of the given rules or files.
        snakemake -s preprocess_RNA.smk -p -j 20 -R construct_TI_models visualize_TI_models

        ## Cluster
        snakemake -s preprocess_RNA.smk -p --cluster "sbatch" --jobs 4
        snakemake -s preprocess_RNA.smk -p --cluster "sbatch --mem=2G --partition=hi_mem --nodelist=node3 --output=slurm_job_files/slurm-%j.out" --jobs 100
        snakemake -s preprocess_RNA.smk -p --cluster "sbatch --mem=2G --partition=full_high --nodelist=master --output=slurm_job_files/slurm-%j.out" --jobs 1
        snakemake -s preprocess_RNA.smk -p --cluster "sbatch --mem=2G --partition=full_high --output=slurm_job_files/slurm-%j.out" --jobs 1000 -k
        snakemake -s preprocess_RNA.smk -p --cluster "sbatch --mem=2G --time=0-03:00:00 --partition=full_high --output=slurm_job_files/slurm-%j.out" --jobs 500 -k -R unzip_fastq_files

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
#ambiguous_assignment_mode = 'fractional'
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

# use `simulated_data` mode                 
simulated_data = False
#simulated_data = True
                 
## ---------------
## - RNA Datasets: | S2-RIP | S5-RIP | pA-RNA | total-RNA |
## ---------------

# import DataFrame containing `sample` annotation and filter for samples that belong to the RNA Pipeline
datasets_df = pd.read_csv(config["sample_annotation_file"], sep = '\t')

if not simulated_data:
  datasets_df = datasets_df[datasets_df['pipeline_type'] == 'RNA']
else:
  datasets_df = datasets_df[(datasets_df['pipeline_type'] == 'simulated-data') & (datasets_df['mutant_id'].str.contains("rna")) ]
  datasets_df['pipeline_type'] = 'RNA'

# Ignore S5-samples for now!
#datasets_df = datasets_df[datasets_df['seq_type'] != 'S5-RIP']

# rename 'sample_id' column to 'id'
datasets_df = datasets_df.rename(columns = {'sample_id' : 'id'})
#datasets_df = datasets_df.iloc[0:2, :]

## -----------------
## - Select Data Set
## -----------------
#import pdb; pdb.set_trace()

# by `dataset_id`
select_dataset = ""
#select_dataset = "80_pA-RNA_1"
#select_dataset = "WT_pA-RNA_1"
#select_dataset = "rna-fake-reads_simulated-data_2"
#datasets_df = datasets_df[datasets_df.id == select_dataset]

# by `dataset_id`s
select_datasets = ["80_pA-RNA_1", "WT_pA-RNA_1"]
#all_samples_df = datasets_df[datasets_df.id.isin(select_datasets)]

# by `mutant_id`
select_mutant = "WT_"
#select_mutant = "80_"
#datasets_df = datasets_df[datasets_df['id'].str.contains(select_mutant)]

## ------------------
## - Ignore Data Sets
## ------------------

# EXITING because of fatal ERROR: not enough memory for BAM sorting:
# SOLUTION: re-run STAR with at least --limitBAMsortRAM 2821915860
ignore_datasets = []
#ignore_datasets = ['302_S2RIP_2']
datasets_df = datasets_df[~datasets_df.id.isin(ignore_datasets)]

## --------------------
## - Truncate Data Sets
## --------------------

## Error - EXITING because of FATAL ERROR in reads input: short read sequence line: 1
## Need to trim reads:
# truncate_datasets = ['302_S2RIP_2']
# truncate_datasets.extend(['283_RNA_pA_4', '301_RNA_pA_3', '301_S2RIP_3',
#                         '302_S2RIP_3', '324_RNA_pA_3', '324_S2RIP_3', '491_S2RIP_3',
#                         '504_RNA_pA_1', '504_RNA_pA_2', '530_RNA_pA_1',
#                         '530_RNA_pA_2', '638_RNA_pA_1', '638_RNA_pA_2',
#                         '63_RNA_pA_3',  '63_RNA_pA_4', '63_S2RIP_2'])
truncate_datasets = datasets_df[datasets_df['trimmed']].id.values
#datasets_df = datasets_df[datasets_df.id.isin(truncate_datasets)]

#import pdb; pdb.set_trace()


#######################################
# Convienient rules to define targets #
#######################################

## Mark a rule as local, so that it is not submitted to the cluster and instead executed on the host node
localrules: all

rule all:
    input:
        # config["rna_star_idx"],
        # ## ---------------------
        # ## 1. Unzip fastq files:
        # ## ---------------------
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/fastq/{dataset_id}.fastq"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # ## ---------------------
        # ## 2. Align fastq files: using `STAR`
        # ## ---------------------
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.out.bam"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # ## -----------------------------------
        # ## 3.  Sort and Index read alignments: using `samtools`
        # ## -----------------------------------
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # ## -------------------------------------------
        # ## 4. Remove rRNA/Tag bam files with gene-id's
        # ## -------------------------------------------
        # # by sample: tagged/filtered bam files
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # ## summarize samples by 'seq_category': (naive) ribosomal RNA and other RNA features Counts as .csv files:
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/ribosomal_rna_pombe_gene_count_matrix.csv"), seq_category=datasets_df['seq_category'].unique()),
        # expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/rna_pombe_gene_count_matrix.csv"), seq_category=datasets_df['seq_category'].unique()),
        ## -------------------
        ## 5. Compute Coverage
        ## -------------------
        # RNA: Coverage `+` strand
        expand(os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/int_plus_{dataset_id}_coverage.wig"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/frac_plus_{dataset_id}_coverage.wig"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # RNA: Coverage `-` strand
        expand(os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/int_minus_{dataset_id}_coverage.wig"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/frac_minus_{dataset_id}_coverage.wig"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        ## -------------------------------------------------
        ## 6. Compute (raw) and TPM-normed Gene Count Tables
        ## -------------------------------------------------
        ## ---------
        ## 6a. Genes
        ## ---------
        # by sample: raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # summarize samples by 'seq_category':
        #expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_gene_read_count_matrix.csv"), seq_category=datasets_df['seq_category']),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_gene_count_matrix.csv"), seq_category=datasets_df['seq_category'].unique()),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_tpm_matrix.csv"), seq_category=datasets_df['seq_category'].unique()),
        # summarize all `RNA` samples:
        os.path.join(config["data_dir"], "results/xp_data/RNA/rna_pombe_gene_count_matrix.csv"),
        os.path.join(config["data_dir"], "results/xp_data/RNA/rna_pombe_tpm_matrix.csv"),
        ## -----------
        ## 6b. Introns
        ## -----------
        # # by sample: raw Intron Counts Table (gxt) and TPM-normed Intron Expression Table (tpm_gxt) as .csv files
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_intron_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        # summarize samples by 'seq_category': raw Intron Counts Table (gxt) and TPM-normed Intron Expression Table (tpm_gxt) as .csv files
        #expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/pombe_intron_read_count_matrix.csv"), seq_category=datasets_df['seq_category']),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/pombe_intron_count_matrix.csv"), seq_category=datasets_df['seq_category'].unique()),
        expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/pombe_tpm_matrix.csv"), seq_category=datasets_df['seq_category'].unique()),

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

# Attention! might not be automatically up-dated when switching from one gtf version to a new one! (manually update)

rule star_index_rna:
    input:
        # reference FASTA file
        #fa = os.path.join(config["genome_dir"], 'spombe_v2/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa'),
        fa = config["fasta_file"],
        # reference GTF/GFF file
        #gff = os.path.join(config["annotation_dir"], 'gff_v2/Schizosaccharomyces_pombe_all_chromosomes.extended.gff3'),
        gff = config["gtf_file"],

    output:
        ## diffferent from ChIP analysis!!
        #directory(os.path.join(config["genome_dir"], "spombe_v2/star_idx"))
        directory(config["rna_star_idx"])

    threads: 20 # set the maximum number of available cores

    shell:
        'mkdir {output} && '
        '{config[star]} --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fa} '
        '--sjdbGTFfile {input.gff} '
        '--sjdbGTFfeatureExon CDS '
        '--sjdbGTFtagExonParentTranscript Parent ' ## necessary for gff
        #'--genomeSAindexNbases 11 ' ## small genomes: spombe ~ 12.82 Mb, min(14, log2(GenomeLength/2) - 1) ~ 11
        ## WARNING: --genomeSAindexNbases 14 is too large for the genome size=13893632, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 10
        '--genomeSAindexNbases 10 ' ## small genomes: spombe ~ 12.82 Mb, min(14, log2(GenomeLength/2) - 1) ~ 11
        ## WARNING: --genomeSAindexNbases 11 is too large for the genome size=13893632, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 10
        #'--readMapNumber 100000' ##  map the first ~100,000 reads
        '--sjdbOverhang 49' ## should we use another value, eg. length read - 1 = 50/49?
        #' 2> {log}'

rule align_reads:
    input:
        ## STAR index
        #ref_dir = os.path.join(config["genome_dir"], "spombe_v2/star_idx"),
        ref_dir = config["rna_star_idx"],
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
        ## Attention! for RNA-seq we run the alignment without the --alignIntronMax 1 arguments!!,
        ## => we still use --alignEndsType EndToEnd, to force end-to-end read alignment, do not soft-clip!
        '{config[star]} --runThreadN {threads} '
        '--genomeDir {input.ref_dir} '
        '--readFilesIn {input.fastq} '
        '--outFileNamePrefix {params.bam_prefix} '
        '--outSAMtype BAM Unsorted '  # sort manually with samtools, due to error
        #'--outSAMtype BAM SortedByCoordinate '
        ## max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped.
        '--outFilterMultimapNmax {params.max_nh} ' 
        #'--outFilterMultimapNmax 16 ' 
        #'--outFilterMultimapNmax 25 '
        #'--outFilterMultimapNmax 99 '
        ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC275578/
        '--alignIntronMin 20 ' ## minimun intron length, 20 nucleotides shortest allowed intron length for S.pombe.
        '--alignIntronMax 2000 ' ## maximun intron length, 2000 nucleotides as the maximum acceptable intron length for S.pombe,
                                ## because this is about twice the size of the longest identified intron in the fungus S.cerevisiae (20)
        '--alignEndsType EndToEnd '
        #'2> {log}'


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
        repeat("benchmarks/1_preprocess/{seq_category}/3b_index_aligned_reads/{dataset_id}.txt", N_BENCHMARKS)

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

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/4a_remove_rrna_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/log/remove_rrna.log")

    shell:
        "python {params.rem_rrna_script} --in_bam {input.bam} --in_gdf {params.gdf} -o {params.tagged_bam_dir} 2> {log}"

rule group_gene_count_tables:
    input:
        ## Ribosomal RNA (rRNA) and mRNA Gene Counts computed when parsing BAM file as .csv files
        all_rrna_counts_df = expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/rrna_counts.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        all_gene_counts_df = expand(os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/gene_counts.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),

    output:
        ## raw and TPM-normed Gene Counts Data (gxt)  as .csv files
        rrna_counts_df = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/ribosomal_rna_pombe_gene_count_matrix.csv"),
        gene_counts_df = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/rna_pombe_gene_count_matrix.csv"),

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
- 5. Compute Coverage
'''

## Attention! Need to decide wether to use `.bam` or `.tagged.bam`
## - If we use .tagged.bam, inter-gene regions will be empty!
## - If we use (raw) .bam, coverage files won't represent exactly what is used in counts!
## => I think we should use (raw) .bam

rule coverage_files:
    input:
        ## .bam file
        bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        #bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"),
        ## (indexed) .bam.bai file
        bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),
        #bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam.bai"),

    output:
        ## raw coverage files
        int_plus_coverage = os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/int_plus_{dataset_id}_coverage.wig"),
        frac_plus_coverage = os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/frac_plus_{dataset_id}_coverage.wig"),
        ## raw coverage files
        int_minus_coverage = os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/int_minus_{dataset_id}_coverage.wig"),
        frac_minus_coverage = os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/frac_minus_{dataset_id}_coverage.wig"),

    params:
        ## This script calculates Coverage (cvg) given an alignment file.
        coverage_script = 'htseq/scripts/coverage.py',
        ## Reference Annotation: GTF/GFF file
        gff = config["gtf_file"],
        ## Parameters:
        ## - `seq_type`: distinguish between ChIP (Unstranded) and RNA (Stranded)
        seq_type = lambda wildcards: datasets_df[datasets_df['id'] == wildcards.dataset_id].pipeline_type.values[0],
        #seq_type = 'RNA',
        out_dir = os.path.join(config["data_dir"], 'results/{seq_category}/coverage'),

    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/5_coverage_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "results/{seq_category}/coverage/{dataset_id}/log/coverage.log")

    shell:
        "python {params.coverage_script} --in_bam {input.bam} --in_gtf {params.gff} --seq_type {params.seq_type} -o {params.out_dir} 2> {log}"


'''
- 6. Compute (raw) and TPM-normed Gene Count tables
'''

## Differences with respect to ChIP-pipeline:
## - we use the tagged.bam files.
## - use strandedness of reads alignments.
## - for pA-RNA & total-RNA samples:
##     - use 'transcript' feature (EXCLUDING intronic region).
##     - norm by 'transcript_length'
## - for RIP:
##     - use whole 'gene' feature (INCLUDING intronic region).
##     - norm by 'gene_length'

## --------
## A. Genes
## --------

rule gene_expression_table:
    input:
        ## (tagged) .bam file
        tagged_bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"),
        #tagged_bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        ## tagged (indexed) .bam.bai file
        tagged_bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam.bai"),
        #tagged_bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),

    output:
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"),
        ## - Summary of counting data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/summary_{dataset_id}_pombe_gene_count_matrix.csv"),

    params:
        ## This script calculates Gene Expression Tables (gxt) given an alignment file(s).
        gxt_script = 'htseq/scripts/gene_counts.py', ## my script
        ## read alignments dir: (in case of RNA tagged_bams)
        tagged_bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/tagged_bam'),
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
        seq_type = lambda wildcards: datasets_df[datasets_df['id'] == wildcards.dataset_id].pipeline_type.values[0],
        #prefix = 'rna_',
        ## - Feature type: (default uses features in 'gdf')
        #feature_type = expand("--feature_type {feature_type} ", feature_type = ['gene', 'region']),
        ## - Feature type: used for naming output_files
        feature_type = 'gene',
        ## - Prefix - with Count Parameters:
        #prefix = 'rna_' + prefix,
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
        ## - How to deal with reads multiple-overlapping FEATURES? (according to count_mode)
        ambiguous_assignment_mode = ambiguous_assignment_mode,
        #ambiguous_assignment_mode = 'all', # ("none", "all", "fractional")
        ## - How to deal with multi-mapped reads?
        multimapped_mode = multimapped_mode,
        #multimapped_mode = 'all', # ("none", "all", "fractional", "ignore_secondary")


    benchmark:
        repeat("benchmarks/1_preprocess/{seq_category}/5a_gene_expression_tables_{dataset_id}.txt", N_BENCHMARKS)

    log:
        os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/log/gene_expression_tables.log")

    shell:
        ## htseq's script
        "python {params.gxt_script} "
        "--in_bam {input.tagged_bam} "
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
        gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        ## - Summary of counting data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/summary_{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
    
    output:
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_gene_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/xp_data/pombe_tpm_matrix.csv"),
        ## - Summary of counting data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
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
        grouped_df = concatenate_dataframes(input['summary_tpm_gxt'], count_col= params['count_col'])
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
        grouped_df = concatenate_dataframes(input['tpm_gxt'], list_df, script_mode=params['script_mode'], count_col= params['count_col'])
        grouped_df.to_csv(output['tpm_gxt'], sep='\t', index_label = ['gene_id'])

# RNA: [S2-RIP, S5-RIP, pA-RNA, total-RNA samples]
rule rna_gene_expression_tables:
    input:
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        ## - Summary of counting data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/xp_data/{dataset_id}/summary_{dataset_id}_pombe_gene_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
    
    output:
        ## raw Gene Counts Table (gxt) and TPM-normed Gene Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/xp_data/RNA/rna_pombe_gene_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/xp_data/RNA/rna_pombe_tpm_matrix.csv"),
        ## - Summary of counting data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = os.path.join(config["data_dir"], "results/xp_data/RNA/rna_summary_pombe_tpm_matrix.csv"),

    params:
        script_mode = 'htseq',
        #script_mode = 'parastou',
        ## Gene Information Table GDF file
        gdf = config["gdf_file"],
        ## Column to select from Data Frame:
        count_col = 'count',

    benchmark:
        repeat("benchmarks/1_preprocess/6_rna_gene_expression_table.txt", N_BENCHMARKS)

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
        ## (tagged) .bam file
        tagged_bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam"),
        #tagged_bam = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam"),
        ## tagged (indexed) .bam.bai file
        tagged_bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/tagged_bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.tagged.bam.bai"),
        #tagged_bam_bai = os.path.join(config["data_dir"], "seq_data/{seq_category}/bam/{dataset_id}/{dataset_id}.Aligned.sortedByCoord.out.bam.bai"),

    output:
        ## raw Intron Counts Table (gxt) and TPM-normed Intron Expression Table (tpm_gxt) as .csv files
        gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_intron_count_matrix.csv"),
        tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"),
        ## - Summary of counting data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/summary_{dataset_id}_pombe_intron_count_matrix.csv")

    params:
        ## This script calculates Intron Expression Tables (gxt) given an alignment file(s).
        gxt_script = 'htseq/scripts/gene_counts.py', ## my script
        ## read alignments dir: (in case of RNA tagged_bams)
        tagged_bam_dir = os.path.join(config["data_dir"], 'seq_data/{seq_category}/tagged_bam'),
        ## Intron Information Table: GDF file
        #gdf = config["gdf_file"],
        gdf = config["intron_df_file"],
        ## Intron Counts dir
        gxt_dir = os.path.join(config["data_dir"], 'results/{seq_category}/intron_xp_data'),
        ## -----------
        ## Parameters:
        ## -----------
        ## - `seq_type`: distinguish between ChIP (Unstranded, gene), RIP (Stranded, gene) and pA-RNA & total-RNA (Stranded, transcript)
        seq_type = lambda wildcards: datasets_df[datasets_df['id'] == wildcards.dataset_id].pipeline_type.values[0],
        ## - Feature type: used for naming output_files
        feature_type = 'intron',
        ## - Prefix - with Count Parameters:
        #prefix = 'rna_' + prefix,
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
        "--in_bam {input.tagged_bam} "
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
        gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_intron_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/{dataset_id}_pombe_tpm_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
        ## Summary of Count Data: [_total, _alignment_not_unique, _no_feature, _ambiguous]
        summary_tpm_gxt = expand(os.path.join(config["data_dir"], "results/{seq_category}/intron_xp_data/{dataset_id}/summary_{dataset_id}_pombe_intron_count_matrix.csv"), zip, seq_category=datasets_df['seq_category'], dataset_id=datasets_df['id']),
    
    ## need to be in '/results' in the project_dir otherwise jupyter can't access the files.
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
        repeat("benchmarks/1_preprocess/{seq_category}/4b_intron_expression_table.txt", N_BENCHMARKS)

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

