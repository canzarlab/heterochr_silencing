
# 1. PreProcessing

In this initial step, we start with the **raw read (`.fastq`)** files, align the reads to a reference genome using `STAR`, count the number of mapped reads (`.bam`) to each gene and summarize everything in raw and TPM-normed **Gene Count/Expression Tables** using `gene_expression_table.py`.

**Observations**:

* In the Pre-processing step the main differences between ChIP-seq (DNA) and RNA-seq pipelines are:
  * The `STAR` alignment configuration, in the ChIP-seq (DNA) alignments we won't allow for **splicing**.
  * In the RNA-seq pipeline there is an additional step where we remove **ribosomal RNA (rRNA)** and tag bam files with gene IDs (using `RemrRNA.py`).
  
* How to count reads `gene_expression_table.py`? For each read check the `NH` (**number of hits**) column in `.bam` file and which **region of the genome** it mapps to:
  * **Euchromatic**/Protein Coding (PC) Genes: `NH = 1` (Only count unique genes)
  * **Heterochromatic** (HTC) Genes: `NH < 15` (select and count 1 of the genes)

* ChIP-seq File types: each IP file must have a corresponding Input file (Reviewers)
  * ChIP
  * Input
  * IP

## A. ChIP-seq:  `PreProcess_ChIP.ipynb`

**1.**  Prepare output folders
```
if not os.path.isdir(dir):    
    !mkdir -p $dir
```
**2.**  Unzip fastq files
```
 !bzip2 -kdv $raw_file
```
  * (Optional step): remove truncated, 0-length, reads from fastq files.
  ```
  !$bioawk -cfastx 'length($seq) > 0 {print "@"$name"\n"$seq"\n+\n"$qual} $filepath  >> $outfilepath'
  ```
 
**3.**  Align fastq files: `STAR` configuration different from RNA-seq, because reads come from DNA. Don't allow for splicing (`--alignIntronMax 1 --alignEndsType EndToEnd`)

  * Run `STAR` aligner
  ```
  !star + ' --runThreadN ' + str(n_threads) + ' --genomeDir ' + genome_dir + ' --readFilesIn ' + fastq_file + ' --outFileNamePrefix ' + bam_prefix + ' --outSAMtype BAM SortedByCoordinate --alignIntronMax 1 --alignEndsType EndToEnd'
  ```
  
  * Index alignment (bam) files with `samtools`
  ```
  !samtools index $filepath
  ```
  
**4.**  Compute raw and TPM-normed Gene Count Tables: `gene_expression_table.py`

```
python + ' ' + 'gene_expression_table.py' + ' -d ' + out_bam + ' -g ' + annotation + ' -o ' + xp_data + ' -x ' + 'chip_'
```

## B. RNA-seq: `PreProcess_RNA.ipynb`

**1.**  Prepare output folders
```
if not os.path.isdir(dir):    
    !mkdir -p $dir
```

**2.**  Unzip fastq files
```
 !bzip2 -kdv $raw_file
```
  * (Optional step): remove truncated, 0-length, reads from fastq files.
  ```
  !$bioawk -cfastx 'length($seq) > 0 {print "@"$name"\n"$seq"\n+\n"$qual} $filepath  >> $outfilepath'
  ```
 
**3.**  Align fastq files: `STAR` regular RNA-seq configuration

  * Run `STAR` aligner
  ```
  !star + ' --runThreadN ' + str(n_threads) + ' --genomeDir ' + genome_dir + ' --readFilesIn ' + fastq_file  + ' --outFileNamePrefix ' + bam_prefix + ' --outSAMtype BAM SortedByCoordinate'
  ```
  
  * Index alignment (bam) files with `samtools`
  ```
  !samtools index $filepath
  ```
  
**4.** Remove ribosomal RNA (rRNA) and Tag bam files with gene IDs: `RemrRNA.py`
```
python + ' ' + 'RemrRNA.py' + ' -d ' + out_bam + ' -g ' + annotation + ' -o ' + out_tagged + ' -i ' + ' '.join(ignore_files)
```

**5.**  Compute raw and TPM-normed Gene Count Tables: `gene_expression_table.py`
```
python + ' ' + 'gene_expression_table.py' + ' -d ' + out_tagged + ' -g ' + annotation + ' -o ' + xp_data + ' -x ' + 'rna_'
```
