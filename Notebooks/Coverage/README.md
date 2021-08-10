
# 2. Coverage

Similarly to the previous step, we will start from the read alignment (`.bam`) files, and:

  * **1.** Create corresponding scaffold coverage (`.bed`) files using `bedtools`.

  * **2.** Use the scaffold coverage (`.bed`) to calculate the actual coverage (`.coverage.txt`) files using some custom functions defined in the notebook.

  * **3.** Normalize the actual coverage (`.coverage.txt`) files generating read-per-million (RPM) normed coverages (`.normed.coverage.txt`) files

**Observations**:

* In the **Coverage** step the main difference between ChIP-seq (DNA) and RNA-seq pipelines is going to be **how we count reads**:

  * In the RNA-seq pipeline, we are going to have to split the read alignment (`.bam`) files into **separate strands** bam (`*.forward.bam`, `.reverse.bam`) files, and the counting will be done in a strand-spefic way.
  * Additionally, **aligned reads** computed using an **splice-aware alignment tool** (e.g. `STAR`) will contain gapped alignments, hence alignment to **to two (or more) disjunct intervals** on the genome might happen (tyically with an intron in between). Therefore, the functions used to compute the coverage (i.e. `generate_coverage_dict`) will be slightly different.
  
* Notebooks seem to be **incomplete**, and under developmment with unfinished ideas (...). It seems to me that the `BamCoveragePlot_ChIP.ipynb` is older than the `BamCoveragePlot-RNA.ipynb`:
  * In both cases, we should create **scaffold coverage** (`.bed`) files using `bedtools` from `bam` files. In the `BamCoveragePlot-RNA.ipynb` these are missing I included the necessary steps to compute them.
  * On the other hand, in the next steps similar **functions** (with same names) are defined in each notebook but are slightly different, and **Step 2: calculate coverage (`.txt`)** and **Step 3: calculate normed coverage (`.normed.txt`)** are grouped in the same `for` loop in `BamCoveragePlot-RNA.ipynb`.
  
  * `bed_coverage_dict()`: takes a scaffold `bed` file as input to create a `bed_cov_dict` **object**.
    * In `BamCoveragePlot_ChIP.ipynb` the `bed_cov_dict` **object** is initialized using both the `key`'s and `value`'s from the `bed` file. 
    * Instead, in `BamCoveragePlot-RNA.ipynb` the `bed_cov_dict` **object** is initialized using only the `key`'s and all `value`'s are set to `0`.

  * `generate_coverage_dict()`: takes a `bam` file as input 
    * In `BamCoveragePlot_ChIP.ipynb`:
      * initializes an **empty dictionary**: `coverage_dict = {}`
      * *"manually"*  computes a `coverage_dict` by using the `update_coverage_dict()` **function** that checks if a `key` already exists otherwise **generates** it.
      
     * Instead, in `BamCoveragePlot-RNA.ipynb` 
       * the `generate_coverage_dict()` **function**, now take takes a `bam` **and**  a scaffold `bed` file.
       * initializes a **non empty dictionary** using the: `coverage_dict = bed_coverage_dict(in_bed)`. As metioned before, this dictionary is non empty although all `value`'s are zero.
       * *"manually"*  computes a `coverage_dict` by using the `update_coverage_dict()` **function** that now only **updates** the `key`'s that already exists because were obtained from the scaffold `bed` file.
  
  * Finally, the resulting is `coverage_dict` is writen to one coverage (`.txt`) output file.
    * In `BamCoveragePlot_ChIP.ipynb`, before outputting the result:
      * create a `bed_cov_dict` **object** from the scaffold `bed` file. As mentioned before, in this case the `bed_cov_dict` **object** is initialized using both the `key`'s **and** `value`'s from the `bed` file.
      * both **objects** `bed_cov_dict` and `coverage_dict` are merged into one coverage (`.txt`) output file.
    * Instead, in `BamCoveragePlot-RNA.ipynb` there is no merge, just write otput into one coverage (`.txt`) file.

       
Additionally, there are corresponding python scripts (which I assume are even newer) in `pyRNAdeg`:
  * ChIP_coverage_files.py
  * RNA_coverage_files.py
  
  these assume that you have the scaffold coverage (`.bed`) files computed using `bedtools`, but are pretty much the same that the notebooks. Finally, there is another python script:
  * `coverage.py`
 
 that seems to be more general and already include both analysis, but seems to be **totally incomplete**.


## A. ChIP-seq: `BamCoveragePlot_ChIP.ipynb`

**1.**  Create the scaffold coverage (`.bed`) files using [`bedtools`](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) from `bam`
```
!$bedtools genomecov -d -ibam $in_bam > $out_path
```
![alt text](https://bedtools.readthedocs.io/en/latest/_images/genomecov-glyph.png)

`bedtools genomecov` computes histograms (default), **per-base reports (-d)** and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.

* Option	Description:
  * `-ibam`: BAM file as input for coverage. Each BAM alignment in A added to the total coverage for the genome.
  * `-d`:	Report the depth at each genome position with 1-based coordinates.
  * `-split`:	Treat “split” BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR “N” and “D” operations to infer the blocks for computing coverage.
  * `-strand`:	Calculate coverage of intervals from a specific strand. With BED files, requires at least 6 columns (strand is column 6).
  * `-scale`: Scale the coverage by a constant factor. Each coverage value is multiplied by this factor before being reported. Useful for normalizing coverage by, e.g., reads per million (RPM). Default is 1.0; i.e., unscaled.

**2.** Calculate coverage (`.txt`) files from `bam` and scaffold coverage (`.bed`) files

Main steps are: (I'm still not sure what this is doing)
  * load bam file 
  ```
  st = pysam.AlignmentFile(open(in_bam, 'rb'))
  ```
  
  * create coverage information for mapped regions
  ```
  coverage_dict = generate_coverage_dict(st)
  ```
  
  * load bedtools coverage file as scafold, create coverage dictionary:
  ```
  bed_cov_dict = bed_coverage_dict(in_bed)
  assert len(bed_cov_dict) >= len(coverage_dict)
  ```
  
  * merge coverage information into a new file
  ```
  with open(out_path, 'w+') as out:

   for coord, item in bed_cov_dict.items():
   
     coverage = item
     if coord in coverage_dict:
         coverage = coverage_dict[coord]

     ch = coord.split(':')[0]
     pos = coord.split(':')[1]
     line = ch + '\t' + str(pos) + '\t' + str(coverage) + '\n'

     out.write(line)
  ```

**3.**  Calculate read-per-million (RPM) normed coverage (`.normed.txt`) files

This is a trivial step only involves:
  * Load in: coverage file (.txt)
  ```
  df = pd.read_csv(coverage_file, sep='\t', names=range(3))
  ```
  
  * divide by read length: 50 bases
  ```
  total = sum(df[2]) / 50
  ```
  
  * divide this number by 1,000,000. ('per million' scaling factor)
  ```
  scaling_factor = (total / 1000000)
  ```
  
## B. RNA-seq: `SeparateStrands-RNA.ipynb` + `BamCoveragePlot-RNA.ipynb`

First of all one needs to execute the: `SeparateStrands-RNA.ipynb` Notebook
* Should this be moved to the pre-processing? Maybe after removing rRNA and mentioning that will be needed for the coverage analysis, why is it not needed for gene counts?
* Alternative to: `SeparateStrands-RNA.ipynb`
  * [Splitting reads in a BAM file by strands](https://josephcckuo.wordpress.com/2016/11/18/splitting-reads-in-a-bam-file-by-strands/)
  * Using [htseq](https://htseq.readthedocs.io/en/release_0.11.1/tour.html#genomic-intervals-and-genomic-arrays) library to compute coverage



