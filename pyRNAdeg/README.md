# pyRNA


## Overview

`pyRNA` is a python **library** that integrates **modules/scripts** used during our analysis. 
These **modules/scripts** are normally called from the `Notebooks` during some step of the pipeline, but contain useful function definitions and can also act as standalones.

We can classify the scripts in 3 main categories (modules?):
* **Pre-process**:
  * Filter Features: in this case remove ribosomal RNA (rRNA) and Tag bam files with gene IDs
      * `RemrRNA.py` (only for RNA-seq pipeline)
  * Gene Counts: compute raw and TPM-normed Gene Count Tables
      * `GeneExpresssionTableChIP.py` (exactly the same as `gene_expression_table.py`)
      * `gene_expression_table.py`
      * `gene_expression_table_extended.py` (include strand-specific analysis)
  * Replicates: functions related to replicates analysis
      * RepTools.py
* **Coverage**:
  * `ChIP_coverage_files.py`
  * `RNA_coverage_files.py`
  * `coverage.py`
  
* **Visualization**:
  * `Viz.py`
  * `viz_strands.py` (include strand-specific analysis)

Additionally, this directory contains an `env.txt` needed to reproduce the exact original conda enviroment used in our analysis.
(I think that a minimal version of this should be re-created Parastou's enviroment contained lots of stuff totally unnecessary for the analysis).

Faltan:
`Util.py
