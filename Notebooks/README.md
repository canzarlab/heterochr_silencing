# Notebooks

This directory contains the actual computational pipeline as a set of jupyter Notebooks that need to be executed individually in a particular order.

## Overview

Remember that before running the pipeline one should activate the proper [conda enviroment](https://github.com/pablommesas/RNAdeg/blob/master/README.md#conda-enviroment).

The main structure of the pipeline is the following:

1.  **PreProcess**
2.  **Coverage** (Partially)
3.  **RatioAnalyses**
4.  **Plots**

there is an additional Notebook used to organize the samples/replicates under investigation:

* `sample_names.ipynb`

and two additional directories extending the analysis:

* **S5**: I think this is the particular analysis of mutant samples **S5**
* **Stranded**: strand specific analysis??

As previously stated, in our paper we analyzed **next-generation sequencing (NGS) data** for:

* RNA Polymerase II (**ChIP-seq**)
* Pol II bound nascent RNA (**RIP-seq**)
* steady state RNA levels (**pA RNA-seq**)

in S. pombe fission yeast.

Each technology is describing a different biological process, but from the analysis perspective the main difference will be to distinguish between the data's origin:
* **DNA**:
  * ChIP-seq
* **RNA**
  * RIP-seq
  * pA RNA-seq

This means that on the **PreProcessing** and **Coverage** analysis we will have two distinct pipelines, one for the **ChIP-seq** data and the other for **RIP-seq** and **pA RNA-seq**.
