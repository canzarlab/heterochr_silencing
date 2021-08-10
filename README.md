# HetchromSilencing

This repository contains the computational pipeline of the paper **Ccr4-Not complex reduces transcription efficiency in heterochromatin**. It contains the scripts and software necessary for reproducing the results in the paper.

## Overview

Current data suggest that heterochromatic silencing is a combination of transcriptional silencing and RNA degradation, however, the contribution of each pathway to silencing is not known. In this study we analyzed RNA Polymerase II (Pol II) occupancy, nascent RNA levels and steady state RNA levels to quantify the contribution of these pathways to heterochromatic silencing in fission yeast.

We found that transcriptional silencing consists of two components, reduced RNA Pol II accessibility and reduced transcriptional efficiency. Using our data we determined and quantified which heterochromatic complexes contribute to reduced RNA Pol II occupancy, reduced transcriptional efficiency and reduced RNA stability.

In our paper we analyzed next-generation sequencing data for:

* RNA Polymerase II (**ChIP-seq**)
* Pol II bound nascent RNA (**RIP-seq**)
* steady state RNA levels (**pA RNA-seq**)

in S. pombe fission yeast.

The [notebooks](https://github.com/canzarlab/heterochr_silencing/tree/master/Notebooks) folders contain Jupyter notebooks that make up the actual pipe-line: [preprocessing](https://github.com/canzarlab/heterochr_silencing/tree/master/Notebooks/PreProcess), [analysing](https://github.com/canzarlab/heterochr_silencing/tree/master/Notebooks/RatioAnalyses), [visualizing results](https://github.com/canzarlab/heterochr_silencing/tree/master/Notebooks/Plots) and [computing coverage](https://github.com/canzarlab/heterochr_silencing/tree/master/Notebooks/Coverage) files.

All necessary auxiliary functions and scripts called within the pipe-line can be found in the [pyRNAdeg](https://github.com/canzarlab/heterochr_silencing/tree/master/pyRNAdeg) library.

## Software requirements

The following programs are required to run the scripts in this repository.

* [STAR](https://github.com/alexdobin/STAR)
* [Samtools](http://www.htslib.org/download/)
* [bioawk](https://anaconda.org/bioconda/bioawk)
* [bedtools](https://anaconda.org/bioconda/bedtools)

To run the scripts in the Jupyter notebooks, the following Python modules are required.

* [numpy](http://docs.scipy.org/doc/numpy-1.10.1/user/install.html)
* [pandas](https://pandas.pydata.org/)
* [scikit-learn](http://scikit-learn.org/stable/install.html)
* [matplotlib](http://matplotlib.org/users/installing.html)
* [pysam](https://github.com/pysam-developers/pysam)

## Conda Enviroment

In order to make our pipe-line more reproducible, one should always run the Notebooks within the same **conda enviroment** using:

`source activate heterochromatin`

the enviroment used in our analysis can be reproduced with the `py36.yml` file.

In Jupyterhub one will need to create a Jupyter Notebook kernel to launch this new enviromment, [see](https://ipython.readthedocs.io/en/stable/install/kernel_install.html) for more details. After sourcing the desired enviroment:

```
source activate heterochromatin
python -m ipykernel install --user --name heterochromatin --display-name "Python (heterochromatin)"
```

Go back to the Jupyterhub dashboard, reload the page, now you should have another option in the `New` menu that says `heterochromatin`. In order to use your new kernel with an existing notebook, click on the notebook file in the dashboard, it will launch with the default kernel, then you can change kernel from the top menu `Kernel` > `Change kernel`.

In order to up-date/remove a Jupyterhub kernel. List the paths of all your kernels:
```
jupyter kernelspec list
````
Then simply uninstall your unwanted-kernel
```
jupyter kernelspec uninstall unwanted-kernel
````

## Pipeline

**1.** Prepare the environment
  * Root of directory structure
  * Set the environment variables
  
**2.** Prepare the necessary files needed for further processing:
  * Genome tables
  * STAR index files
  
**3.** Align files and rename and move them, use /tmp folder first.

**4.** Index the files in the new folder.

**5.** Calculate gene expression data.
  * Calculate Pearson-r correlation scores.
  * Merge the replicates
  
**6.** Calculate RS, TE, and PO change tables.

**7.** Mass produce all figures, possibly show a compressed version in the notebook.

**8.** Produce coverage files and save them.

**9.** Run the Rscript and save the figures.

