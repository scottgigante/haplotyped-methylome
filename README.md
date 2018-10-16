## Haplotyped Methylome

Reproducibility instructions for Gigante et al., 2018.

_Note: this repository is under construction!_

Data avilable at [ENA Accession PRJEB27157](https://www.ebi.ac.uk/ena/data/view/PRJEB27157).

[![Directed Acyclic Dependency Graph](dependency_graph.svg)](http://htmlpreview.github.io/?https://github.com/scottgigante/haplotyped-methylome/blob/master/dependency_graph.svg)

### System requirements

* R
* Python>=3.5
* [Albacore](https://community.nanoporetech.com/downloads) (optional)
* Lots of RAM (tested on 500GB)
* Lots of disk space

#### Dependencies

To install with `conda`, run the following command.

```
conda env create -f environment.yml
source active haplotyped_methylome
```

To install without conda, see the list of dependencies at the bottom of this README.

### Required data

For the standard workflow, `snakemake` will download all the necessary files.

If you wish to avoid running `albacore`, `bwa` and `nanopolish` on the raw nanopore data, you can run the following command, which downloads the output of these programs and tricks `snakemake` into thinking you have run the pipeline from the beginning:

```
snakemake intermediate_download
```

### Running the workflow

To generate all plots, tables and notebooks, simply run from the root directory:

```
snakemake --cores 16
```

If you don't wish to run the full analysis, you can run specific rules from the Snakefile by running, for example:

```
snakemake --cores 16 rnaseq_analysis
```

### Installation without `conda`

Software dependencies:

* [SAMtools](http://www.htslib.org/download/)
* [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* [SNPsplit](https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/)
* [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Pandoc](https://pandoc.org/installing.html)
* [BWA](https://sourceforge.net/projects/bio-bwa/files/) (optional)
* [Nanopolish](https://nanopolish.readthedocs.io/en/latest/installation.html) (optional)

Python package dependencies:

```
pip install --user -r requirements.txt
```

R package dependencies:

```
Rscript install_R_deps.R
```

[![Directed Acyclic Dependency Graph: Methylation](methylation_dependency_graph.svg)](http://htmlpreview.github.io/?https://github.com/scottgigante/haplotyped-methylome/blob/master/methylation_dependency_graph.svg)

[![Directed Acyclic Dependency Graph: Haplotyping](haplotype_dependency_graph.svg)](http://htmlpreview.github.io/?https://github.com/scottgigante/haplotyped-methylome/blob/master/haplotype_dependency_graph.svg)

[![Directed Acyclic Dependency Graph: RNA-seq](rnaseq_dependency_graph.svg)](http://htmlpreview.github.io/?https://github.com/scottgigante/haplotyped-methylome/blob/master/rnaseq_dependency_graph.svg)
