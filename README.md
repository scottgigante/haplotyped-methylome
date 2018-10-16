## Haplotyped Methylome

Reproducibility instructions for Gigante et al., 2018.

_Note: this repository is under construction!_

Data avilable at [ENA Accession PRJEB27157](https://www.ebi.ac.uk/ena/data/view/PRJEB27157).

[![Directed Acyclic Dependency Graph](dependency_graph.svg)](http://htmlpreview.github.io/?https://github.com/scottgigante/haplotyped-methylome/blob/master/dependency_graph.svg)

### System requirements

* R
* Python>=3.5

#### Dependencies

To install with `conda`, run the following command.

```
conda env create -f environment.yml
```

To install without conda, see the list of dependencies at the bottom of this README.

### Required data

For the standard workflow, `snakemake` will download all the necessary files.

If you wish to avoid running `nanopolish` on the raw read data, you can download these files from our website and store them in the following directory structure.

```
+── README.md
├── genome_data
├── notebooks
├── scripts
└── nanopore
    ├──albacore_1.2.2.b6.sorted.bam
    ├──albacore_1.2.2.b6.phased_sorted.bam
    ├──albacore_1.2.2.b6.methylation.tsv
    ├──albacore_1.2.2.cast.sorted.bam
    ├──albacore_1.2.2.cast.phased_sorted.bam
    ├──albacore_1.2.2.cast.methylation.tsv
    ├──albacore_1.2.2.b6xcast.sorted.bam
    ├──albacore_1.2.2.b6xcast.phased_sorted.bam
    ├──albacore_1.2.2.b6xcast.methylation.tsv
    ├──albacore_2.7.7.castxb6.promethion.sorted.bam
    ├──albacore_2.7.7.castxb6.promethion.phased_sorted.bam
    └──albacore_2.7.7.castxb6.promethion.methylation.tsv
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
* [BWA](https://sourceforge.net/projects/bio-bwa/files/)
* [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* [SNPsplit](https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/)
* [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Nanopolish](https://nanopolish.readthedocs.io/en/latest/installation.html) (optional)
* [Pandoc](https://pandoc.org/installing.html)

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
