#!/bin/bash
snakemake --dag all | dot -Tsvg > dependency_graph.svg
snakemake --dag haplotype_analysis | dot -Tsvg > haplotype_dependency_graph.svg
snakemake --dag methylation_analysis | dot -Tsvg > methylation_dependency_graph.svg
snakemake --dag rnaseq_analysis | dot -Tsvg > rnaseq_dependency_graph.svg
