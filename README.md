# HPC-Bash-Workflows

HPC-Bash-Workflows is a collection of bash and shell scripts designed to streamline job submissions and manage workflows on High-Performance Computing (HPC) systems. These scripts provide a robust, efficient, and customizable solution for automating repetitive tasks, handling job dependencies, and optimizing resource utilization on HPC clusters. Whether you're a researcher, bioinformatician, or data scientist, this repository aims to simplify your HPC job management, enhancing productivity and reproducibility.

## Choice of Reference Genome and Indexing

For many genomic workflows, choosing the correct reference genome and properly indexing it is crucial for accurate and efficient analysis. This repository includes scripts and configurations to assist with these tasks.

### Reference Genome

Ensure that you have the correct reference genome for your analysis. Reference genomes can be downloaded from various sources such as NCBI, Ensembl, or UCSC.

### Indexing the Reference Genome

Indexing the reference genome is a necessary step before running many bioinformatics tools. This repository provides a script to automate the indexing process.

To index a reference genome, use the following script:

```bash
bin/index_reference.sh -g path/to/genome.fasta -t tool_name
```