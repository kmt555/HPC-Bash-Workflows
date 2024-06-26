# Reference Genome and Indexing

## Choosing the Reference Genome

When selecting a reference genome, ensure it matches the species and version required for your analysis. Common sources for downloading reference genomes include:
- [NCBI](https://www.ncbi.nlm.nih.gov/)
- [Ensembl](https://www.ensembl.org/)
- [UCSC Genome Browser](https://genome.ucsc.edu/)

After downloading, verify the integrity of the genome file using checksums provided by the source.

## Indexing the Reference Genome

Indexing the reference genome is necessary for many bioinformatics tools. This repository includes a script to automate this process for tools such as `bwa` and `bowtie2`.

### Usage

To index a reference genome, run the following command:

```bash
bin/index_reference.sh -g path/to/genome.fasta -t tool_name
```