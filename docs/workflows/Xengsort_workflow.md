# Xengsort workflow PDX data

```
source ~/tools/enable_conda.sh
conda activate xengsort
```


```
host_fasta=/pathto/mouse_mm39.vM30/GRCm39.primary_assembly.genome.fa
graft_fasta=/pathto/human_GRCh38.v41/GRCh38.primary_assembly.genome.fa
xengsort_idx_name=GRCm39vM30_GRCh38v41

## indexing human and mouse reference genomes
xengsort index \
    --index ${xengsort_idx_name} \
    -H ${host_fasta} \
    -G ${graft_fasta} \
    --fill 0.88 \
    -k 25 \
    -n 4_500_000_000 \
    -W ${cpus}
```


```
xengsort classify \
        --index ${xengsort_index}/${xengsort_idx_name} \
        --fastq ${FASTA_R1} \
        --pairs ${FASTA_R2} \
        --prefix ${sampleID} \
        --mode count \
        --threads ${cpus} \
        --out ${sampleID} \
        --chunksize 32.0 \
        --compression gz &> ${sampleID}_xengsort_log.txt
```


```
conda deactivate
```


