# WGS workflow using bwa-mem and GATK

## Alignment

Often times FASTA file are generated using multiple sequencing lanes (e.g., L1_R1.fq with L1_R2.fq and L2_R1.fq with L2_R2.fq). It is important to do the alignment step on each lane separately (making sure specific RG tags to each lane are added), remove or mark the duplicates, and only then merge the resulting BAM files for donwstream processing. This strategy ensures we can correctly track the "read group" (RG tag) information.

### Indexing 

JAX: [`bwa:0.7.17--hed695b0_6.`](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/bwa/bwa_index.nf)

```
bwa index ${fasta}
```

KMT: [bwa version 0.7.18-r1243-dirty cloned from GitHub on 07/05/2024](https://github.com/lh3/bwa)
```
~/tools/bwa/bwa index -p GRCh38.v41 -a bwtsw ${graft_fasta}
```

### Alignment

Adding read group information is an improtant step while aligning the reads. This is later used by GATK pipeline when calling variants.

[JAX read groups](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/utility_modules/read_groups.nf)
[JAX python script to collect read group info from FASTQ files](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/bin/shared/read_group_from_fastq.py)

Our fasta has this info:
```
@LH00516:106:22C7F5LT4:1:1109:17508:16773 1:N:0:GTAAGCTCCA+TGTGCGGTAT
```

JAX:
```
# Example id: HISEQ2000:190:D19U8ACXX:5:TAGCTT
[...]
line = '@RG\\tID:{0}{1}\\tLB:{0}{2}\\tSM:{3}\\tPL:ILLUMINA'.\
```



[JAX:](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/bwa/bwa_mem_hla.nf)
```
container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'

rg=\$(cat $read_groups)
run-bwamem -t $task.cpus -R \${rg} -o ${sampleID}_${index} -H ${params.ref_fa_indices} $inputfq | sh
```

ME:
```
```

## Variant calling

In this workflow, GATK's HaplotypeCaller in GVCF mode for single-sample genotype calling is run. While the GVCF mode is particularly powerful for joint genotyping across multiple samples, it also offers several benefits when applied to single-sample variant calling. in particular, by using GVCF mode for single-sample analyses, one achieves both consistent methodology and the ability to efficiently scale your analysis to include additional samples in the future.

The benefits of running this workflow in GVCF mode are:

 * Consistency for Future Analyses: By generating a GVCF for single samples, one has an option to later combine this single sample with additional samples for joint genotyping, if project objectives change;
 
 * Improved Variant Calling Accuracy: The GVCF format captures more detailed information about the variant calling process, including per-site and per-sample metrics that can improve the accuracy and confidence of calls, even for single samples.
 
 * Scalability: If you start with single-sample analyses but plan to scale up to multi-sample analyses in the future, having GVCF files allows you to easily integrate new samples into a joint genotyping workflow.
 
 * Variant Recalibration: Generating GVCFs can facilitate more accurate variant recalibration and filtering because of the additional information captured in the GVCF.


[JAX GATK workflow](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/tree/main/modules/gatk)
[JAX germline WGS nf](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/workflows/wgs.nf)
[JAX somatic WES nf](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/workflows/somatic_wes.nf)
 
## Usage

To run WGS workflow look follow the following:

```bash
bin/wgs_workflow.sh -t tool_name -g path/to/indexed_genome -r path/to/reads.fastq -o output_directory ...s

```
 
