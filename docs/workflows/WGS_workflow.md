# WGS workflow using bwa-mem and GATK

## Alignment

Often times FASTA file are generated using multiple sequencing lanes (e.g., L1_R1.fq with L1_R2.fq and L2_R1.fq with L2_R2.fq). It is important to do the alignment step on each lane separately (making sure specific RG tags to each lane are added), remove or mark the duplicates, and only then merge the resulting BAM files for donwstream processing. This strategy ensures we can correctly track the "read group" (RG tag) information.

### Indexing 

KMT: [bwa version 0.7.18-r1243-dirty cloned from GitHub on 07/05/2024](https://github.com/lh3/bwa)
```
~/tools/bwa/bwa index -p GRCh38.v41 -a bwtsw ${graft_fasta}
```

JAX: [`bwa:0.7.17--hed695b0_6.`](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/bwa/bwa_index.nf)
```
bwa index ${fasta}
```

### Alignment

#### Read Groups '@RG'

When aligning NGS data, setting read group information is crucial, especially for the GATK pipeline used in variant calling. Read groups, indicated by '@RG', are not output by mappers and must be specified by the user during the alignment step. A typical read group format looks like this:

```
'@RG\tID:'$ID'\tLB:'$LB'\tSM:'$SM'\tPL:'$PL
```

In this format, the variables store information about the sample and the sequencing run. The tags ID, LB, SM, and PL are separated by colons.

- ID: This typically stores the flowcell name and lane number, and/or library information, ensuring global uniqueness across all sequencing runs.
- LB: This is a library identifier. If DNA samples were sequenced across multiple lanes, this will be encoded by LB. This information is utilized by PicardTools' MarkDuplicates to mark potential molecular duplicates.
- SM: This represents the sample name and will appear in the final VCF column. GATK processes all read groups with the same SM value as reads from the same biological sample. When sequencing the same sample multiple times, use a general name for the sample rather than individual replicate names.
- PL: per GATK specification, vaild values are: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.

While read groups can be defined according to personal preference when working privately, it is essential to adopt a clear and consistent schema when collaborating with others. This ensures clarity and consistency in the data processing and analysis pipeline.

**Related references:**

- https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups

- https://www.biostars.org/p/280837/

- [JAX read groups](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/utility_modules/read_groups.nf)

- [JAX python script to collect read group info from FASTQ files](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/bin/shared/read_group_from_fastq.py)

#### BWAMEM

KMT:
```
header=$(gunzip -c ${FASTA_R1} | head -n 1)

ID=$(echo $header | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g' )
LB=$(echo $header | cut -f 4 -d":")
SM=$(echo ${sampleID} | cut -f 1-4 -d"_")
PL=ILLUMINA

~/tools/bwa/bwa mem -K 100000000 -Y -t ${cpus} -R '@RG\tID:'$ID'\tLB:'$LB'\tSM:'$SM'\tPL:'$PL ${ref_indx} ${FASTA_R1} ${FASTA_R2} | ~/tools/samblaster/samblaster -a --addMateTags | samtools view -b -S - > ./BAMs/${sampleID}.bam
```

[JAX:](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/bwa/bwa_mem_hla.nf)
```
container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'

rg=\$(cat $read_groups)
run-bwamem -t $task.cpus -R \${rg} -o ${sampleID}_${index} -H ${params.ref_fa_indices} $inputfq | sh
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
 
