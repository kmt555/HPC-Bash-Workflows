- [WGS workflow using bwa-mem and GATK](#wgs-workflow-using-bwa-mem-and-gatk)
  * [Genome Indexing](#genome-indexing)
  * [Alignment](#alignment)
    + [Read Groups '@RG'](#read-groups---rg-)
    + [Map and mark duplicates](#map-and-mark-duplicates)
    + [Merge BAMs from multiple lanes](#merge-bams-from-multiple-lanes)
  * [Variant calling](#variant-calling)
  * [Usage](#usage)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


# WGS workflow using bwa-mem and GATK

FASTA files are often generated using multiple sequencing lanes (e.g., L1_R1.fq with L1_R2.fq and L2_R1.fq with L2_R2.fq). It is crucial to perform the alignment step on each lane separately, ensuring that specific read group (RG) tags are added to each lane. After alignment, duplicates should be removed or marked before merging the resulting BAM files for downstream processing. This strategy ensures accurate tracking of read group information, which is vital for subsequent analyses.

NOTE: Theoretically, it is possible to use Gencode GRCh38.primary_assembly with GATK hg38 bundle resources. However, one should ensure compatibility by verifying that the reference sequences, naming conventions, and coordinates match across all resources used in the analysis.

**Example:** make sure that chromosome lengths are the same. Use `samtools faidx` to create fai files.

```
$ cut -f 1,2 /PATHTO/human_GRCh38.v41/GRCh38.primary_assembly.genome.fa.fai | head
chr1	248956422
chr2	242193529
chr3	198295559
chr4	190214555
chr5	181538259
chr6	170805979
chr7	159345973
chr8	145138636
chr9	138394717
chr10	133797422
...

$ cut -f 1,2 /PATHTO/GATK_Hg38_bundle/Homo_sapiens_assembly38.fasta.fai | head
chr1	248956422
chr2	242193529
chr3	198295559
chr4	190214555
chr5	181538259
chr6	170805979
chr7	159345973
chr8	145138636
chr9	138394717
chr10	133797422
...
```

A more comprehensive analysis is also recommended, such as sequence comparison or using specialzied tools like [MUMmer](https://mummer.sourceforge.net)

```
diff <(grep -v '^>' GRCh38.primary_assembly.genome.fa) <(grep -v '^>' Homo_sapiens_assembly38.fasta)
```

## Genome indexing for read mapping

KMT: [bwa version 0.7.18-r1243-dirty cloned from GitHub on 07/05/2024](https://github.com/lh3/bwa)
```
~/tools/bwa/bwa index -p GRCh38.v41 -a bwtsw ${graft_fasta}
```

JAX: [`bwa:0.7.17--hed695b0_6.`](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/bwa/bwa_index.nf)
```
bwa index ${fasta}
```

## Alignment

### Read Groups '@RG'

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

### Map and mark duplicates

The most critical aspect of setting up the alignment run is properly configuring the read group information. For example, given the following FASTA header `@LH00516:106:22C7F5LT4:1:1109:17508:16773 1:N:0:GTAAGCTCCA+TGTGCGGTAT`, we can set the RG tag as follows:

```
header=$(gunzip -c ${FASTA_R1} | head -n 1)

ID=$(echo $header | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g' )
LB=$(echo $header | cut -f 4 -d":")
SM=$(echo ${sampleID} | cut -f 1-4 -d"_")
PL=ILLUMINA

~/tools/bwa/bwa mem -K 100000000 -Y -t ${cpus} -R '@RG\tID:'$ID'\tLB:'$LB'\tSM:'$SM'\tPL:'$PL ${ref_indx} ${FASTA_R1} ${FASTA_R2} | ~/tools/samblaster/samblaster -a --addMateTags | samtools view -b -S - > ./BAMs/${sampleID}.bam
```

This initial alignment step is followed by `java -Xmx16g -jar $picardtools SortSam` and `MarkDuplicates`, and `samtools sort` and `index.`

```
# sort by query name:
java -Xmx16g -jar $picardtools SortSam I=./BAMs/${sampleID}.bam O=./BAMs/${sampleID}.qname.bam SORT_ORDER=queryname

# mark duplicates:
java -Xmx16g -jar $picardtools MarkDuplicates I=./BAMs/${sampleID}.qname.bam O=./BAMs/${sampleID}.mdups.bam M=./METRICS/${sampleID}.mdups_metrics.txt ASSUME_SORT_ORDER=queryname

# and finally sort and index with samtools
samtools sort -m 7680MiB -o ./BAMs/${sampleID}.mdups.sorted.bam ./BAMs/${sampleID}.mdups.bam
samtools index ./BAMs/${sampleID}.mdups.sorted.bam ./BAMs/${sampleID}.mdups.sorted.bam.bai
```

[JAX:](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/bwa/bwa_mem_hla.nf) has developed a Python script for read group assignment, which is a valuable resource to consider if the above method does not fully meet specific project needs:

- [JAX read groups](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/modules/utility_modules/read_groups.nf)

- [JAX python script to collect read group info from FASTQ files](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/blob/main/bin/shared/read_group_from_fastq.py)

```
container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'

rg=\$(cat $read_groups)
run-bwamem -t $task.cpus -R \${rg} -o ${sampleID}_${index} -H ${params.ref_fa_indices} $inputfq | sh
```

### Merge BAMs from multiple lanes

As specified above, we merge BAM files from multiple lanes to streamline downstream processing and ensure accurate variant calling.

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
 
