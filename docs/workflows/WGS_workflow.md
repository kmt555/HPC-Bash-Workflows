# WGS workflow using GATK

In this workflow, GATK's HaplotypeCaller in GVCF mode for single-sample genotype calling is run. While the GVCF mode is particularly powerful for joint genotyping across multiple samples, it also offers several benefits when applied to single-sample variant calling. in particular, by using GVCF mode for single-sample analyses, one achieves both consistent methodology and the ability to efficiently scale your analysis to include additional samples in the future.

The benefits of running this workflow in GVCF mode are:

 * Consistency for Future Analyses: By generating a GVCF for single samples, one has an option to later combine this single sample with additional samples for joint genotyping, if project objectives change;
 
 * Improved Variant Calling Accuracy: The GVCF format captures more detailed information about the variant calling process, including per-site and per-sample metrics that can improve the accuracy and confidence of calls, even for single samples.
 
 * Scalability: If you start with single-sample analyses but plan to scale up to multi-sample analyses in the future, having GVCF files allows you to easily integrate new samples into a joint genotyping workflow.
 
 * Variant Recalibration: Generating GVCFs can facilitate more accurate variant recalibration and filtering because of the additional information captured in the GVCF.
 
### Usage

To run WGS workflow look follow the following:

```bash
bin/wgs_workflow.sh -t tool_name -g path/to/indexed_genome -r path/to/reads.fastq -o output_directory ...s

```
 