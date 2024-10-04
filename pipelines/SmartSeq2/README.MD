# Overview of the folder

    script_dir=${hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2}

This folder contains script to perform SmartSeq2 sequencing analysis of single cell plate sorted.


The overall workflow is:
Run FastQC for fastq files:
    
    bash ${script_dir}/FastQC.sh

Run MultiQC to generate single QC report:

    bash ${script_dir}/MultiQC.sh

If needed, trim adapter from sequencing read with Trim Galore:

    bash ${script_dir}/trim-galore.sh

Run FastQC for trimmed fastq files:

    bash ${script_dir}/FastQC_post_trimming.sh

Run MultiQC to generate single QC report:

    bash ${script_dir}/MultiQC_post_trimming.sh

Run Star alignment:

    bash ${script_dir}/STAR.sh

Sort the output .bam files by name:

    bash ${script_dir}/Sorting_bam.sh

Perform quantification with subread(featureCounts):

    bash ${script_dir}/featureCounts.sh


The modules and softwares to perform the analysis are deposited in the PMC module directories (CentOS7), MultiQC and Trim galore are deposited in the ss2 environment in the shared space of the group.