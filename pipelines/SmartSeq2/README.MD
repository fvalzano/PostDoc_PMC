# Overview of the folder

    SCRIPT_DIR=${YourFavouriteDirectory}
This folder contains script to perform SmartSeq2 sequencing analysis of single cell plate sorted.


The overall workflow is:
Run FastQC for fastq files:
    
    bash ${SCRIPT_DIR}/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/FastQC.sh

Run MultiQC to generate single QC report:

    bash ${SCRIPT_DIR}/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/MultiQC.sh

If needed, trim adapter from sequencing read with Trim Galore:

    bash ${SCRIPT_DIR}/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/trim-galore.sh

Run FastQC for trimmed fastq files:

    bash ${SCRIPT_DIR}/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/FastQC_post_trimming.sh

Run MultiQC to generate single QC report:

    bash ${SCRIPT_DIR}/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/MultiQC_post_trimming.sh

Run Star alignment:

    bash ${SCRIPT_DIR}/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/STAR.sh

Sort the output .bam files by name:

    bash ${SCRIPT_DIR}/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/Sorting_bam.sh

Perform quantification with subread(featureCounts):

    bash ${SCRIPT_DIR}/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/featureCounts.sh


The modules and softwares to perform the analysis are deposited in the PMC module directories (CentOS7), MultiQC and Trim galore are deposited in the ss2 environment in the shared space of the group.