Workflow for smart seq2 analysis:

Fastq -> FastQC -> MultiQC  -> STAR  -> featureCounts -> counts matrix
      -> Trim_galore -> FastQC_post_trimming -> MultiQC_post_trimming  -> STAR  -> featureCounts -> counts matrix      

FastQC -> /hpc/local/CentOS7/pmc_research/etc/modulefiles
STAR -> /hpc/local/CentOS7/pmc_research/etc/modulefiles
subread -> /hpc/local/CentOS7/pmc_research/etc/modulefiles
MultiQC -> /hpc/local/Rocky8/pmc_kool/miniconda3/envs/ss2
Trim-galore -> /hpc/local/Rocky8/pmc_kool/miniconda3/envs/ss2