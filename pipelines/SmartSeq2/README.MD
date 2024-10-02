Workflow for smart seq2 analysis:

Fastq -> (Trim_galore) -> FastQC -> MultiQC  -> STAR  -> Sorting -> featureCounts -> counts matrix

FastQC -> /hpc/local/CentOS7/pmc_research/etc/modulefiles
STAR -> /hpc/local/CentOS7/pmc_research/etc/modulefiles
subread -> /hpc/local/CentOS7/pmc_research/etc/modulefiles
MultiQC -> /hpc/local/Rocky8/pmc_kool/miniconda3/envs/ss2
Trim-galore -> /hpc/local/Rocky8/pmc_kool/miniconda3/envs/ss2