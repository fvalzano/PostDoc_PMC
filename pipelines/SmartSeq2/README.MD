This folder contains script to perform SmartSeq2 sequencing analysis of single cell plate sorted.
The overall workflow is:
Fastq -> (Trim_galore) -> FastQC -> MultiQC  -> STAR  -> Sorting -> featureCounts -> counts matrix

The modules and softwares to perform the analysis are deposited in the PMC module directories (CentOS7), MultiQC and Trim galore are deposited in the ss2 environment in the shared space of the group.