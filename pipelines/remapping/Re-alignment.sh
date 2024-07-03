module use --append /hpc/local/Rocky8/pmc_research/etc/modulefiles
#load cellranger version used for previous runs
module load cellranger/6.1.1
cellranger count --localcores=6 --localmem=32 --localvmem=96 --include-introns --transcriptome=/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/remapping/refdata-gex-GRCh38-2020-A --fastqs=/hpc/pmc_kool/fvalzano/PostDoc_PMC/remapping/data_to_remap/YourFastqFiles --sample=YourSample --id YourID