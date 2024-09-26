#!/bin/bash
#SBATCH --time=24:00:00 # Walltime
#SBATCH --mem=64G
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=8 # number of threads we want to run on
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ref_genome_star/log.out
#SBATCH --job-name=STAR_genome_generation

#Do not send the job on the head node but reserve a node with salloc...
module use --append /hpc/local/CentOS7/pmc_research/etc/modulefiles
module load STAR
GENOME_SIZE=3000000000   # Replace with approximate genome size in base pairs (e.g., 3,000,000,000 for human)
THREADS=8  
STAR --runMode genomeGenerate  \
     --runThreadN 8 \
     --genomeDir /hpc/pmc_kool/fvalzano/Ref_genome_star/star \
     --genomeFastaFiles /hpc/pmc_kool/fvalzano/Ref_genome_star/fasta/genome.fa \
     --sjdbGTFfile /hpc/pmc_kool/fvalzano/Ref_genome_star/genes/genes.gtf \
     --sjdbOverhang 100