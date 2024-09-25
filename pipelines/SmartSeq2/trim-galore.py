#before executing this script, activate the ss2 environment cntaining multiqc -> conda activate ss2
import os
import subprocess

# Define the base directory containing FASTQ files in subdirectories
base_dir = "/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/"

# Define the output directory for trimmed FASTQ files
output_dir = "/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/trimmed_fastq"
os.makedirs(output_dir, exist_ok=True)

# Walk through the base directory and find all FASTQ files
for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith(".fastq.gz"):
            fastq_file = os.path.join(root, file)
            print(f"Running Trim Galore on {fastq_file}...")
            
            # Run Trim Galore
            subprocess.run([
                            "trim_galore", 
                            fastq_file, 
                            "-o", 
                            output_dir
            ])

print("Trim Galore analysis completed for all FASTQ files.")
