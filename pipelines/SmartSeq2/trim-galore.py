import os
import subprocess

# Define the base directory containing FASTQ files in subdirectories
base_dir = "/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/fastq"

# Define the output directory for trimmed FASTQ files
output_dir = "/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/trimmed_fastq"
os.makedirs(output_dir, exist_ok=True)

# Walk through the base directory and find all FASTQ files
fastq_pairs = {}
for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith(".fastq.gz"):
            # Identify if the file is R1 or R2
            if "_R1.fastq.gz" in file:
                sample_name = file.replace("_R1.fastq.gz", "")
                fastq_pairs.setdefault(sample_name, {})['R1'] = os.path.join(root, file)
            elif "_R2.fastq.gz" in file:
                sample_name = file.replace("_R2.fastq.gz", "")
                fastq_pairs.setdefault(sample_name, {})['R2'] = os.path.join(root, file)

# Run Trim Galore on paired-end FASTQ files
for sample, fastq_files in fastq_pairs.items():
    if 'R1' in fastq_files and 'R2' in fastq_files:
        r1_fastq = fastq_files['R1']
        r2_fastq = fastq_files['R2']
        print(f"Running Trim Galore on paired files: {r1_fastq} and {r2_fastq}...")

        # Run Trim Galore for paired-end reads
        subprocess.run([
                        "trim_galore", 
                        "--paired", 
                        r1_fastq, 
                        r2_fastq, 
                        "-o", 
                        output_dir
        ])

print("Trim Galore analysis completed for all paired-end FASTQ files.")
