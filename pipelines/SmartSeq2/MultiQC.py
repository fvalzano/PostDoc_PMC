#before executing this script, activate the ss2 environment containing multiqc -> conda activate ss2
import subprocess
import os

# Define the directory where FastQC reports are stored
fastqc_output_dir = "/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/fastqc_reports"
multiqc_output_dir = "/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/multiqc_reports"

# Create the output directory for MultiQC report
if not os.path.exists(multiqc_output_dir):
    os.makedirs(multiqc_output_dir)

# Command to run MultiQC
cmd = f"multiqc {fastqc_output_dir} -o {multiqc_output_dir}"

try:
    print(f"Running MultiQC on the FastQC output in '{fastqc_output_dir}'...")
    subprocess.run(cmd, shell=True, check=True)
    print(f"MultiQC report generated successfully in '{multiqc_output_dir}'")
except subprocess.CalledProcessError as e:
    print(f"Error running MultiQC: {e}")

