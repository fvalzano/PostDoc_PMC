cd /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fasta_inputs
for dir in */; do
    cd $dir
    for sample in `ls`; do
        mkdir -p /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fastq2ubam/fastq2ubams_ITCCP4/$dir$sample
        output_file=/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fastq2ubam/fastq2ubams_ITCCP4/$dir$sample/fastq_ubam_workflow_fv.inputs.json
        json_content='{
    "MultiFastqUbamWorkflow.image_picard_digest": "docker.io/princessmaximacenter/picard@sha256:d561453a4a25dfac7933ec2237259d1c5adc974287cc639418bd3c2af9881019",
    "MultiFastqUbamWorkflow.sample_details": [
        {
            "sample_name": "", #SampleID
            "read_group_name": "001",
            "library_name": "", #SampleID
            "barcode": "", #path to barcode file
            "lanes_list": "" #path to .fasta.input file
        }
        # Add as many samples as you need ---Important - change the library_name and read_group_name---
    ],
    "MultiFastqUbamWorkflow.data_host": "hpct05",
    "MultiFastqUbamWorkflow.target_dir": "/hpc/pmc_kool/fvalzano/",  # Change Target directory to your favourite directory
    "MultiFastqUbamWorkflow.run_name": "", #SampleID
    "MultiFastqUbamWorkflow.machine_name": "Some_sequencer",
    "MultiFastqUbamWorkflow.platform_unit": "", #SampleID
    "MultiFastqUbamWorkflow.platform": "Illumina",
    "MultiFastqUbamWorkflow.sequencing_center": "Heidelberg",
    "MultiFastqUbamWorkflow.tmpspace_gb": 140,
    "MultiFastqUbamWorkflow.wallclock": "06:00:00",
    "MultiFastqUbamWorkflow.cache_path": "/hpc/pmc_kool/fvalzano/ubams/"  # Change cache_path to your favourite directory
}'

        # Write the JSON content to the file
        echo "$json_content" > "$output_file"
    done
    cd ..
done