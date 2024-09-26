#!/usr/bin/env bash

##This belongs to the orginal script - flag it and explicit directories 
#echo "${inputs_json}"
#
#if [ "${inputs_json}" == "" ]; then
#        echo "ERROR:No inputs.json provided\n"
#        echo "inputs_json=<INPUTS.JSON> $0"
#        exit 1
#fi

#This outputs the directory names to feed in the subsequent loop of submission to the CaaS
directories=()
cd /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fastq2ubam/fastq2ubams_AZ
for dir in */; do
    # Add the full path of each directory to the 'directories' array
    directories+=("$(pwd)/${dir}")
done 

cd /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/cromwell
workflow_source='/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/multi_fastq_ubam_workflow.wdl'
workflow_dependencies="*.wdl"

set -o errexit
set -o pipefail

source ./utils/get_cromwell_credentials.sh
source ./utils/starting_cromwell_jobs.sh

cd ../wdl
rm -f modules.zip
zip -q modules.zip $workflow_dependencies

for dir in "${directories[@]}"; do
    workflow_inputs=${dir}/fastq_ubam_workflow_fv.inputs.json
    _caas_submit \
        "$workflow_source" \
        "$workflow_inputs" \
        "modules.zip"
done