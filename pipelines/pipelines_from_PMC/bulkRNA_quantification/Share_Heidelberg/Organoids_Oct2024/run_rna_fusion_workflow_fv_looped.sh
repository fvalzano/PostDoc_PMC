#!/usr/bin/env bash

directories=()
cd /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/ubam2counts/ubam2counts_HD_Oct2024
for dir in */; do
    # Add the full path of each directory to the 'directories' array
    directories+=("$(pwd)/${dir}")
done 

cd /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/cromwell
workflow_source='rna_fusions_germline_snv_workflow.wdl'
workflow_dependencies="*.wdl"

set -o errexit
set -o pipefail

source ./utils/get_cromwell_credentials.sh
source ./utils/starting_cromwell_jobs.sh

cd ../wdl
rm -f modules.zip
zip -q modules.zip $workflow_dependencies

for dir in "${directories[@]}"; do
    workflow_inputs=${dir}/rna_fusions_germline_snv_inputs_no_molgenis_fv.json
    echo "Submitting workflow for ${workflow_inputs}"
    _caas_submit \
        "$workflow_source" \
        "$workflow_inputs" \
        "modules.zip"
done