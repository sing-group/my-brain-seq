#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- multiqc]: Starting the creation of the MultiQC report..."

# prepare the paths for the report
results_path="${workingDir}/${outDir}/"
output_dir="${workingDir}/${outDir}/${mqcOut}/"

# run multiqc
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/multiqc:${multiqcVersion} \
	multiqc "${results_path}" --outdir "${output_dir}"

echo "[PIPELINE -- multiqc]: Done"