#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- fastqc-qc]: Performing a quality control with FastQC..."
echo "[PIPELINE -- fastqc-qc]: Processing file: ${file}"

fqc_input="${fastqDir}/${file}"
fqc_output="${workingDir}/${outDir}/${fqcOut}/"

docker run --rm \
	-v ${workingDir}:${workingDir} \
	-v ${fastqDir}:${fastqDir} \
	pegi3s/fastqc:${fastqcVersion} \
	-o ${fqc_output} \
	${fqc_input}
