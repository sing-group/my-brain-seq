#!/bin/bash
set -o nounset
set -o errexit

echo "[MBS | fastqc-qc]: Performing a quality control with FastQC..."
echo "[MBS | fastqc-qc]: Processing file: ${file}"

fqc_input="${fastqDir}/${file}"
fqc_output="${workingDir}/${outDir}/${fqcOut}/"

docker run --rm \
	-v ${workingDir}:${workingDir} \
	-v ${fastqDir}:${fastqDir} \
	pegi3s/fastqc:${fastqcVersion} \
	-o ${fqc_output} \
	${fqc_input}
