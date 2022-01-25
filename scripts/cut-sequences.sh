#!/bin/bash
set -o nounset
set -o errexit

if [ "${adapter}" != "NA" ]; then
	echo "[PIPELINE -- cut-sequences]: Removing the adapter with Cutadapt..."
	echo "[PIPELINE -- cut-sequences]: Processing file: ${file}"
	
	docker run --rm \
		-v ${workingDir}:${workingDir} \
		-v ${fastqDir}:${fastqDir} \
		pegi3s/cutadapt:${cutadaptVersion} \
		-m 1 \
		-a ${adapter} \
		-o ${workingDir}/${outDir}/${ctdOut}/trimmed_${file} \
		${fastqDir}/${file}
else
	echo "[PIPELINE -- cut-sequences]: No adapter specified, adapter removal skipped"
fi
