#!/bin/bash
set -o nounset
set -o errexit

if [ "${adapter}" != "NA" ]; then
	echo "[MBS | cut-sequences]: Removing the adapter with Cutadapt..."
	echo "[MBS | cut-sequences]: Processing file: ${file}"
	
	docker run --rm \
		-v ${workingDir}:${workingDir} \
		-v ${fastqDir}:${fastqDir} \
		pegi3s/cutadapt:${cutadaptVersion} \
		-m 1 \
		-a ${adapter} \
		-o ${workingDir}/${outDir}/${ctdOut}/trimmed_${file} \
		${fastqDir}/${file}
else
	echo "[MBS | cut-sequences]: No adapter specified, adapter removal skipped"
fi
