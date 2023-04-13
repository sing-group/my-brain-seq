#!/bin/bash
set -o nounset
set -o errexit

echo "[MBS | build-genome-index]: Building the genome index for Bowtie..."
echo "[MBS | build-genome-index]: Genome path: ${genome}"

genome_name="$(basename ${genome} | cut -d '.' -f1)"

#makes the path to the directory using the genome filename
bgi_output="${workingDir}/input/genome/${genome_name}/bowtie-index_${genome_name}"

#second argument passed to bowtie
bgi_index_path="${bgi_output}/${genome_name}"

#if the genome directory exists skip the index creation
if [[ -d ${bgi_output} ]]
then
	echo "[MBS | build-genome-index]: Genome index already exists in: ${bgi_output}"
	echo "[MBS | build-genome-index]: Skipping Bowtie1 index creation"
else
	echo "[MBS | build-genome-index]:   Path to genome file:           ${genome}"
	echo "[MBS | build-genome-index]:   Output path for the index:     ${bgi_output}"
	echo "[MBS | build-genome-index]:   Root name for the index files: ${genome}"

	echo "[MBS | build-genome-index]: Running bowtie for the genome index creation"
	docker run --rm \
		-v ${workingDir}:${workingDir} \
		-v ${genome}:${genome} \
		-v ${bgi_output}:${bgi_output} \
		pegi3s/bowtie1:${bowtieVersion} \
		sh -c "mkdir -p ${bgi_output} \
		&& cd ${bgi_output} \
		&& bowtie-build -f ${genome} ${bgi_index_path}"
fi
