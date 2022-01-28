#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- build-genome-index]: Building the genome index for Bowtie..."
echo "[PIPELINE -- build-genome-index]: Genome path: ${genome}"

bgi_genome="${workingDir}/input/genome/${genome}"
bgi_wd_genome="${workingDir}/input/genome/${genome}"

#makes the path to the directory using the genome filename as part of the directory name
bgi_output="${workingDir}/input/genome/${genome}/bowtie-index_${genome}"
bgi_output_wd="${workingDir}/input/genome/${genome}/bowtie-index_${genome}"

#first argument passed to bowtie (files on ${genome} separated by comma)
bgi_files="$(find ${bgi_wd_genome} -maxdepth 1 -type f | rev | cut -d\/ -f 1 | rev | tr '\n' ',' | sed 's/,$//')"

#second argument passed to bowtie
bgi_index_path="${bgi_output}/${genome}"

#if the genome directory exists skip the index creation
if [[ -d ${bgi_output_wd} ]]
then
	echo "[PIPELINE -- build-genome-index]: Genome index already exists in: ${bgi_output_wd}"
	echo "[PIPELINE -- build-genome-index]: Skipping Bowtie1 index creation"
else
	echo "[PIPELINE -- build-genome-index]:   Path to genome file:           ${bgi_genome}"
	echo "[PIPELINE -- build-genome-index]:   Output path for the index:     ${bgi_output}"
	echo "[PIPELINE -- build-genome-index]:   Root name for the index files: ${genome}"

	echo "[PIPELINE -- build-genome-index]: Running bowtie for the genome index creation"
	docker run --rm \
		-v ${workingDir}:${workingDir} \
		pegi3s/bowtie1:${bowtieVersion} \
		sh -c "rm -rf ${bgi_output}\
		&& mkdir ${bgi_output} \
		&& cd ${bgi_genome} \
		&& bowtie-build -f ${bgi_files} ${bgi_index_path}"
fi
