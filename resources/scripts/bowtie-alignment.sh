#!/bin/bash
set -o nounset
set -o errexit

echo "[MBS | bowtie-alignment]: Performing genome alignment with Bowtie..."

#choose the path to the bowtie index (depending if the user provides a 
#pre-built index or if it was built in the task "build-genome-index").
#if genome flag
if [ ! -z ${genome} ]; then
	genome_name="$(basename ${genome} | cut -d '.' -f1)"
	bw_index_path="${workingDir}/input/genome/${genome_name}/bowtie-index_${genome_name}"
	bwtIndex="${bw_index_path}/${genome_name}"
elif [ ! -z ${bwtIndex} ]; then
	bw_index_path="$(dirname ${bwtIndex})"
else
	echo "[MBS | bowtie-alignment]: [ERROR] No genome or bowtie index provided."
fi

# if adapter specified
if [ "${adapter}" != "NA" ]; then
	bw_fastq="${workingDir}/${outDir}/${ctdOut}/${file}"
	
# if no adapter specified
else
	bw_fastq="${fastqDir}/${file}"
fi

bw_fastq=$(echo $bw_fastq | tr -s '/')  #convert the path to single slash

bw_output="${workingDir}/${outDir}/${bwtOut}/${file}.sam"
bw_output=$(echo $bw_output | tr -s '/')  #convert the path to single slash

echo "[MBS | bowtie-alignment]: Index to align with: ${bwtIndex}"
	
#test if the file exist, if so removes it
if [ -f "${bw_output}" ]; then
	echo "[MBS | bowtie-alignment]: ${file}.sam already exists, removing ${file}.sam."
	rm -f "${bw_output}"
fi

echo "[MBS | bowtie-alignment]: Aligning file: ${file}.sam"

docker run --rm \
	-v ${workingDir}:${workingDir} \
	-v ${bw_index_path}:${bw_index_path} \
	-v ${fastqDir}:${fastqDir} \
	pegi3s/bowtie1:${bowtieVersion} \
	bowtie -S \
	${bwtIndex} \
	${bw_fastq} \
	${bw_output}
