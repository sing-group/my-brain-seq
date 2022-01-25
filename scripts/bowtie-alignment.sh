#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- bowtie-alignment]: Performing genome alignment with Bowtie..."

#choose the path to the bowtie index (depending if the user provides a pre-built index or if it was built in the task "build-genome-index").
#if genome flag
if [ ! -z ${genome} ]; then
	bgi_genome="${workingDir}/input/genome/${genome}/"
	bw_index_path="${workingDir}/input/genome/${genome}/bowtie-index_${genome}/"
else
	bw_index_path=${bwtIndex}
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

#gets the root of the index file name
bw_index=$(ls -1 ${bw_index_path} | head -n1 | sed 's/\..\.ebwt//' | sed 's@\(.*\)@'${bw_index_path}'\/\1@')
bw_index=$(echo $bw_index | tr -s '/')
#bw_index=$(echo $bw_index | cut -d / -f 3- | awk '{print "${workingDir}/" $0}')
echo "[PIPELINE -- bowtie-alignment]: Index to align with: ${bw_index}"
	
#test if the file exist, if so removes it
if [ -f "${bw_output}" ]; then
	echo "[PIPELINE -- bowtie-alignment]: ${file}.sam already exists, removing ${file}.sam."
	rm -f "${bw_output}"
fi

echo "[PIPELINE -- bowtie-alignment]: Aligning file: ${file}.sam"

docker run --rm \
	-v ${workingDir}:${workingDir} \
	-v ${bwtIndex}:${bwtIndex} \
	-v ${fastqDir}:${fastqDir} \
	pegi3s/bowtie1:${bowtieVersion} \
	bowtie -S \
	${bw_index} \
	${bw_fastq} \
	${bw_output}
