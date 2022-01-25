#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- initialization]: Building the directory tree..."

#Function for the directory creation: first arg=${xOut}, second arg=name of the software
build_dir () {
	if [ ! -d "${workingDir}/${outDir}/$1" ]; then
		echo "[PIPELINE -- initialization]: Creating the directory for the $2 results"
		mkdir -p ${workingDir}/${outDir}/$1/
	else
		echo "[PIPELINE -- initialization]: $1 already exist"; echo "[PIPELINE -- initialization]: Removing ${workingDir}/${outDir}/$1 directory and files..."
		rm -Rf ${workingDir}/${outDir}/$1/
		echo "[PIPELINE -- initialization]: Creating the directory for the $2 results"
		mkdir -p ${workingDir}/${outDir}/$1/
	fi
}

#FastQC folder
build_dir ${fqcOut} "FastQC"

#Cutadapt folder
build_dir ${ctdOut} "Cutadapt"

#Bowtie folder
build_dir ${bwtOut} "Bowtie"

#Samtools stats folder
build_dir ${bamstOut} "BAM stats"

#FeatureCounts folder
build_dir ${ftqOut} "FeatureCounts"

#DESeq2 folder
build_dir ${dsqOut} "DESeq2"

#EdgeR folder
build_dir ${edgOut} "EdgeR"

#DEA integration folder
build_dir ${deaIntOut} "DEA integration"

rm -rf ${workingDir}/compi_scripts
mkdir -p ${workingDir}/compi_scripts
