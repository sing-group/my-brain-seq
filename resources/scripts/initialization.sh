#!/bin/bash
set -o nounset
set -o errexit

echo "[MBS | initialization]: Building the directory tree..."

# function for the directory creation
build_dir () {
	# $1 : ${xOut}, pipeline parameter with the output directory.
	# $2 : Name of the software, e.g.: DESeq2.
	if [ ! -d "${workingDir}/${outDir}/$1" ]; then
		echo "[MBS | initialization]: Creating the directory for the $2 results"
		mkdir -p ${workingDir}/${outDir}/$1/
	else
		echo "[MBS | initialization]: $1 already exist"; echo "[MBS | initialization]: Removing ${workingDir}/${outDir}/$1 directory and files..."
		rm -Rf ${workingDir}/${outDir}/$1/
		echo "[MBS | initialization]: Creating the directory for the $2 results"
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

#MultiQC folder
build_dir ${mqcOut} "MultiQC"

rm -rf ${workingDir}/compi_scripts
mkdir -p ${workingDir}/compi_scripts

rm -rf ${workingDir}/compi_databases
mkdir -p ${workingDir}/compi_databases
