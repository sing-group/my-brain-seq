#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- bam-stats]: Analyzing bam files with SAMtools stats..."

bam_input="${workingDir}/${outDir}/${bwtOut}/${bamFile}"
bs_output="${workingDir}/${outDir}/${bamstOut}/$(echo ${bamFile} | cut -d'.' -f1)"
bs_stats_file="${bs_output}/${bamFile}.stats"
bs_plot_names="${bs_output}/${bamFile}_PLOT"

mkdir "${bs_output}"

docker run --rm\
	-v ${workingDir}:${workingDir} \
	pegi3s/samtools_bcftools:${samtoolsVersion} \
	sh -c "samtools stats ${bam_input} > ${bs_stats_file}"

docker run --rm\
	-v ${workingDir}:${workingDir} \
	pegi3s/samtools_bcftools:${samtoolsVersion} \
	sh -c "plot-bamstats -p ${bs_plot_names} ${bs_stats_file}"
