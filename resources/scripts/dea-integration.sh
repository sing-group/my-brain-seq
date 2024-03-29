#!/bin/bash
set -o nounset
set -o errexit

echo "[MBS | dea-integration]: Integrating DESeq2 and EdgeR results..."

SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh

# lock Rscript before copying to avoid errors when parallel tasks are running
cp_and_lock ${deaIntRscript} 'dea-integration' ${scriptsDir}

# find the filenames of DESeq2 and EdgeR results
echo "[MBS | dea-integration]: Finding DESeq2 and EdgeR results"
contrast_filename=$(echo "${comparison}" | cut -d'=' -f1 | tr -d \" | sed -e 's/ /\\ /g')
vp_comparison_label="$(echo ${contrast_filename} | xargs)"

pipel_dir_name='pipel'

## Script input parameters:
##  1.- path_deseq        : the file resulting from the DESeq2 analysis.
##  2.- path_edger        : the file resulting from the EdgeR analysis.
##  3.- input_contrast    : a line of the contrast_file.
##  4.- path_output       : output path to save the results.
##  5.- path_output_pipel : output path to save additional results.

# Inputs of the R script
path_deseq="${workingDir}/${outDir}/${dsqOut}/${vp_comparison_label}/DESeq2_${vp_comparison_label}.tsv"
path_edger="${workingDir}/${outDir}/${edgOut}/${vp_comparison_label}/EdgeR_${vp_comparison_label}.tsv"
input_contrast="${comparison}"
path_output="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}"
path_output_pipel="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}/${pipel_dir_name}"

# make the directory for the results of the contrast and the pipeline files
mkdir -p "${path_output_pipel}"

echo "[MBS | dea-integration]:	DESeq2 file: $path_deseq"
echo "[MBS | dea-integration]:	EdgeR  file: $path_edger"

echo "[MBS | dea-integration]: Running DESeq2-EdgeR integration..."

docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_data-analysis \
		Rscript ${workingDir}/compi_scripts/${deaIntRscript} "${path_deseq}" "${path_edger}" "${input_contrast}" "${path_output}" "${path_output_pipel}" "${qvalue}" "${log2FC}"

echo "[MBS | dea-integration]: Done!"
