#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- venn]: Building a file for the Venn diagram"

SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh

# lock Rscript before copying to avoid errors when parallel tasks are running
cp_and_lock ${vennRscript} 'venn'

# find the filenames of filtered DESeq2 and EdgeR results
echo "[PIPELINE -- venn]: Finding filtered DESeq2 and EdgeR results"
vp_comparison_label=$(echo "${comparison}" | cut -d'"' -f2)
vp_comparison_label_file=$(echo "${comparison}" | cut -d'"' -f4)
#vp_comparison_label="$(echo $contrast_filename | xargs)"

# output dir
pipel_dir_name='pipel'
output="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}/"
output_pipel="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}/${pipel_dir_name}/"

## Script input parameters:
##  1.- venn_path: path to a file with the DE features returned by DESeq2 + EdgeR, one per line.
##  2.- venn_output: path to the output directory for the graph.
##  3.- venn_output_format: file format of the output image (png, svg or tiff).
##  4.- input_contrast : a line of the contrast_file.
#
# Inputs
venn_path=$(echo "${output_pipel}/DEmiRNAs_${vp_comparison_label_file}_deseq-edger_venn.tsv")
venn_output="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}/"
venn_output_format="${vennFormat}"
input_contrast="${comparison}"

# Run the R script
echo "[PIPELINE -- venn]: Building the Venn diagram"
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_venn-diagram \
		Rscript "${workingDir}/compi_scripts/${vennRscript}" "${venn_path}" "${venn_output}" "${venn_output_format}" "${input_contrast}"

