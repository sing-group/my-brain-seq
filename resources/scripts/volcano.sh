#!/bin/bash
set -o nounset
set -o errexit

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${enhancedVolcanoRscript} ${workingDir}/compi_scripts/${enhancedVolcanoRscript}

# get the contrast name to build the output filename
contrast_label=$(echo "${comparison}" | cut -d'"' -f2 )

# function to perform the volcano analysis
function run_volcano {
# $1 : ${vp_path_counts}
# $2 : ${path_output_docker}
# $3 : ${contrast_label}
# $4 : ${vp_software}
echo "[PIPELINE -- volcano]: Building the Volcano Plot of ${4} results..."
  docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_enhanced-volcano:${rEnhancedVolcanoVersion} \
		Rscript "${workingDir}/compi_scripts/${enhancedVolcanoRscript}" "${1}" "${2}" "${3}" "${4}"
}

# perform the volcano plot on the corresponding data
if [[ ${selectDEAsoftware} == 'edger' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker="${workingDir}/${outDir}/${edgOut}/"
	vp_path_counts="${path_output_docker}/$(echo EdgeR_${contrast_label} | xargs).tsv"
	vp_software="edger"
	run_volcano "${vp_path_counts}" "${path_output_docker}" "${contrast_label}" "${vp_software}"
fi

if [[ ${selectDEAsoftware} == 'deseq' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker="${workingDir}/${outDir}/${dsqOut}/"
	vp_path_counts="${path_output_docker}/$(echo DESeq2_${contrast_label} | xargs).tsv"
	vp_software="deseq"
	run_volcano "${vp_path_counts}" "${path_output_docker}" "${contrast_label}" "${vp_software}"
fi

if [[ ${selectDEAsoftware} = 'both' ]]
then
	pipel_dir_name='pipel'
	path_output_docker="${workingDir}/${outDir}/${deaIntOut}/${contrast_label}/"
	vp_path_counts="${path_output_docker}/DEmiRNAs_${contrast_label}_deseq-edger_integrated.tsv"
	vp_software="both"
	run_volcano "${vp_path_counts}" "${path_output_docker}" "${contrast_label}" "${vp_software}"
fi

#vp_path_output=${path_output_docker}

