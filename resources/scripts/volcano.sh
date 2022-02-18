#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- volcano]: Building the Volcano Plot..."

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${enhancedVolcanoRscript} ${workingDir}/compi_scripts/${enhancedVolcanoRscript}

# get the contrast name to build the output filename
contrast_label=$(echo "${comparison}" | cut -d'=' -f1 | tr -d \" | xargs)

if [[ ${selectDEAsoftware} == 'edger' ]]
then
	path_output_docker="${workingDir}/${outDir}/${edgOut}/"
	vp_path_counts="${path_output_docker}/$(echo EdgeR_${contrast_label} | xargs).tsv"
	vp_software="edger"
elif [[ ${selectDEAsoftware} == 'deseq' ]]
then
	path_output_docker="${workingDir}/${outDir}/${dsqOut}/"
	vp_path_counts="${path_output_docker}/$(echo DESeq2_${contrast_label} | xargs).tsv"
	vp_software="deseq"
elif [[ ${selectDEAsoftware} = 'both' ]]
then
	pipel_dir_name='pipel'
	path_output_docker="${workingDir}/${outDir}/${deaIntOut}/${contrast_label}/"
	vp_path_counts="${path_output_docker}/DEmiRNAs_${contrast_label}_deseq-edger_integrated.tsv"
	vp_software="both"
fi

vp_path_output=${path_output_docker}

docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_enhanced-volcano:${rEnhancedVolcanoVersion} \
		Rscript ${workingDir}/compi_scripts/${enhancedVolcanoRscript} ${vp_path_counts} ${path_output_docker} ${contrast_label} ${vp_software}
