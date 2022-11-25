#!/bin/bash
set -o nounset
set -o errexit

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${hclustMakeTableRscript} ${workingDir}/compi_scripts/${hclustMakeTableRscript}

# get the contrast name to build the output filename
contrast_label=$(echo "${comparison}" | cut -d'"' -f2 )

# function to perform the hierarchical clustering analysis
function run_hclust_make-table {
# $1 : ${dea_results}
# $2 : ${path_counts}
# $3 : ${path_conditions}
# $4 : ${input_contrast}
# $5 : ${path_output_docker}
# $6 : ${software}

echo "[PIPELINE -- hclust]: Building the hierarchical clustering table for ${4} results..."
  docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_data-analysis:${rdatanalysisVersion} \
		Rscript "${workingDir}/compi_scripts/${hclustMakeTableRscript}" "${1}" "${2}" "${3}" "${4}" "${5}" "${6}"
}

function run_hclust {
# $1 : ${path_hclust_file}
# $2 : ${input_contrast}
# $3 : ${path_output}
# $4 : ${software}

echo "[PIPELINE -- hclust]: Running the hierarchical clustering analysis for ${4} results..."
  docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_data-analysis:${rdatanalysisVersion} \
		Rscript "${workingDir}/compi_scripts/${hclustRscript}" "${1}" "${2}" "${3}" "${4}"
}

path_conditions="${conditions}"
input_contrast="${comparison}"

# perform the hierarchical clustering on the corresponding data
if [[ ${selectDEAsoftware} == 'edger' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker="${workingDir}/${outDir}/${edgOut}/"
	dea_results="${path_output_docker}/$(echo EdgeR_${contrast_label} | xargs).tsv"
	path_counts="${path_output_docker}/pipel/$(echo all-counts_edger_${contrast_label} | xargs).txt"
	software="EdgeR"
	path_hclust_file="${path_output_docker}/$(echo hclust_EdgeR_${contrast_label} | xargs).tsv"
	# prepare tables
	run_hclust_make-table "${dea_results}" "${path_counts}" "${path_conditions}" "${input_contrast}" "${path_output_docker}" "${software}"
	# hierarchical clustering
	run_hclust "${path_hclust_file}" "${input_contrast}" "${path_output_docker}" "${software}"
fi

if [[ ${selectDEAsoftware} == 'deseq' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker="${workingDir}/${outDir}/${dsqOut}/"
	dea_results="${path_output_docker}/$(echo DESeq2_${contrast_label} | xargs).tsv"
	path_counts="${path_output_docker}/pipel/$(echo all-counts_deseq_${contrast_label} | xargs).txt"
	software="DESeq2"
	path_hclust_file="${path_output_docker}/$(echo hclust_DESeq2_${contrast_label} | xargs).tsv"
	# prepare tables
	run_hclust_make-table "${dea_results}" "${path_counts}" "${path_conditions}" "${input_contrast}" "${path_output_docker}" "${software}"
	# hierarchical clustering
	run_hclust "${path_hclust_file}" "${input_contrast}" "${path_output_docker}" "${software}"
fi

if [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker="${workingDir}/${outDir}/${deaIntOut}/${contrast_label}/"
	dea_results="${path_output_docker}/DEmiRNAs_${contrast_label}_deseq-edger_integrated.tsv"
	# uses the counts of EdgeR DESeq2 and
	path_counts="${workingDir}/${outDir}/${edgOut}/pipel/$(echo all-counts_edger_${contrast_label} | xargs).txt"
	software="DESeq2-EdgeR"
	path_hclust_file="${path_output_docker}/$(echo hclust_DESeq2-EdgeR_${contrast_label} | xargs).tsv"
	# prepare tables
	run_hclust_make-table "${dea_results}" "${path_counts}" "${path_conditions}" "${input_contrast}" "${path_output_docker}" "${software}"
	# hierarchical clustering
	run_hclust "${path_hclust_file}" "${input_contrast}" "${path_output_docker}" "${software}"
fi
