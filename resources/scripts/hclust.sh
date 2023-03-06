#!/bin/bash
set -o nounset
set -o errexit

SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh

# lock Rscript before copying to avoid errors when parallel tasks are running
cp_and_lock ${hclustMakeTableRscript} 'hclust' ${scriptsDir}
cp_and_lock ${hclustRscript} 'hclust' ${scriptsDir}
cp_and_lock ${deseq2NormalizationRscript} 'hclust' ${scriptsDir}

# get the contrast name to build the output filename
contrast_label=$(echo "${comparison}" | cut -d'"' -f2 )
contrast_label_samples=$(echo "${comparison}" | cut -d'"' -f4 ) # FIXME: change the name of the integrated dea results, modify the pipeline accordingly

# function to perform the global count normalization with DESeq2
function run_deseq2_normalization {
# $1 : ${path_counts}       $4 : ${software}
# $2 : ${path_conditions}   $5 : ${input_contrast}
# $3 : ${path_output} 
echo "[PIPELINE -- hclust]: Performing the global count normalization with DESeq2..."
  docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_deseq2:${rDeseq2Version} \
		Rscript "${workingDir}/compi_scripts/${deseq2NormalizationRscript}" "${1}" "${2}" "${3}" "${4}" "${5}"
}

# function to perform the hierarchical clustering analysis
function run_hclust_make-table {
# $1 : ${dea_results}     $4 : ${input_contrast}
# $2 : ${path_counts}     $5 : ${path_output_docker}
# $3 : ${path_conditions} $6 : ${software}
echo "[PIPELINE -- hclust]: Building the hierarchical clustering table for ${4} results..."
  docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_data-analysis:${rdatanalysisVersion} \
		Rscript "${workingDir}/compi_scripts/${hclustMakeTableRscript}" "${1}" "${2}" "${3}" "${4}" "${5}" "${6}"
}

# function to run the hclust analysis with the Rscript
function run_hclust {
# $1 : ${path_hclust_file}  $3 : ${path_output}  $5 : ${conditions_file}
# $2 : ${input_contrast}    $4 : ${software}
echo "[PIPELINE -- hclust]: Running the hierarchical clustering analysis for ${4} results..."
  docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_data-analysis:${rdatanalysisVersion} \
		Rscript "${workingDir}/compi_scripts/${hclustRscript}" "${1}" "${2}" "${3}" "${4}" "${5}"
}

# function to test if dea file exists and run the "hclust" Rscript
function test_and_run_hclust {
# $1 : ${path_hclust_file}  $3 : ${path_output}  $5 : ${conditions_file}
# $2 : ${input_contrast}    $4 : ${software}
if [[ -f "$1" ]]
then
	run_hclust "${1}" "${2}" "${3}" "${4}" "${5}"
else
	echo "[PIPELINE -- hclust > $4]: No hclust table."
	echo "[PIPELINE -- hclust > $4]: Skipping."
fi
}

# function to test if dea file exists and run all the Rscripts
function test_and_run() {
# $1 : ${dea_results}         $6 : ${software}
# $2 : ${path_counts}         $7 : Text "DESeq2", "EdgeR" or "integrated"
# $3 : ${path_conditions}     $8 : ${path_hclust_file}
# $4 : ${input_contrast}      $9 : ${path_output_pipel}
# $5 : ${path_output_docker}
if [[ -f "$1" ]]
then
	run_hclust_make-table "${1}" "${2}" "${3}" "${4}" "${9}" "${6}"
	test_and_run_hclust "${8}" "${4}" "${5}" "${6}" "${3}"
else
	echo "[PIPELINE -- hclust > $6]: No $7 results."
	echo "[PIPELINE -- hclust > $6]: Done."
fi
}

path_conditions="${conditions}"
input_contrast="${comparison}"

# check if a DESeq2 global count normalization should be performed
if [[ -v ${1} ]]
then
	path_counts="${workingDir}/${outDir}/${ftqOut}/all-counts.txt"

	if [[ ${selectDEAsoftware} == 'edger' ]] || [[ ${selectDEAsoftware} == 'both' ]]
	then
		path_output_normalized_counts="$(get_output_dir edger ${comparison})/pipel/"
		software="EdgeR"
		run_deseq2_normalization "${path_counts}" "${path_conditions}" "${path_output_normalized_counts}" "${software}" "${input_contrast}"
	fi

	if [[ ${selectDEAsoftware} == 'deseq' ]] || [[ ${selectDEAsoftware} == 'both' ]]
	then
		path_output_normalized_counts="$(get_output_dir deseq ${comparison})/pipel/"
		software="DESeq2"
		run_deseq2_normalization "${path_counts}" "${path_conditions}" "${path_output_normalized_counts}" "${software}" "${input_contrast}"
	fi

	if [[ ${selectDEAsoftware} == 'both' ]]
	then
		path_output_normalized_counts="${workingDir}/${outDir}/${deaIntOut}/${contrast_label}/pipel/"
		software="DESeq2-EdgeR"
		run_deseq2_normalization "${path_counts}" "${path_conditions}" "${path_output_normalized_counts}" "${software}" "${input_contrast}"
	fi

fi

# perform the hierarchical clustering on the corresponding data
if [[ ${selectDEAsoftware} == 'edger' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker=$(get_output_dir edger ${comparison})
	path_output_pipel="${path_output_docker}/pipel/"
	dea_results="${path_output_docker}/$(echo EdgeR_${contrast_label} | xargs).tsv"
	software="EdgeR"

	if [[ -v ${1} ]]; then
		path_counts="${path_output_pipel}/all-counts_${software}_${contrast_label}_normalized.txt"
	else
		path_counts="${path_output_docker}/pipel/$(echo all-counts_edger_${contrast_label} | xargs).txt"
	fi

	path_hclust_file="${path_output_pipel}/$(echo hclust_EdgeR_${contrast_label} | xargs).tsv"
	# run the analysis
	test_and_run "${dea_results}" "${path_counts}" "${path_conditions}" "${input_contrast}" "${path_output_docker}" "${software}" 'EdgeR' "${path_hclust_file}" "${path_output_pipel}"
fi

if [[ ${selectDEAsoftware} == 'deseq' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker=$(get_output_dir deseq ${comparison})
	path_output_pipel="${path_output_docker}/pipel/"
	dea_results="${path_output_docker}/$(echo DESeq2_${contrast_label} | xargs).tsv"
	software="DESeq2"

	if [[ -v ${1} ]]; then
		path_counts="${path_output_pipel}/all-counts_${software}_${contrast_label}_normalized.txt"
	else
		path_counts="${path_output_docker}/pipel/$(echo all-counts_deseq_${contrast_label} | xargs).txt"
	fi

	path_hclust_file="${path_output_pipel}/$(echo hclust_DESeq2_${contrast_label} | xargs).tsv"
	# run the analysis
	test_and_run "${dea_results}" "${path_counts}" "${path_conditions}" "${input_contrast}" "${path_output_docker}" "${software}" 'EdgeR' "${path_hclust_file}" "${path_output_pipel}"
fi

if [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker="${workingDir}/${outDir}/${deaIntOut}/${contrast_label}/"
	path_output_pipel="${path_output_docker}/pipel/"
	dea_results="${path_output_docker}/DEmiRNAs_${contrast_label_samples}_deseq-edger_integrated.tsv"
	software="DESeq2-EdgeR"

	# uses the counts of EdgeR (data is the same for both software)
	if [[ -v ${1} ]]; then
		path_counts="${path_output_pipel}/all-counts_${software}_${contrast_label}_normalized.txt"
	else
		path_counts="${path_output_docker}/pipel/$(echo all-counts_edger_${contrast_label} | xargs).txt"
	fi

	path_hclust_file="${path_output_pipel}/$(echo hclust_DESeq2-EdgeR_${contrast_label} | xargs).tsv"
	# run the analysis
	test_and_run "${dea_results}" "${path_counts}" "${path_conditions}" "${input_contrast}" "${path_output_docker}" "${software}" 'EdgeR' "${path_hclust_file}" "${path_output_pipel}"
fi
