#!/bin/bash
set -o nounset
set -o errexit

SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh

# lock Rscript before copying to avoid errors when parallel tasks are running
cp_and_lock ${enhancedVolcanoRscript} 'volcano' ${scriptsDir}

# get the contrast name to build the output filename
contrast_label=$(echo "${comparison}" | cut -d'"' -f2 )

# function to perform the volcano analysis
function run_volcano {
# $1 : ${vp_path_counts}
# $2 : ${path_output_docker}
# $3 : ${contrast_label}
# $4 : ${vp_software}
echo "[MBS | volcano]: Building the Volcano Plot of ${4} results..."
  docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_enhanced-volcano:${rEnhancedVolcanoVersion} \
		Rscript "${workingDir}/compi_scripts/${enhancedVolcanoRscript}" "${1}" "${2}" "${3}" "${4}"
}

# function to test if dea file exists and run analysis
function test_and_run() {
# $1 : ${vp_path_counts}
# $2 : ${path_output_docker}
# $3 : ${contrast_label}
# $4 : ${vp_software}
# $5 : Text "DESeq2", "EdgeR" or "integrated".
if [[ -f "$1" ]]
then
	run_volcano "$1" "$2" "$3" "$4"
else
	echo "[MBS | volcano | $4]: No $5 results."
	echo "[MBS | volcano | $4]: Done."
fi
}

# perform the volcano plot on the corresponding data
if [[ ${selectDEAsoftware} == 'edger' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker=$(get_output_dir edger ${comparison})
	vp_path_counts="${path_output_docker}/$(echo EdgeR_${contrast_label} | xargs).tsv"
	vp_software="edger"
	test_and_run "${vp_path_counts}" "${path_output_docker}" "${contrast_label}" "${vp_software}" "EdgeR"
fi

if [[ ${selectDEAsoftware} == 'deseq' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output_docker=$(get_output_dir deseq ${comparison})
	vp_path_counts="${path_output_docker}/$(echo DESeq2_${contrast_label} | xargs).tsv"
	vp_software="deseq"
	test_and_run "${vp_path_counts}" "${path_output_docker}" "${contrast_label}" "${vp_software}" "DESeq2"
fi

if [[ ${selectDEAsoftware} == 'both' ]]
then
	contrast_name=$(echo "${comparison}" | cut -d'"' -f4 )
	pipel_dir_name='pipel'
	path_output_docker="${workingDir}/${outDir}/${deaIntOut}/${contrast_label}/"
	vp_path_counts="${path_output_docker}/DEmiRNAs_${contrast_name}_deseq-edger_integrated.tsv"
	vp_software="both"
	test_and_run "${vp_path_counts}" "${path_output_docker}" "${contrast_label}" "${vp_software}" "integrated"
fi

#vp_path_output=${path_output_docker}
