#!/bin/bash
set -o nounset
set -o errexit

## Script input parameters:
##  1.- path_tarbase_db           : a file downloaded from Tarbase.
##  2.- path_reactome_db          : a file downloaded from Reactome database.
##  3.- path_enrichment_table     : the file resulting from a DEA of the pipeline.
##  4.- path_reactome_interaction : a file downloaded from Reactome database.
##  5.- input_contrast            : a line of the contrast_file.
##  6.- software                  : the software that produced the "dea_result_file" (DESeq2/EdgeR/DESeq2-EdgeR).
##  7.- db_organism               : the name of the organism, e.g.: 'Homo sapiens'
##  8.- path_output               : path to save the output files.

echo "[MBS | network]: Starting the network creation..."

SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh

# lock Rscript before copying to avoid errors when parallel tasks are running
cp_and_lock ${networkRscript} 'network' ${scriptsDir}
cp_and_lock ${tarbaseDB} 'network' ${databasesDir}
cp_and_lock ${reactomeDB} 'network' ${databasesDir}
cp_and_lock ${reactomeInteractionsDB} 'network' ${databasesDir}

# function to perform the network analysis
function run_functional_enrichment {
# $1 : ${path_tarbase_db}             $5 : ${input_contrast}
# $2 : ${path_reactome_db}            $6 : ${software}
# $3 : ${path_enrichment_table}       $7 : ${db_organism}
# $4 : ${path_reactome_interaction}   $8 : ${path_output}
echo "[MBS | network | ${6}]: Creating the network of ${6} results..."
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_network:${rNetworkVersion} \
	Rscript ${workingDir}/compi_scripts/${networkRscript} "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"
rm -rf "${8}"/network__R-HSA-*[!\.html]
}

# function to test if dea file exists and run analysis
function test_and_run() {
# $1 : ${path_tarbase_db}             $5 : ${input_contrast}   $9 : Text "DESeq2", "EdgeR" or "integrated".
# $2 : ${path_reactome_db}            $6 : ${software}
# $3 : ${path_enrichment_table}       $7 : ${db_organism}
# $4 : ${path_reactome_interaction}   $8 : ${path_output}
if [[ -f "$3" ]]
then
	run_functional_enrichment "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"
else
	echo "[MBS | network | $6]: No $9 results."
	echo "[MBS | network | $6]: Done."
fi
}

# get the contrast label and name to build the filenames
contrast_label=$(echo "${comparison}" | cut -d'"' -f2 )
contrast_name=$(echo "${comparison}" | cut -d'"' -f4 )

# inputs
path_tarbase_db="${workingDir}/compi_databases/${tarbaseDB}"
path_reactome_db="${workingDir}/compi_databases/${reactomeDB}"
path_reactome_interaction="${workingDir}/compi_databases/${reactomeInteractionsDB}"
input_contrast="${comparison}"
db_organism="${organism}"

# Build the inputs and run
if [[ ${selectDEAsoftware} == 'edger' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output=$(get_output_dir edger ${comparison})
	software="EdgeR"
	path_enrichment_table="${path_output}/$(echo enrichment_table_${software}_${contrast_name} | xargs).tsv"
	test_and_run "${path_tarbase_db}" "${path_reactome_db}" "${path_enrichment_table}" "${path_reactome_interaction}" "${input_contrast}" "${software}" "${db_organism}" "${path_output}" 'EdgeR'
fi

if [[ ${selectDEAsoftware} == 'deseq' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output=$(get_output_dir deseq ${comparison})
	software="DESeq2"
	path_enrichment_table="${path_output}/$(echo enrichment_table_${software}_${contrast_name} | xargs).tsv"
	test_and_run "${path_tarbase_db}" "${path_reactome_db}" "${path_enrichment_table}" "${path_reactome_interaction}" "${input_contrast}" "${software}" "${db_organism}" "${path_output}" 'DESeq2'
fi

if [[ ${selectDEAsoftware} == 'both' ]]
then
	#get the contrast name to find the integrated file
	path_output="${workingDir}/${outDir}/${deaIntOut}/${contrast_label}/"
	path_enrichment_table="${path_output}/enrichment_table_DESeq2-EdgeR_${contrast_name}.tsv"
	software='DESeq2-EdgeR'
	test_and_run "${path_tarbase_db}" "${path_reactome_db}" "${path_enrichment_table}" "${path_reactome_interaction}" "${input_contrast}" "${software}" "${db_organism}" "${path_output}" 'integrated'
fi
