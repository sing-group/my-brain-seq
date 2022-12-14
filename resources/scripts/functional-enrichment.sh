#!/bin/bash
set -o nounset
set -o errexit

## Script input parameters
##  1.- tarbase_db_file  : a file downloaded from Tarbase.
##  2.- reactome_db_file : a file downloaded from Reactome database.
##  3.- dea_result_file  : the file resulting from a DEA of the pipeline.
##  4.- input_contrast   : a line of the contrast_file.
##  5.- software         : the software that produced the "dea_result_file" (DESeq2/EdgeR/DESeq2-EdgeR).
##  6.- db_organism      : the name of the organism, e.g.: 'Homo sapiens'
##  7.- path_output      : path to save the output files.

echo "[PIPELINE -- functional-enrichment]: Starting the functional enrichment analysis..."

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${functionalEnrichmentRscript} ${workingDir}/compi_scripts/${functionalEnrichmentRscript}
cp ${databasesDir}/${tarbaseDB} ${workingDir}/compi_databases/${tarbaseDB}
cp ${databasesDir}/${reactomeDB} ${workingDir}/compi_databases/${reactomeDB}

# function to perform the enrichment analysis
function run_functional_enrichment {
# $1 : ${tarbase_db_file}    $5 : ${software}
# $2 : ${reactome_db_file}   $6 : ${db_organism}
# $3 : ${dea_result_file}    $7 : ${path_output}
# $4 : ${input_contrast}
echo "[PIPELINE -- functional-enrichment > ${5}]: Performing the functional enrichment analysis of ${5} results..."
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_data-analysis:${rdatanalysisVersion} \
	Rscript ${workingDir}/compi_scripts/${functionalEnrichmentRscript} "$1" "$2" "$3" "$4" "$5" "$6" "$7"
}

# function to test if dea file exists and run analysis
function test_and_run() {
# $1 : ${tarbase_db_file}    $5 : ${software}
# $2 : ${reactome_db_file}   $6 : ${db_organism}
# $3 : ${dea_result_file}    $7 : ${path_output}
# $4 : ${input_contrast}     $8 : Text "DESeq2", "EdgeR" or "integrated".
if [[ -f "$3" ]]
then
	run_functional_enrichment "$1" "$2" "$3" "$4" "$5" "$6" "$7"
else
	echo "[PIPELINE -- functional-enrichment > $5]: No $8 results."
	echo "[PIPELINE -- functional-enrichment > $5]: Done."
fi
}

# get the contrast name to build the output filename
contrast_label=$(echo "${comparison}" | cut -d'"' -f2 )

# Inputs
tarbase_db_file="${workingDir}/compi_databases/${tarbaseDB}"
reactome_db_file="${workingDir}/compi_databases/${reactomeDB}"
input_contrast="${comparison}"
db_organism="${organism}"

# Build the inputs
if [[ ${selectDEAsoftware} == 'edger' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output="${workingDir}/${outDir}/${edgOut}/"
	dea_result_file="${path_output}/$(echo EdgeR_${contrast_label} | xargs).tsv"
	software="EdgeR"
	test_and_run "${tarbase_db_file}" "${reactome_db_file}" "${dea_result_file}" "${input_contrast}" "${software}" "${db_organism}" "${path_output}" 'EdgeR'
fi

if [[ ${selectDEAsoftware} == 'deseq' ]] || [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output="${workingDir}/${outDir}/${dsqOut}/"
	dea_result_file="${path_output}/$(echo DESeq2_${contrast_label} | xargs).tsv"
	software="DESeq2"
	test_and_run "${tarbase_db_file}" "${reactome_db_file}" "${dea_result_file}" "${input_contrast}" "${software}" "${db_organism}" "${path_output}" 'DESeq2'
fi

if [[ ${selectDEAsoftware} == 'both' ]]
then
	path_output="${workingDir}/${outDir}/${deaIntOut}/${contrast_label}/"
	dea_result_file="${path_output}/DEmiRNAs_${contrast_label}_deseq-edger_integrated.tsv"
	software='DESeq2-EdgeR'
	test_and_run "${tarbase_db_file}" "${reactome_db_file}" "${dea_result_file}" "${input_contrast}" "${software}" "${db_organism}" "${path_output}" 'integrated'
fi
