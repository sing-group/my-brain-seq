#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- edger]: Performing differential expression analysis with EdgeR..."

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${edgerRscript} ${workingDir}/compi_scripts/${edgerRscript}

#Inputs
er_comparison="$(echo ""${comparison}"" | cut -d= -f1 | tr -d \" | xargs)"
er_path_pipel="${workingDir}/${outDir}/${edgOut}/pipel"
er_path_counts="${er_path_pipel}/all-counts_edger_${er_comparison}.txt"
er_path_cond="${conditions}"
er_path_contrast="${contrast}"
er_path_output="${workingDir}/${outDir}/${edgOut}/"

mkdir -p "${er_path_pipel}"
touch "${er_path_counts}"
cat "${workingDir}/${outDir}/${ftqOut}/all-counts.txt" | tail -n +2 > "${er_path_counts}"


echo "[PIPELINE -- edger]: Running EdgeR analysis..."
docker run --rm \
	-v ${workingDir}:${workingDir} \
	-v ${er_path_contrast}:${er_path_contrast} \
	-v ${er_path_cond}:${er_path_cond} \
	pegi3s/r_edger:${rEdgerVersion} \
		Rscript ${workingDir}/compi_scripts/${edgerRscript} "${er_path_counts}" "${er_path_cond}" "${comparison}" "${er_path_output}"
