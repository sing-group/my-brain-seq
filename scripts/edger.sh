#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- edger]: Performing differential expression analysis with EdgeR..."

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${edgerRscript} ${workingDir}/compi_scripts/${edgerRscript}

touch "${workingDir}/${outDir}/${edgOut}/all-counts_edger.txt"
cat "${workingDir}/${outDir}/${ftqOut}/all-counts.txt" | tail -n +2 > "${workingDir}/${outDir}/${edgOut}/all-counts_edger.txt"

#Inputs
er_path_counts="${workingDir}/${outDir}/${edgOut}/all-counts_edger.txt"
er_path_cond=${conditions}
er_path_contrast=${contrast}
er_path_output="${workingDir}/${outDir}/${edgOut}/"

echo "[PIPELINE -- edger]: Running EdgeR analysis..."
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_edger:${rEdgerVersion} \
		Rscript ${workingDir}/compi_scripts/${edgerRscript} ${er_path_counts} ${er_path_cond} ${er_path_contrast} ${er_path_output}