#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- deseq]: Performing differential expression analysis with DESeq2..."

SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh

# lock Rscript before copying to avoid errors when parallel tasks are running
cp_and_lock ${deSeq2Rscript} 'deseq' ${scriptsDir}

#INPUTS
#common
dea_comparison="$(echo ""${comparison}"" | cut -d= -f1 | tr -d \" | xargs)"
dea_path_cond="${conditions}"
dea_path_contrast="${contrast}"
#deseq
dea_path_pipel="${workingDir}/${outDir}/${dsqOut}/pipel"
dea_path_counts="${dea_path_pipel}/all-counts_deseq_${dea_comparison}.txt"
dea_path_output="${workingDir}/${outDir}/${dsqOut}/"

mkdir -p "${dea_path_pipel}"
touch "${dea_path_counts}"
cat "${workingDir}/${outDir}/${ftqOut}/all-counts.txt" | tail -n +2 > "${dea_path_counts}"


#echo "[PIPELINE -- edger]: Running EdgeR analysis..."
echo "[PIPELINE -- deseq]: Running DESeq2 analysis..."
docker run --rm \
	-v ${workingDir}:${workingDir} \
	-v ${dea_path_contrast}:${dea_path_contrast} \
	-v ${dea_path_cond}:${dea_path_cond} \
	pegi3s/r_deseq2:${rDeseq2Version} \
		Rscript ${workingDir}/compi_scripts/${deSeq2Rscript} "${dea_path_counts}" "${dea_path_cond}" "${comparison}" "${dea_path_output}"
