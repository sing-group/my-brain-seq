#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- venn]: Building a file for the Venn diagram"

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${vennRscript} ${workingDir}/compi_scripts/${vennRscript}

#Build a file with the DE features found by DESeq2 (col 1) and EdgeR (col 2)
venn_output_local="${workingDir}/${outDir}/${deaIntOut}/"
venn_des005_path=$(echo "${venn_output_local}/deseq_qval-0.05.tsv")
venn_edg005_path=$(echo "${venn_output_local}/edger_qval-0.05.tsv")
paste <(cut -f1 ${venn_des005_path}) <(cut -f1 ${venn_edg005_path}) > ${venn_output_local}/DEmiRNAs_deseq-edger_venn.tsv

venn_path="${workingDir}/${outDir}/${deaIntOut}/DEmiRNAs_deseq-edger_venn.tsv"
venn_output="${workingDir}/${outDir}/${deaIntOut}/"
venn_output_format=${vennFormat}

echo "[PIPELINE -- venn]: Building the Venn diagram"
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_venn-diagram \
		Rscript ${workingDir}/compi_scripts/${vennRscript} ${venn_path} ${venn_output} ${venn_output_format}
