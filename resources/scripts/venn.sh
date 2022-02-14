#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- venn]: Building a file for the Venn diagram"

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${vennRscript} ${workingDir}/compi_scripts/${vennRscript}

#Output dir
venn_output_local="${workingDir}/${outDir}/${deaIntOut}/"

# find the filenames of filtered DESeq2 and EdgeR results
echo "[PIPELINE -- venn]: Finding filtered DESeq2 and EdgeR results"
contrast_filename=$(tail ${contrast} --lines=+2 | cut -d'=' -f1 | tr -d \")
vp_comparison_label="$(echo $contrast_filename | xargs)"
venn_des005_path=$(echo "${venn_output_local}/deseq_${vp_comparison_label}_qval-0.05.tsv")
venn_edg005_path=$(echo "${venn_output_local}/edger_${vp_comparison_label}_qval-0.05.tsv")
echo "[PIPELINE -- venn]:	DESeq2 file: $venn_des005_path"
echo "[PIPELINE -- venn]:	EdgeR  file: $venn_edg005_path"


#Build a file with the DE features found by DESeq2 (col 1) and EdgeR (col 2)
join_file="DEmiRNAs_deseq-edger_${vp_comparison_label}_venn.tsv"
echo "[PIPELINE -- venn]: Joining DESeq2 and EdgeR results in one file"
paste <(cut -f1 ${venn_des005_path}) <(cut -f1 ${venn_edg005_path}) > ${venn_output_local}/${join_file}
echo "[PIPELINE -- venn]:	Done, saved in: ${venn_output_local}/${join_file}"

venn_path="${workingDir}/${outDir}/${deaIntOut}/${join_file}"
venn_output="${workingDir}/${outDir}/${deaIntOut}/"
venn_output_format=${vennFormat}

echo "[PIPELINE -- venn]: Building the Venn diagram"
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_venn-diagram \
		Rscript ${workingDir}/compi_scripts/${vennRscript} ${venn_path} ${venn_output} ${venn_output_format} ${vp_comparison_label}
