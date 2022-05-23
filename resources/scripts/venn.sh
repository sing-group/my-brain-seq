#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- venn]: Building a file for the Venn diagram"

# makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${vennRscript} ${workingDir}/compi_scripts/${vennRscript}

# find the filenames of filtered DESeq2 and EdgeR results
echo "[PIPELINE -- venn]: Finding filtered DESeq2 and EdgeR results"
vp_comparison_label=$(echo "${comparison}" | cut -d'"' -f2)
#vp_comparison_label="$(echo $contrast_filename | xargs)"

# output dir
pipel_dir_name='pipel'
output="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}/"
output_pipel="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}/${pipel_dir_name}/"

venn_des005_path=$(echo "${output_pipel}/deseq_${vp_comparison_label}_qval-0.05.tsv")
venn_edg005_path=$(echo "${output_pipel}/edger_${vp_comparison_label}_qval-0.05.tsv")
echo "[PIPELINE -- venn]:	DESeq2 file: $venn_des005_path"
echo "[PIPELINE -- venn]:	EdgeR  file: $venn_edg005_path"


# build a file with the DE features found by DESeq2 (col 1) and EdgeR (col 2)
join_file="DEmiRNAs_deseq-edger_${vp_comparison_label}_venn.tsv"
echo "[PIPELINE -- venn]: Joining DESeq2 and EdgeR results in one file"
paste <(cut -f1 "${venn_des005_path}") <(cut -f1 "${venn_edg005_path}") > "${output_pipel}/${join_file}"
echo "[PIPELINE -- venn]:	Done, saved in: ${output_pipel}/${join_file}"

venn_path="${output_pipel}${join_file}"
venn_output="${output}"
venn_output_format="${vennFormat}"

echo "[PIPELINE -- venn]: Building the Venn diagram"
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_venn-diagram \
		Rscript "${workingDir}/compi_scripts/${vennRscript}" "${venn_path}" "${venn_output}" "${venn_output_format}" "${vp_comparison_label}"
