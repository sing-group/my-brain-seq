#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- dea-integration]: Integrating DESeq2 and EdgeR results..."

# makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${deaIntRscript} ${workingDir}/compi_scripts/${deaIntRscript}

# find the filenames of DESeq2 and EdgeR results
echo "[PIPELINE -- dea-integration]: Finding DESeq2 and EdgeR results"
contrast_filename=$(echo "${comparison}" | cut -d'=' -f1 | tr -d \")
vp_comparison_label="$(echo ""${contrast_filename}"" | xargs)"
deseq_file="${workingDir}/${outDir}/${dsqOut}/DESeq2_${vp_comparison_label}.tsv"
edger_file="${workingDir}/${outDir}/${edgOut}/EdgeR_${vp_comparison_label}.tsv"
echo "[PIPELINE -- dea-integration]:	DESeq2 file: $deseq_file"
echo "[PIPELINE -- dea-integration]:	EdgeR  file: $edger_file"

# output_pipel directories
pipel_dir_name='pipel'
output="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}"
output_pipel="${workingDir}/${outDir}/${deaIntOut}/${vp_comparison_label}/${pipel_dir_name}"

# make a directory for the pipeline intermediate files
mkdir -p ${output_pipel}

# remove header and filters the results by q-value
echo '[PIPELINE -- dea-integration]: Filtering results by q-value < 0.05'
echo "[PIPELINE -- dea-integration]:	Filtering DESeq2 results"
cat "${deseq_file}" | tail -n +2 | awk '$7<5e-2' | sort > "${output_pipel}/deseq_${vp_comparison_label}_qval-0.05.tsv"
dnum=$(cat "${deseq_file}" | tail -n +2 | awk '$7<5e-2' | wc -l)
echo "[PIPELINE -- dea-integration]:		Done, found $dnum coincidences"
echo "[PIPELINE -- dea-integration]:	Filtering EdgeR results"
cat "${edger_file}" | tail -n +2 | awk '$5<5e-2' | sort > "${output_pipel}/edger_${vp_comparison_label}_qval-0.05.tsv"
enum=$(cat "${edger_file}" | tail -n +2 | awk '$5<5e-2' | wc -l)
echo "[PIPELINE -- dea-integration]:		Done, found $enum coincidences"

# build one file with both results: first column "Feature", second "log2FC" and third "q-value"
echo "[PIPELINE -- dea-integration]: Saving filtered results"
cat <(cut -f1,3,7 "${output_pipel}/deseq_${vp_comparison_label}_qval-0.05.tsv") \
	<(cut -f1,2,4 "${output_pipel}/edger_${vp_comparison_label}_qval-0.05.tsv") > "${output_pipel}/all_${vp_comparison_label}.tsv"
echo "[PIPELINE -- dea-integration]:	Saved in: ${output_pipel}/all_${vp_comparison_label}.tsv"

# find duplicated lines in "Feature" column and saves them in "coincidences.tsv"
echo "[PIPELINE -- dea-integration]: Findind coincidences"
cat "${output_pipel}/edger_${vp_comparison_label}_qval-0.05.tsv" "${output_pipel}/deseq_${vp_comparison_label}_qval-0.05.tsv" | cut -f 1 | sort | uniq -d > "${output_pipel}/coincidences_${vp_comparison_label}.tsv"

# if coincidences.tsv is not empty, perform the integration
if [ -s "${output_pipel}/coincidences_${vp_comparison_label}.tsv" ]; then
	# Search the duplicates in all.tsv and adds them to a file, this file will have two rows per miRNA (one for DESeq2 results, the other for EdgeR)
	grep -f "${output_pipel}/coincidences_${vp_comparison_label}.tsv" "${output_pipel}/all_${vp_comparison_label}.tsv" > "${output_pipel}/DEmiRNAs_${vp_comparison_label}_deseq-edger_integrated.tsv"
	sed -i '1s/^/Feature	log2FC	q-value\n/' "${output_pipel}/DEmiRNAs_${vp_comparison_label}_deseq-edger_integrated.tsv"

	# Build the results of the proces by averaging the DESeq2 and EdgeR values
	echo "[PIPELINE -- dea-integration]: Averaging values and building the results..."
	docker run --rm \
		-v ${workingDir}:${workingDir} \
		pegi3s/r_data-analysis \
			Rscript ${workingDir}/compi_scripts/${deaIntRscript} ${output_pipel}/DEmiRNAs_${vp_comparison_label}_deseq-edger_integrated.tsv ${output} ${vp_comparison_label}

else
	echo "[PIPELINE -- dea-integration]:	[WARNING]: No coincidences between EdgeR and DESeq2 results for a q-value < 0.05"
	echo "[PIPELINE -- dea-integration]:	[WARNING]: Venn and Volcano plots skipped"
fi

echo "[PIPELINE -- dea-integration]: Done!"
