#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- dea-integration]: Integrating DESeq2 and EdgeR results..."

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${deaIntRscript} ${workingDir}/compi_scripts/${deaIntRscript}

dint_output="${workingDir}/${outDir}/${deaIntOut}"
dint_output_dk="${workingDir}/${outDir}/${deaIntOut}"

# find the filenames of DESeq2 and EdgeR results
contrast_filename=$(tail ${contrast} --lines=+2 | cut -d'=' -f1 | tr -d \")
vp_comparison_label="$(echo $contrast_filename | xargs)"
deseq_file="${workingDir}/${outDir}/${dsqOut}/DESeq2_${vp_comparison_label}.tsv"
edger_file="${workingDir}/${outDir}/${edgOut}/EdgeR_${vp_comparison_label}.tsv"

echo "[PIPELINE -- dea-integration]: Filtering results by q-value"
#remove header and filters the results by q-value
cat "${deseq_file}" | tail -n +2 | awk '$7<5e-2' | sort > "${dint_output}/deseq_qval-0.05.tsv"
cat "${edger_file}" | tail -n +2 | awk '$4<5e-2' | sort > "${dint_output}/edger_qval-0.05.tsv"

echo "[PIPELINE -- dea-integration]: Saving filtered results"
#build one file with both results: first column "Feature", second "log2FC" and third "q-value"
cat <(cut -f1,3,7 "${dint_output}/deseq_qval-0.05.tsv") \
	<(cut -f1,2,4 "${dint_output}/edger_qval-0.05.tsv") > "${dint_output}/all.tsv"

echo "[PIPELINE -- dea-integration]: Findind coincidences"
#find duplicated lines in "Feature" column and saves them in "coincidences.tsv"
cat "${dint_output}/edger_qval-0.05.tsv" "${dint_output}/deseq_qval-0.05.tsv" | cut -f 1 | sort | uniq -d > "${dint_output}/coincidences.tsv"

#if coincidences.tsv is not empty, perform the integration
if [ -s "${dint_output}/coincidences.tsv" ]; then
	#Search the duplicates in all.tsv and adds them to a file, this file will have two rows per miRNA (one for DESeq2 results, the other for EdgeR)
	grep -f "${dint_output}/coincidences.tsv" "${dint_output}/all.tsv" > "${dint_output}/DEmiRNAs_deseq-edger_integrated.tsv"
	sed -i '1s/^/Feature	log2FC	q-value\n/' "${dint_output}/DEmiRNAs_deseq-edger_integrated.tsv"

	#Build the results of the proces by averaging the DESeq2 and EdgeR values
	echo "[PIPELINE -- dea-integration]: Averaging values and building the results..."
	docker run --rm \
		-v ${workingDir}:${workingDir} \
		pegi3s/r_data-analysis \
			Rscript ${workingDir}/compi_scripts/${deaIntRscript} ${dint_output_dk}/DEmiRNAs_deseq-edger_integrated.tsv ${dint_output_dk}

else
	echo "[PIPELINE -- dea-integration]: [WARNING]: No coincidences between Edger and DESeq2 results for a q-value < 0.05"
	echo "[PIPELINE -- dea-integration]: [WARNING]: Venn and Volcano plots skipped"
fi

echo "[PIPELINE -- dea-integration]: Done!"
