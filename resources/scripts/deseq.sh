#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- deseq]: Performing differential expression analysis with DESeq2..."

#Makes a copy of the scripts used in the analysis to working-dir
cp "${scriptsDir}/${deSeq2Rscript}" "${workingDir}/compi_scripts/${deSeq2Rscript}"
cp "${scriptsDir}/${filterCtsRscript}" "${workingDir}/compi_scripts/${filterCtsRscript}"

#--------------------------------------------
# IDENTIFY THE FACTORS OF THE CONTRAST
#--------------------------------------------
# get the reference factor for the comparison (eg.: control)
ref_factor=$(echo "${comparison}" | cut -d '=' -f2 | xargs | cut -d'-' -f2)
# get the condition name to adapt the conditions_file.txt
cond_factor=$(echo "${comparison}" | cut -d '=' -f2 | xargs | cut -d'-' -f1)
# get the contrast name to build the output filename
contrast_filename=$(echo "${comparison}" | head -n1 | cut -d'=' -f1 | xargs | tr -d \")

#--------------------------------------------
# PATHS
#--------------------------------------------
export path_output_results="${workingDir}/${outDir}/${dsqOut}"
export path_output_pipel="${workingDir}/${outDir}/${dsqOut}/pipel"
export path_scriptR_docker="${workingDir}/${scriptsDir}/${deSeq2Rscript}"
conversion_file="${path_output_pipel}/conversion_file_${contrast_filename}.txt"

# creates the directory for the pipeline files
mkdir -p ${path_output_pipel}

#--------------------------------------------
# ADDING FULL PATHS TO FILTERERD CONDITIONS
#--------------------------------------------
echo "[PIPELINE -- deseq]: Converting $(basename ${conditions}) format to full paths"
# creates an empty file and adds the header of conditions_file.
head -n1 ${conditions} > "${conversion_file}"
# converts the names of the conversion_file to a format recognized by the pipeline
# if adapter specified, adds trimmed at the beginning of the rootname and ".fastq.bam" at the end
if [ "${adapter}" != "NA" ]; then
	tail --lines=+2 ${conditions} | sed 's/^/trimmed_/' | sed 's/\t/.fastq.bam\t/' >> "${conversion_file}"
#if adapter not specified only adds ".fastq.bam" at the end
else
	tail --lines=+2 ${conditions} | sed 's/\t/.fastq.bam\t/' >> "${conversion_file}"
fi

#Converts the "conditions" file to absolute paths
echo "[PIPELINE -- deseq]: Adding full paths to the condition file"
export dqCond="${path_output_pipel}/deseq_conditions_file.txt"
export dqCond_docker="${path_output_pipel}/deseq_conditions_file.txt"

# Adds "${workingDir}/..." at the beginning of the row names
cat "${conversion_file}" | sed '2,$s+^+'"${workingDir}/${outDir}/${bwtOut}/"'+' | tr -s '/' > "${dqCond}"

#--------------------------------------------
#        FILTERING ALL_COUNTS.TXT
#--------------------------------------------
echo '[PIPELINE -- deseq]: Filtering conditions based on: "${comparison}"'
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_deseq2:${rDeseq2Version} \
		Rscript "${workingDir}/compi_scripts/${filterCtsRscript}" "${workingDir}/${outDir}/${ftqOut}/all-counts.txt" "${dqCond}" "${comparison}" "${path_output_pipel}"
		
counts_filtered="${path_output_pipel}/counts_${cond_factor}-${ref_factor}.tsv"
conditions_filtered="${path_output_pipel}/conditions_${cond_factor}-${ref_factor}.tsv"

#--------------------------------------------
#       RUNNING THE DESEQ ANALYSIS
#--------------------------------------------
echo "[PIPELINE -- deseq]: Contrast to perform: ${contrast_filename}"
echo "[PIPELINE -- deseq]: Running DESeq2 analysis..."
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_deseq2:${rDeseq2Version} \
		Rscript "${workingDir}/compi_scripts/${deSeq2Rscript}" "${counts_filtered}" "${conditions_filtered}" "${ref_factor}" "${path_output_results}" "${contrast_filename}"
