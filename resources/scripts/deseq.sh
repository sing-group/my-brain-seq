#!/bin/bash
set -o nounset
set -o errexit

echo "[PIPELINE -- deseq]: Performing differential expression analysis with DESeq2..."

#Makes a copy of the scripts used in the analysis to working-dir
cp ${scriptsDir}/${deSeq2Rscript} ${workingDir}/compi_scripts/${deSeq2Rscript}

export path_output="${workingDir}/${outDir}/${dsqOut}"
export path_output_docker="${workingDir}/${outDir}/${dsqOut}"
export path_scriptR_docker="${workingDir}/${scriptsDir}/${deSeq2Rscript}"

conversion_file="${path_output}/conversion_file.txt"

# get the reference factor for the comparison (eg.: control)
ref_factor=$(tail --lines=+2 ${contrast} | sed 's/".*" = //' | sed 's/".*-//' | sed 's/"$//')

# get the contrast name to build the output filename
contrast_filename=$(tail ${contrast} --lines=+2 | cut -d'=' -f1 | tr -d \")

# creates an empty file and adds the header of conversion_file.
head -n1 ${conditions} > ${conversion_file}

# converts the names of the conversion_file to a format recognized by the pipeline
echo "[PIPELINE -- deseq]: Converting condition_file format"
# if adapter specified, adds trimmed at the beginning of the rootname and ".fastq.bam" at the end
if [ "${adapter}" != "NA" ]; then
	tail --lines=+2 ${conditions} | sed 's/^/trimmed_/' | sed 's/\t/.fastq.bam\t/' >> ${conversion_file}
#if adapter not specified only adds ".fastq.bam" at the end
else
	tail --lines=+2 ${conditions} | sed 's/\t/.fastq.bam\t/' >> ${conversion_file}
fi

#Converts the "conditions" file to absolute paths
echo "[PIPELINE -- deseq]: Adding full paths to the condition file"
export dqCond="${path_output}/deseq_conditions_file.txt"
export dqCond_docker="${path_output_docker}/deseq_conditions_file.txt"

# Adds "${workingDir}/..." at the beginning of the row names
cat "${conversion_file}" | sed '2,$s+^+'"${workingDir}/${outDir}/${bwtOut}/"'+' | tr -s '/' > "${dqCond}"

echo "[PIPELINE -- deseq]: Contrast to perform: ${contrast_filename}"

echo "[PIPELINE -- deseq]: Running DESeq2 analysis..."
docker run --rm \
	-v ${workingDir}:${workingDir} \
	pegi3s/r_deseq2:${rDeseq2Version} \
		Rscript ${workingDir}/compi_scripts/${deSeq2Rscript} ${workingDir}/${outDir}/${ftqOut}/all-counts.txt ${dqCond_docker} ${ref_factor} ${path_output_docker} ${contrast_filename}
