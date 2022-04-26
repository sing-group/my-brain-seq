# myBrain-Seq [![license](https://img.shields.io/badge/license-MIT-brightgreen)](https://github.com/sing-group/my-brain-seq) [![dockerhub](https://img.shields.io/badge/hub-docker-blue)](https://hub.docker.com/r/singgroup/my-brain-seq) [![compihub](https://img.shields.io/badge/hub-compi-blue)](https://www.sing-group.org/compihub/explore/625e719acc1507001943ab7f)
> **myBrain-Seq** is a [Compi](https://www.sing-group.org/compi/) pipeline for miRNA-Seq analysis of neuropsychiatric data. A Docker image is available for this pipeline in [this Docker Hub repository](https://hub.docker.com/r/singgroup/my-brain-seq).

## myBrain-Seq repositories

- [GitHub](https://github.com/sing-group/my-brain-seq)
- [Docker Hub](https://hub.docker.com/r/singgroup/my-brain-seq)
- [CompiHub](https://www.sing-group.org/compihub/explore/625e719acc1507001943ab7f)

# What does myBrain-Seq do?

**MyBrain-Seq** is a [Compi](https://www.sing-group.org/compi/) pipeline to automatically identify differentially expressed microRNAs (miRNAs) between two conditions. It combines two differential expression analysis software, namely DESeq2 and EdgeR, to offer an integrated results suitable for experimental validation. Its features and analysis are designed and tuned to improve its performance with neuropsychiatric data, but it can also be used with any other type of miRNA data.

MyBrain-Seq preprocess each input FastQ file separately before the differential expression analysis. This process comprises:

1. Quality control of the sequences using FastQC.
2. Trimming of the adapter sequences using Cutadapt (optional).
3. Alignment to the reference genome with Bowtie. 
4. Format conversion of the Bowtie output files to BAM using sam-tools.
5. Quality control of the alignments with sam-tools. 
6. Quantification and annotation with featureCounts.

After all this steps were completed for each sample, myBrain-Seq performs the differential expression analysis. This process comprises:

1. Differential expression analysis with DESeq2. 
2. Differential expression analysis with EdgeR.
3. Intersection of the DESeq2 and EdgeR results and averagement of their q-values and fold change (optional).
4. Creation of a venn diagram with the integrated results using VennDiagram.
5. Creation of a volcano plot with the results using EnhancedVolcano.

In addition, the user can instruct myBrain-Seq to generate a genome index for the Bowtie alignment; this index will be built in parallel with the initialization tasks.



# Using the myBrain-Seq image in Linux
## Building the directory tree

To start a new analysis, the first thing to do is build the directory tree in your local file system (`mbs_project` in the example). This directory tree will be refered as the "working directory" and its structure is recognized and used by the pipeline during the analysis. It has the following structure: 

```
/path/to/mbs-project 
	|-- input 
	|   |-- compi.parameters 
	|   |-- conditions_file.txt 
	|   `-- contrast_file.txt 
	`-- output
```

Where:

- **compi.parameters** contains the paths and parameters needed for the analysis.
- **conditions_file.txt** contains the names and conditions of each fastQ file. 
- **contrast_file.txt** contains the names and labels of the conditions to compare in the differential expression analysis.

 The creations of these files is detailed in the following sections. You may find it convenient to create additional directories and files within the working directory to group all the data related to a particular study.

## Writting the `compi.parameters` file

The `compi.parameters` file is used by myBrain-Seq to locate the files needed for the analysis as well as to define which optional tasks will be run. Here is an example of a `compi.parameters` file using the working directory created in the previous example:

```
workingDir=/path/to/mbs-project
fastqDir=/path/to/study_1/data/
gffFile=/path/to/study_1/refs/mirbase_hsa.gff3
conditions=/path/to/mbs-project/input/conditions_file_study_1.txt
contrast=/path/to/mbs-project/input/contrast_file_study_1.txt
bwtIndex=/path/to/study_1/refs/bowtie-index_GRCh38
adapter=TGGAATTCTCGGGTGCCAAGG
```

This file contains the following mandatory parameters:

- **workingDir**: the path to the myBrain-Seq working directory of the analysis (first example of this section).
- **fastqDir**: the path to the directory with the fastQ files to analyse.
- **gffFile**: the path to the the GFF3 file with the miRNA annotations. This file could be obtained from [miRBase](https://www.mirbase.org/ftp.shtml) or [NCBI Genomes](https://www.ncbi.nlm.nih.gov/genome/).
- **conditions**: the path to conditions_file.txt.
- **contrast**: the path to contrast_file.txt.
- **genome** *(optional if bwtIndex is provided)*: the path to the reference genome in FASTA from which the Bowtie index will be built.
- **bwtIndex** *(optional if genome is provided)*: the path to a directory containing a Bowtie index. If this parameter is omitted myBrain-Seq will build a new index using a genome in FASTA provided in the genome parameter.

And the following optional parameters:

- **adapter**: the sequence of the adapter to remove. If this parameter is omitted myBrain-Seq will skip the adapter removal step.
- **gffFeature**: the name of the feature of the GFF3 file from which the attributes will be obtained; the default value is "miRNA".
- **gffAttribute**: the name of the attribute in the GFF3 file to use in the annotations; the default value is "Name".

## Writting the `conditions_file.txt` file

The `conditions_file.txt` is a TSV file used by myBrain-Seq to link each fastQ file with its condition. This information will be used to choose the group of samples to compare in the differential expression analysis. Here is an example of a conditions file:

	name	condition	label
	C019 	control		C_before_treatment
	C020 	control		C_before_treatment
	C021 	control		C_after_treatment
	C022	control		C_after_treatment
	P012D	FE 			FE_before_treatment
	P013A	FE			FE_before_treatment
	P014A	FE			FE_after_treatment
	P015D	FE			FE_after_treatment
	P014A	SEP			SEP_before_treatment
	P015D	SEP			SEP_before_treatment

In order to obtain a file with a valid format, the following considerations must be taken into account:

- Columns must be separated by single tabulations.
- The first row must be the header: “name”, “condition” and “label”.
- The first row must be the file rootnames of the fastQ files, i.e.: C019.fastq --> C019.
- The second row must be the conditions.
 - The third column is the label, which is only used so that the user can identify each sample in case there is more than one condition. It has no impact on the analysis result and can be omitted.

## Writting the `contrast_file.txt` file

The `contrast_file.txt` is used by myBrain-Seq to perform the comparisons between samples of two different conditions in the differential expression analysis. Each line on this file corresponds with a contrast that myBrain-Seq has to perform. Here is an example of a contrast file:

```
name
"Control_vs_first_episode" = "control-FE"
"Control_vs_second_episode" = "control-SEP"
"First_episode_vs_second_episode" = "FE-SEP"
```

In order to obtain a file with a valid format, the following considerations must be taken into account:

- The first row should be "name", in lowercase.
- In the second row, the contrast label must appear double quoted, then a space, an equal, another space, and between double quotes the contrast: factor to compare - reference factor (the reference factor is usually the control) .
- The factors to be compared must be the same as those specified in the "condition" column of the conditions_file.

## Running the myBrain-Seq analysis

Once all the required files were built you need to **build a runner** in order to perform the myBrain-Seq analysis. This can be done using the script `make_run-sh.sh`. This script uses the `compi.parameters` file as reference and will mount all the needed Docker volumes and build a directory for the logs. To use this script **adapt the first line** on the following code:

```bash
COMPI_PARAMETERS=/path/to/compi.parameters

WD=$(cat $COMPI_PARAMETERS | grep 'workingDir' | cut -d'=' -f2)
docker run --rm --entrypoint /init-working-dir/make_run-sh.sh -v ${WD}:${WD} -v ${COMPI_PARAMETERS}:${WD}/compi.parameters singgroup/my-brain-seq ${WD}/compi.parameters
```

The resulting file will be an script saved in the working directory and named as `run_<name of the workingDir>.sh`. To begin with the analysis just run it, for example:

```
./run_mbs-project.sh
```

## Find out tasks with errors

Some tasks may produce errors that do not cause the pipeline to fail, but they can be important. Such errors are reported in the log files produced in the `logs` directory of the pipeline working directory. Inside this directory myBrain-Seq will create additional directories with the logs of each execution, they will be named with the date and hour of the analysis. Files containing the errors are saved with extension ".err.log", whereas normal output is saved with extension ".out.log".



# Pipeline parameters

These are the pipeline parameters:

- `workingDir`: The working directory of the project.

- `fastqDir`: The directory containing the fastq files (default is relative to workingDir).

- `outDir`: The directory containing the pipeline outputs (relative to workingDir).

- `adapter`: The sequence of the adapter to remove.

- `genome`: The directory path to the genome to align.

- `bwtIndex`: The directory path containing the Bowtie index.

- `gffFile`: The path to the .gff file of the reference genome.

- `gffFeature`: Feature of the .gff file to use for the annotations (eg.: miRNA, gene, transcript...), default miRNA.

- `gffAttribute`: Attribute of the .gff to use in the annotations (eg.: Name, gene_id, transcript_id...), default Name.

- `conditions`: The path to the .tsv file with the rootnames of the samples, conditions and labels.

- `contrast`: The path to the .tsv file with the contrast DESeq2 has to perform.

- `vennFormat`: The file format of the Venn diagram (png/svg/tiff), default png.

- `fqcOut`: The relative path to the directory containing the FastQC results.

- `ctdOut`: The relative path to the directory containing the Cutadapt results.

- `bwtOut`: The relative path to the directory containing the Bowtie results.

- `bamstOut`: The relative path to the directory containing the Samtools stats and Plot-bamstats results.

- `ftqOut`: The relative path to the directory containing the FeatureCounts results.

- `dsqOut`: The relative path to the directory containing the DESeq2 results.

- `edgOut`: The relative path to the directory containing the EdgeR results.

- `deaIntOut`: The relative path to the directory containing the results of the DESeq2 and EdgeR integration.

- `scriptsDir`: The relative path to the directory containing the R script to run the DESeq2 analysis.

- `testAdapterBashScript`: The relative path to the directory containing the R script to get the path of the aligned/unaligned data.

- `deSeq2Rscript`: The relative path to the directory containing the R script to run the DESeq2 analysis.

- `filterCtsRscript`: The relative path to the directory containing the R script used to filter all-counts.txt and conditions_file.

- `edgerRscript`: The relative path to the directory containing the R script to run the EdgeR analysis.

- `enhancedVolcanoRscript`: The relative path to the directory containing the R script to build the Volcano plot.

- `deaIntRscript`: The relative path to the directory containing the R script to run the DESeq-EdgeR results integration.

- `vennRscript`: The relative path to the directory containing the R script to run the DESeq-EdgeR results integration.

- `rDeseq2Version`: Version of the pegi3s/r_deseq2 Docker image to use.

- `rEdgerVersion`: Version of the pegi3s/r_edger Docker image to use.

- `rEnhancedVolcanoVersion`: Version of the pegi3s/r_enhanced-volcano Docker image to use.		

- `cutadaptVersion`: Version of the pegi3s/cutadapt Docker image to use.

- `fastqcVersion`: Version of the pegi3s/fastqc Docker image to use.

- `bowtieVersion`: Version of the pegi3s/bowtie1 Docker image to use.

- `featureCountsVersion`: Version of the pegi3s/feature-counts Docker image to use.

- `samtoolsVersion`: Version of the pegi3s/samtools_bcftools Docker image to use.

- `samtoolsBamstatsVersion`: Version of the pegi3s/samtools_bcftools Docker image to use for bam analysis.

- `rdatanalysisVersion`: Version of the pegi3s/r_data-analysis Docker image to use.

- `rVennVersion`: Version of the pegi3s/r_venn-diagram Docker image to use.

- `skipPullDockerImages`: Use this flag to skip the pull-docker-images task.

- `selectDEAsoftware`: Use this param to select the differential expression analysis software (deseq, edger or both).

# Test data

The sample data is available [here](http://static.sing-group.org/software/myBrainSeq/downloads/test-data-0.1.0.zip). Download and uncompress it, you will get a directory named `working-dir` that contains an example of a functional working directory, were the data and biological references were grouped within it. Here you can find:

- A directory called `input`, with the `compi.parameters`, `condition_file.txt` and `contrast_file.txt` of this particular study.
- A directory called `data`, with the fastQ files of the study, the Bowtie index and the miRNA annotations. 

To run the pipeline with this test data, edit the `compi.parameters` (at `/working-dir/input`) and modify the paths to adapt them to the absolute location of the working directory in your computer (e.g.: `workingDir=/working-dir` could be `workingDir=/home/user/working-dir`). After doing this, just run the `run.sh` script included.

## Running time

- ≈ 5 minutes - 5 parallel tasks - Ubuntu 20.04.4 LTS, 8 CPUs (Intel® Core™ i7-9700 @ 3.00GHz), 16GB of RAM and SSD disk.

# For Developers

## Building the Docker image

To build the Docker image, [`compi-dk`](https://www.sing-group.org/compi/#downloads) is required. Once you have it installed, simply run `compi-dk build` from the project directory to build the Docker image. The image will be created with the name specified in the `compi.project` file. This file also specifies the version of compi that goes into the Docker image.


# Team 

MyBrain-Seq is a pipeline developed by the SING and NeuroEpigenetics Lab groups.

- Daniel Pérez-Rodríguez [![ORCID](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-8110-3567), daniel.prz.rodriguez@gmail.com 
- Hugo López-Fernández [![ORCID](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-6476-7206), hlfernandez@uvigo.es
- Roberto C. Agís-Balboa [![ORCID](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0001-9899-9569), roberto.carlos.agis.balboa@hotmail.com
