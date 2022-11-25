#-------------------------------------------------------------------------------
#                                  INPUTS
#-------------------------------------------------------------------------------
## Script input parameters:
##  1.- path_dea_results: path to the DEA file, e.g.: 'DESeq2_X-Y.tsv'.
##  2.- path_counts: path to the counts file (without the first line).
##  3.- path_conditions: path to the condition file.
##  4.- input_contrast: A line of the 'contrast_file' with one contrast.
##  5.- path_output: the output dir for the results.
##  6.- software: the software that produced the tsv file (deseq/edger/both).

# to print bare text without "[1]"
pt = function(text){
  cat(text, sep='\n')
}
# to print messages
ptm = function(text){
  header = '[PIPELINE -- hclust -- run_make_hclust_tables.R]:'
  cat(paste(header, text), sep='\n')
}

#INPUTS
args = commandArgs(trailingOnly = TRUE)
path_dea_results=as.character(args[1])
path_counts=as.character(args[2])
path_conditions=as.character(args[3])
input_contrast=as.character(args[4])
path_output=as.character(args[5])
software=as.character(args[6])

#INPUT
args = commandArgs(trailingOnly = TRUE)
pt(''); pt('')
pt("======================================================")
pt('   [PIPELINE -- hclust]: run_make_hclust_table.R      ')
pt('......................................................')
pt(paste0("  DEA file:        ", as.character(args[1])))
pt(paste0("  Counts file:     ", as.character(args[2])))
pt(paste0("  Conditions file: ", as.character(args[3])))
pt(paste0("  Contrast file:   ", as.character(args[4])))
pt(paste0("  Output dir:      ", as.character(args[5])))
pt(paste0("  Software:        ", as.character(args[6])))
pt('......................................................')
pt(''); pt('')

#-------------------------------------------------------------------------------
#                                FUNCTIONS
#-------------------------------------------------------------------------------
ptm('Loading functions')
suppressMessages(library(dplyr))

read_tsv = function(file, skip = 0, rowname = NULL) {
  # Read a tsv file
  table = read.csv(file = file,
                   sep = "\t",
                   header = TRUE,
                   check.names = FALSE,
                   row.names = rowname,
                   skip = skip)}

getLabel = function(sample_name, conditions_table) {
  #Uses conditions_table to get the label corresponding to a value in column$name
  index=match(sample_name, conditions_table$name)
  label = conditions_table$label[index]
  return(as.character(label))
}

fix_cts_colnames = function(counts_table, conditions_table) {
  #Change the full path in the colnames of counts to the name of the sample:
  #  /home/x/working_dir/output/3_bowtie/trimmed_reduced_C001.fastq.gz.bam ->
  #  reduced_C001
  #builds a list with the samples labeled in order
  user_labels = c() #e.g.: C,C,PA,PA
  sampleRepl = c() #e.g.: C01,C02,PA01,PA02
  for (i in conditions_table$name){
    #find the index of each label in "counts_table" by using grep
    index = grep(i, colnames(counts_table))
    #change "counts_table" colnames to its corresponding label
    colnames(counts_table)[index] = i
    #store values on each iteration
    sampleRepl = c(sampleRepl, i)
    user_labels = c(user_labels, getLabel(i, conditions_table))
  }
  return(counts_table)
}

get_contrast_factors = function(input_contrast) {
  #Get the reference and condition factors from the input_contrast (pipeline)
  #Get the current contrast
  contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=', header = F, skip = 1)
  #get the reference factor and the condition (removing trailing/leading white spaces)
  contrast = strsplit((contrast_table$V2), '-')[[1]]
  label_contrast = strsplit((contrast_table$V1), '-')[[1]]
  ref_factor = trimws(contrast[2])
  cond_factor = as.character(trimws(contrast[1]))
  label_reference_factor = trimws(label_contrast[2])
  label_condition_factor = as.character(trimws(label_contrast[1]))
  result = list("reference" = ref_factor, 
                "condition" = cond_factor,
                "label_reference" = label_reference_factor,
                "label_condition" = label_condition_factor)
}

#-------------------------------------------------------------------------------
#                               LOAD TABLES
#-------------------------------------------------------------------------------
ptm('Loading files')
#DEA results
dea_results = read_tsv(path_dea_results)
#Counts
cts = read_tsv(path_counts, rowname="Geneid")
cts = cts[-1:-5]
#Conditions
conditions = read_tsv(path_conditions)

#-------------------------------------------------------------------------------
#                             GET MIRNA LIST
#-------------------------------------------------------------------------------
ptm('Getting miRNA profile')
# if integrated results, then miRNA profile = first column of dea_results
if (software == "DESeq2-EdgeR"){
  mirna_profile = as.list(dea_results$Feature)
# else filter the miRNAs by q_value and FC
}else{
  #filter by q-value and FDR
  q_value = 0.05
  fold_change = 0.5
  dea_results_filtered = filter(dea_results,
                                qvalue<q_value & (log2FC>fold_change | log2FC< -fold_change))
  
  #get the miRNA list
  mirna_profile = dea_results_filtered$Feature
}

#-------------------------------------------------------------------------------
#                          PREPARE HCLUST TABLE
#-------------------------------------------------------------------------------
ptm('Preparing hclust table')
#Change the colnames of cts from paths to the name of the samples
cts = fix_cts_colnames(cts, conditions)
#Get the expression of the miRNA profile
hclust_table = filter(cts, rownames(cts) %in% mirna_profile)

#-------------------------------------------------------------------------------
#                             SAVE RESULTS
#-------------------------------------------------------------------------------
ptm('Saving hclust results')
#Output file name considering the software
contrast_factors = get_contrast_factors(input_contrast)
output_file = paste('hclust_',
                    software,
                    '_',
                    contrast_factors$label_condition,
                    '-',
                    contrast_factors$label_reference,
                    '.tsv',
                    sep = '')
ptm(paste('Results file:', output_file))
#Full output path
path_output_file = paste(path_output, output_file, sep='')
#Save table
write.table(hclust_table,
            path_output_file,
            row.names = TRUE,
            col.names = TRUE,
            sep = '\t')

ptm('Done')