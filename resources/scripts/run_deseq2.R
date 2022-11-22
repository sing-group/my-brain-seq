## Script input parameters:
##  1.- path_counts: path to the counts file (without the first line).
##  2.- path_cond: path to the condition file.
##  3.- input_contrast: A line of the contrast file with one contrast.
##  4.- path_output: the output dir for the results.

suppressMessages(library("DESeq2"))

getLabel = function(sample_name) {
  #Uses conditions_table to get the label corresponding to a value in column$name
  index=match(sample_name, conditions_table$name)
  label = conditions_table$label[index]
  return(as.character(label))
}

#INPUTS
args = commandArgs(trailingOnly = TRUE)
path_counts=as.character(args[1])
path_cond=as.character(args[2])
input_contrast=as.character(args[3])
path_output=as.character(args[4])

#DETECT FACTORS AND LABELS
#read the input files
contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=')
conditions_table = read.delim(path_cond)
rownames(conditions_table) = conditions_table$name
#get the reference factor and the condition (removing trailing/leading whitespaces)
contrast = strsplit(as.character(contrast_table$name), '-')[[1]]
ref_factor = trimws(contrast[2])
cond_factor = as.character(trimws(contrast[1]))
print(paste0('[PIPELINE -- edger]: Reference factor: ', ref_factor))
print(paste0('[PIPELINE -- edger]: Condition factor: ', cond_factor))

#get the label of the reference factor
ref_contrast_label=trimws(strsplit(as.character(rownames(contrast_table)), '-')[[1]][2])
#get the label of the condition factor
cond_contrast_label=trimws(strsplit(as.character(rownames(contrast_table)), '-')[[1]][1])

#import the read counts and remove the featureCount annotations
cts = read.delim(path_counts, row.names="Geneid")
cts = cts[-1:-5]
#builds a list with the samples labeled in order
user_labels = c() #e.g.: C,C,PA,PA
sampleRepl = c() #e.g.: C01,C02,PA01,PA02
for (i in conditions_table$name){
  #find the index of each label in "cts" by using grep
  index = grep(i, colnames(cts))
  #change "cts" colnames to its corresponding label
  colnames(cts)[index] = i
  #store values on each iteration
  sampleRepl = c(sampleRepl, i)
  user_labels = c(user_labels, getLabel(i))
}

#PREPARE COFOUNDING FACTORS IF ANY
# get confounding variables
default_columns = c('name', 'condition', 'label')
adjust_for = subset(colnames(conditions_table),
                    !colnames(conditions_table) %in% default_columns)
print_adjust = paste(adjust_for, collapse=', ')
if (length(adjust_for > 1)){
  print(paste0('[PIPELINE -- edger]: Adjust model for: ', print_adjust))
}else{
  print(paste0('[PIPELINE -- edger]: No model adjustment'))
}

#confounding variables to factors
if (length(adjust_for > 1)){
  for (i in adjust_for) {
    conditions_table[[i]] = factor(conditions_table[[i]])
  }
}

#SUBSET THE DESIRED FACTORS
print('[PIPELINE -- edger]: Subsetting desired factors')
# samples with the desired factors to compare
des_samples = as.vector(subset(conditions_table, conditions_table$condition == ref_factor | conditions_table$condition == cond_factor )$name)
# remove the undesired samples from the counts
cts=cts[,(names(cts) %in% des_samples)]
# remove the undesired samples from the conditions_table
conditions_table=conditions_table[(conditions_table$name %in% des_samples),]


#################### END OF THE COMMON PART (DESEQ2 EDGER) #####################

# Checking that the cts and annotations are correct
test_col = all(rownames(conditions_table) %in% colnames(cts))
print(paste0("Rownames of 'condition_file.txt' are the same as colnames of 'all-counts.txt': ", test_col))

# If annotations are incorrect, show both headers
if (test_col == FALSE){
  print("[ERROR]: Rownames of 'condition_file.txt' are NOT the same as colnames of 'all-counts.txt'")
  print("    ROWNAMES of 'condition_file.txt': ")
  print(paste('         1. ', rownames(conditions_table)[1], '        '))
  print(paste('         2. ', rownames(conditions_table)[2], '        '))
  print(paste('         3. ', '...', '        '))
  print("    COLNAMES of 'all-counts.txt': ")
  print(paste('         1. ', colnames(cts)[1], '        '))
  print(paste('         2. ', colnames(cts)[2], '        '))
  print(paste('         3. ', '...', '        '))
}

#Test if colnames of cts are in the same order that rownames of conditions_table
if (all(rownames(conditions_table) == colnames(cts)) == FALSE) {
  cts = cts[, rownames(conditions_table)] #if not, arranges them
}

#-------------------------------------------------------------------------------
#                             DESEQ2 ANALYSIS
#-------------------------------------------------------------------------------
#DESEQ2 DESIGN
#builds the edger model matrix
if (length(adjust_for > 0)){
  adjustment = paste(adjust_for, collapse=' + ', sep = '')
  design = paste("~", adjustment, ' + condition', sep='')
} else {
  design = "~condition"
}

print(paste0('[PIPELINE -- edger]: Model: ', design))
print('')

#Construct DESeq2 dataset
dds = DESeqDataSetFromMatrix(countData = cts,
                             colData = conditions_table,
                             design = formula(design))

#Keep only rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Performs the differential expression analysis
dds <- DESeq(dds)

#Order results by the smallest p-value and visualize
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
resOrdered

################# BEGINNING OF THE COMMON PART (DESEQ2 EDGER) ##################

#SAVE RESULTS
output_tag = paste('DESeq2_', cond_contrast_label, '-', ref_contrast_label, sep = '')
output_file = paste(output_tag, '.tsv', sep = '')
path_output_file = paste(path_output, output_file, sep='')

# results = topTags(et, nrow(et))$table # Line exclusive of EdgeR analysis
dataframe_save = as.data.frame(resOrdered)
dataframe_save = cbind(Feature = rownames(resOrdered), dataframe_save)
write.table(dataframe_save,
            path_output_file,
            row.names = FALSE,
            col.names = TRUE, 
            sep = '\t')
