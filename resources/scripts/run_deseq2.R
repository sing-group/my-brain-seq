## Script input parameters:
##  1.- path_counts: path to the counts file (without the first line).
##  2.- path_cond: path to the condition file.
##  3.- input_contrast: A line of the contrast file with one contrast.
##  4.- path_output: the output dir for the results.

suppressMessages(library("DESeq2"))

# to print without "[1]"
pt = function(text){
  cat(text, sep='\n')
}
# to print messages
ptm = function(text){
  header = '[PIPELINE -- deseq -- run_deseq2.R]:'
  cat(paste(header, text), sep='\n')
}

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

padj = 0.05
logFC = 0.5

ptm("======================================================")
ptm('        [PIPELINE -- deseq]: run_deseq2.R        ')
ptm('......................................................')
ptm(paste0("  Counts file:     ", as.character(args[1])))
ptm(paste0("  Condition file:  ", as.character(args[2])))
ptm(paste0("  Input contrast:  ", as.character(args[3])))
ptm(paste0("  Output path:     ", as.character(args[4])))
ptm("======================================================")


#DETECT FACTORS AND LABELS
#read the input files
contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=')
conditions_table = read.delim(path_cond)
rownames(conditions_table) = conditions_table$name
#get the reference factor and the condition (removing trailing/leading whitespaces)
contrast = strsplit(as.character(contrast_table$name), '-')[[1]]
ref_factor = trimws(contrast[2])
cond_factor = as.character(trimws(contrast[1]))
ptm(paste0('Reference factor: ', ref_factor))
ptm(paste0('Condition factor: ', cond_factor))

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
  ptm(paste0('Adjust model for: ', print_adjust))
}else{
  ptm(paste0('No model adjustment'))
}

#confounding variables to factors
if (length(adjust_for > 1)){
  for (i in adjust_for) {
    conditions_table[[i]] = factor(conditions_table[[i]])
  }
}

#SUBSET THE DESIRED FACTORS
ptm('Subsetting desired factors')
# samples with the desired factors to compare
des_samples = as.vector(subset(conditions_table, conditions_table$condition == ref_factor | conditions_table$condition == cond_factor )$name)
# remove the undesired samples from the counts
cts=cts[,(names(cts) %in% des_samples)]
# remove the undesired samples from the conditions_table
conditions_table=conditions_table[(conditions_table$name %in% des_samples),]


#################### END OF THE COMMON PART (DESEQ2 EDGER) #####################

# Checking that the cts and annotations are correct
test_col = all(rownames(conditions_table) %in% colnames(cts))
ptm(paste0("Rownames of 'condition_file.txt' are the same as colnames of 'all-counts.txt': ", test_col))

# If annotations are incorrect, show both headers
if (test_col == FALSE){
  ptm("[ERROR]: Rownames of 'condition_file.txt' are NOT the same as colnames of 'all-counts.txt'")
  pt("    ROWNAMES of 'condition_file.txt': ")
  pt(paste('         1. ', rownames(conditions_table)[1], '        '))
  pt(paste('         2. ', rownames(conditions_table)[2], '        '))
  pt(paste('         3. ', '...', '        '))
  pt("    COLNAMES of 'all-counts.txt': ")
  pt(paste('         1. ', colnames(cts)[1], '        '))
  pt(paste('         2. ', colnames(cts)[2], '        '))
  pt(paste('         3. ', '...', '        '))
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

ptm(paste0('Model: ', design))
pt('')

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

#Get a list with the DE miRNAs
DE_miRNAs = res[!is.na(res$padj),]
DE_miRNAs = DE_miRNAs[DE_miRNAs$padj <= padj,]
DE_miRNAs = DE_miRNAs[DE_miRNAs$log2FoldChange <= -logFC | DE_miRNAs$log2FoldChange >= logFC,]

################# BEGINNING OF THE COMMON PART (DESEQ2 EDGER) ##################

#SAVE RESULTS
pt('')
ptm('Saving results')
output_tag = paste('DESeq2_', cond_contrast_label, '-', ref_contrast_label, sep = '')
output_file = paste(output_tag, '.tsv', sep = '')
path_output_file = paste(path_output, output_file, sep='')

# save DE miRNA list
output_file_DEmiRNAs = paste('differentially-expressed-miRNAs_', output_file, sep = '')
path_output_file_DEmiRNAs = paste(path_output, output_file_DEmiRNAs, sep='')
write.table(rownames(DE_miRNAs),
            path_output_file_DEmiRNAs,
            row.names = FALSE,
            col.names = FALSE)

# results = topTags(et, nrow(et))$table # Line exclusive of EdgeR analysis
dataframe_save = as.data.frame(resOrdered)
dataframe_save = cbind(Feature = rownames(resOrdered), dataframe_save)
colnames(dataframe_save) = c('Feature', 'baseMean', 'log2FC', 'lfcSE', 'stat', 'pvalue', 'qvalue')
write.table(dataframe_save,
            path_output_file,
            row.names = FALSE,
            col.names = TRUE, 
            sep = '\t')
ptm('Done')