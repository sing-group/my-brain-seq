## Script input parameters:
##  1.- path_counts: path to the counts file (without the first line).
##  2.- path_cond: path to the condition file.
##  3.- input_contrast: A line of the contrast file with one contrast.
##  4.- path_output: the output dir for the results.

suppressMessages(library('edgeR'))

# to print without "[1]"
pt = function(text){
  cat(text, sep='\n')
}
# to print messages
ptm = function(text){
  header = '[MBS | edger | run_edger_exact-test.R]:'
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
padj=as.numeric(args[5])
logFC=as.numeric(args[6])

ptm("======================================================")
ptm('    [MBS | edger]: run_edger_exact-test.R       ')
ptm('......................................................')
ptm(paste0("  Input counts file:    ", as.character(args[1])))
ptm(paste0("  Input condition file: ", as.character(args[2])))
ptm(paste0("  Contrast:             ", as.character(args[3])))
ptm(paste0("  Output path:          ", as.character(args[4])))
ptm('  --')
ptm(paste0("  qvalue:          ", as.character(args[5])))
ptm(paste0("  log2FC:          ", as.character(args[6])))
ptm("======================================================")

#DETECT FACTORS AND LABELS
#read the input files
contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=')
conditions_table = read.delim(path_cond)
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
x = read.delim(path_counts, row.names="Geneid")
x = x[-1:-5]
#builds a list with the samples labeled in order
user_labels = c() #e.g.: C,C,PA,PA
sampleRepl = c() #e.g.: C01,C02,PA01,PA02
for (i in conditions_table$name){
  #find the index of each label in "x" by using grep
  index = grep(i, colnames(x))
  #change "x" colnames to its corresponding label
  colnames(x)[index] = i
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
x=x[,(names(x) %in% des_samples)]
# remove the undesired samples from the conditions_table
conditions_table=conditions_table[(conditions_table$name %in% des_samples),]

#-------------------------------------------------------------------------------
#                             EDGER ANALYSIS
#-------------------------------------------------------------------------------
ptm('Performing the differential expression with EdgeR')
#build the group using ref_label as reference factor (first position in vector)
group = factor(conditions_table$condition, levels = c(ref_factor, cond_factor))
#build the DGEList object for the EdgeR analysis
y = DGEList(counts=x, group=group)

#FILTERING AND NORMALIZATION
#filter out lowly expressed miRNAs
keep = filterByExpr(y)
y = y[keep,,keep.lib.sizes=FALSE]
#normalization by library size (uses TMM)
y = calcNormFactors(y)

#EDGER DESIGN
#builds the edger model matrix
if (length(adjust_for > 0)){
  adjustment = paste('conditions_table$', adjust_for, collapse=' + ', sep = '')
  design = paste("~group +", adjustment)
  design = model.matrix(formula(design))
  # for printing model
  print_adjustment = paste(adjust_for, collapse=' + ')
  print_design = paste("~group +", print_adjustment)
} else {
  print_design = "~group"
  design = model.matrix(~group)
}
rownames(design) = colnames(x)

ptm(paste0('Model: ', print_design))
pt('')

#estimate dispersion
y = estimateDisp(y,design) 
#Differential expression analysis using glmFit and glmLRT
fit = glmFit(y, design)
# Testing the effect of 'group' with coef=2
lrt = glmLRT(fit, coef=2) 

#RESULTS
#Get a list with the DE miRNAs (p.value cut-off refers to adjusted pvalue)
DE_miRNAs = topTags(lrt, n=Inf, p.value=padj, adjust.method="BH", sort.by="PValue")

if (length(DE_miRNAs) != 0){
  DE_miRNAs = data.frame(DE_miRNAs$table)
  DE_miRNAs = DE_miRNAs[abs(DE_miRNAs$logFC) >= logFC,] 
} else {
  print(paste0('[PIPELINE | edger]: [INFO] No differentially expressed miRNAs on contrast', contrast_table$name))
}

#summary of the results
print(head(DE_miRNAs, 10))

#save results
pt('')
ptm('Saving results')
output_tag = paste('EdgeR_', cond_contrast_label, '-', ref_contrast_label, sep = '')
output_file = paste(output_tag, '.tsv', sep = '')
path_output_file = paste(path_output, output_file, sep='')
#save DE miRNA list
output_file_DEmiRNAs = paste('differentially-expressed-miRNAs_', output_file, sep = '')
path_output_file_DEmiRNAs = paste(path_output, output_file_DEmiRNAs, sep='')
write.table(row.names(DE_miRNAs),
            path_output_file_DEmiRNAs,
            row.names = FALSE,
            col.names = FALSE)
#save the whole EdgeR table
results = topTags(lrt, nrow(lrt))$table
dataframe_save = as.data.frame(results)
dataframe_save = cbind(Feature = rownames(results), dataframe_save)
colnames(dataframe_save) = c('Feature', 'log2FC', 'logCPM', 'LR', 'pvalue', 'qvalue') # logFC = log2FC, see: https://www.biostars.org/p/303806/#303829
dataframe_save = dataframe_save[, !(colnames(dataframe_save) %in% c('LR'))]
write.table(dataframe_save, path_output_file, row.names = FALSE, col.names = TRUE, sep = '\t')
ptm('Done')