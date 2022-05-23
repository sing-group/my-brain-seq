## Script input parameters:
##  1.- path_counts: path to the counts file (without the first line).
##  2.- path_cond: path to the condition file.
##  3.- input_contrast: A line of the contrast file with one contrast.
##  4.- path_output: the output dir for the results.

suppressMessages(library(edgeR))

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

cat("  ------------------------------------------------------------------", '\n')
cat("[PIPELINE -- edger]: Input counts file: ", as.character(args[1]), '\n')
cat("[PIPELINE -- edger]: Input condition file: ", as.character(args[2]), '\n')
cat("[PIPELINE -- edger]: Contrast: ", as.character(args[3]), '\n')
cat("[PIPELINE -- edger]: Output directory: ", as.character(args[4]), '\n')
cat("  ------------------------------------------------------------------", '\n')

#DETECT FACTORS AND LABELS
#read the input files
contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=')
conditions_table = read.delim(path_cond)
#get the reference factor and the condition (removing trailing/leading whitespaces)
contrast = strsplit(as.character(contrast_table$name), '-')[[1]]
ref_factor = trimws(contrast[2])
cond_factor = as.character(trimws(contrast[1]))
print(paste0('[PIPELINE -- edger]: Reference factor: ', ref_factor))
print(paste0('[PIPELINE -- edger]: Condition factor: ', cond_factor))

# #get the label of the reference factor
# ref_index=match(ref_factor, conditions_table$condition)
# ref_label = as.character(conditions_table$label[ref_index])
ref_contrast_label=trimws(strsplit(as.character(rownames(contrast_table)), '-')[[1]][2])
# #get the label of the condition factor
# cond_index = match(cond_factor, conditions_table$condition)
# cond_label = as.character(conditions_table$label[cond_index])
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

#SUBSET THE DESIRED FACTORS
print('[PIPELINE -- edger]: Subsetting desired factors')
# samples with the desired factors to compare
des_samples = as.vector(subset(conditions_table, conditions_table$condition == ref_factor | conditions_table$condition == cond_factor )$name)
# remove the undesired samples from the counts
x=x[,(names(x) %in% des_samples)]
# remove the undesired samples from the conditions_table
conditions_table=conditions_table[(conditions_table$name %in% des_samples),]

#EDGER ANALYSIS
print('[PIPELINE -- edger]: Performing the differential expression with EdgeR')
print('')
#build the group using ref_label as reference factor (first position in vector)
# group = factor(user_labels, levels = c(ref_label, cond_label))
group = factor(conditions_table$condition, levels = c(ref_factor, cond_factor))
#build the DGEList object for the EdgeR analysis
y = DGEList(counts=x, group=group)
#filter out lowly expressed miRNAs
keep = filterByExpr(y)
y = y[keep,,keep.lib.sizes=FALSE]
#normalization by library size (uses TMM)
y = calcNormFactors(y)
#build the model matrix
design = model.matrix(~group)
rownames(design) = colnames(x)
#estimate dispersion
y = estimateDisp(y,design) 
#Differential expression analysis using an exact test
et = exactTest(y)

#RESULTS
#summary of the results
print(topTags(et))
print(summary(decideTests(et)))
#save results
output_label = trimws(as.character(rownames(contrast_table)))
output_tag = paste('EdgeR_', output_label, sep = '')
output_file = paste(output_tag, '.tsv', sep = '')
path_output_file = paste(path_output, output_file, sep='')

results = topTags(et, nrow(et))$table
dataframe_save = as.data.frame(results)
dataframe_save = cbind(Feature = rownames(results), dataframe_save)
write.table(dataframe_save, path_output_file, row.names = FALSE, col.names = TRUE, sep = '\t')