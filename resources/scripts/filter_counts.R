#Command line arguments that will be used from pipeline.xml
args = commandArgs(trailingOnly = TRUE)

#INPUT ARGUMENTS
#counts file (all_counts.txt)
path_counts=as.character(args[1])
#condition file
path_cond=as.character(args[2])
#contrast file
path_contrast=as.character(args[3])
#output directory
output_dir=as.character(args[4])

#LOAD FILES
#read the input files
all_counts = read.delim(path_counts, skip = 1, row.names = 1, check.names=FALSE); counts = all_counts[-1:-5]
conditions_table = read.delim(path_cond)
contrast_table = read.delim(path_contrast, sep = '=')

#DETECT FACTORS AND LABELS
#get the condition and control
contrast = strsplit(as.character(contrast_table$name), '-')[[1]]
ref_factor = trimws(contrast[2])
cond_factor = as.character(trimws(contrast[1]))

#SUBSET DESIRED FACTORS
# samples with the desired factors
des_samples = as.vector(subset(conditions_table, conditions_table$condition == ref_factor | conditions_table$condition == cond_factor )$name)
# desired samples on counts
counts=counts[,(names(counts) %in% des_samples)]
# desired samples on conditions_table
conditions_table=conditions_table[(conditions_table$name %in% des_samples),]

#WRITE THE OUTPUT
# counts
path_output_counts = paste0(output_dir, '/', 'counts_', cond_factor, '-', ref_factor, '.tsv')
counts = cbind(rownames(counts), data.frame(counts, row.names=NULL, check.names=FALSE))
colnames(counts)[1] = 'Feature' 
write.table(counts, path_output_counts, row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
# conditions
path_output_conditions = paste0(output_dir, '/', 'conditions_', cond_factor, '-', ref_factor, '.tsv')
write.table(conditions_table, path_output_conditions, row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
