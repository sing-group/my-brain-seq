#-------------------------------------------------------------------------------
#                                INPUTS
#-------------------------------------------------------------------------------
## Script input parameters:
##  1.- path_counts : path to the "all-counts.txt" file.
##  2.- path_conditions: path to "conditions_file.txt".
##  3.- path_output : the path to save the table (/wd/output/6_x/y-z/pipel/)
##  4.- software : the software of the current iteration (DESeq2/EdgeR/DESeq2-EdgeR).
##  5.- input_contrast : a line of the contrast_file.

args = commandArgs(trailingOnly = TRUE)
path_counts     = as.character(args[1])
path_conditions = as.character(args[2])
path_output     = as.character(args[3])
software        = as.character(args[4])
input_contrast  = as.character(args[5])

#-------------------------------------------------------------------------------
#                           MANUAL PARAMETERS
#-------------------------------------------------------------------------------
# colname of the main factor, e.g.: condition (control, disease)
colname_factor = "condition"

# save an output file
save_output = TRUE

#-------------------------------------------------------------------------------
#                     LIBRARIES, FUNCTIONS, HEADER
#-------------------------------------------------------------------------------
#libraries
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))

# function to print messages
ptm = function(text){
  header = paste0('[MBS | hclust | run_deseq2_normalization.R | ', software, ']: ')
  cat(paste(header, text), sep='\n')
}

# function to change the colnames of "all-counts.txt" from paths to sample names
fix_counts_colnames = function(cts, conditions){
  cts_temp = cts
  for (i in rownames(conditions)){
    index = grep(i, colnames(cts_temp))
    colnames(cts_temp)[index] = i
  }
  return(cts_temp)
}

# function to return DESeq2 normalized data
getNormalizedCounts = function(countData, colData, design) {
  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = design
  )
  dds <- estimateSizeFactors(dds)
  return (counts(dds, normalized=TRUE))
}

# header
ptm("============================================================================")
ptm('             [MBS | hclust | run_deseq2_normalization.R              ')
ptm('............................................................................')
ptm(paste0("  path_counts     :   ", as.character(args[1])))
ptm(paste0("  path_conditions :   ", as.character(args[2])))
ptm(paste0("  path_output     :   ", as.character(args[3])))
ptm("============================================================================")

#-------------------------------------------------------------------------------
#                          READING INPUT FILES
#-------------------------------------------------------------------------------
ptm('Reading input files')
cts = read.csv(file = path_counts,
               sep = "\t",
               row.names = 1,
               skip = 1,
               header = TRUE,
               check.names=FALSE)

conditions = read.csv(file = path_conditions, 
                   sep = "\t",
                   row.names = 1,
                   header = TRUE)

contrast_table = read.delim(text = str_replace_all(input_contrast, " ", ""),
                            sep = '=',
                            col.names = c('label', 'contrast'),
                            header = FALSE)

#-------------------------------------------------------------------------------
#                           PREPARING DATA
#-------------------------------------------------------------------------------
ptm('Preparing count file for normalization')
# fix cts colnames
cts = fix_counts_colnames(cts, conditions)

# keep only cts columns in conditions
cts = cts %>% select(rownames(conditions))

# remove all columns of conditions but the colname with the factor
conditions = as.data.frame(conditions[colname_factor])

# replace NA/empty values for "missing" to treat as group apart                # See https://support.bioconductor.org/p/111785/
conditions[conditions==""]<-NA
conditions = conditions %>% 
  dplyr::mutate_all(~replace(., is.na(.), 'missing'))

# checking that the cts and annotations are correct
test_col = all(rownames(conditions) %in% colnames(cts))
ptm(as.character(paste0("Rownames of 'condition_file.txt' are the same as colnames of 'all-counts.txt': ", test_col)))

# test if colnames of cts are in the same order that rownames of conditions
if (all(rownames(conditions) == colnames(cts)) == FALSE) {
  cts = cts[, rownames(conditions)] #if not, arranges them
}

#-------------------------------------------------------------------------------
#                         DESEQ2 NORMALIZATION
#-------------------------------------------------------------------------------
ptm('Normalizing counts by median of ratios method with DESeq2')
# run the normalization without specifying a model
cts_normalized = getNormalizedCounts(cts, conditions, ~1)
dataframe_save = as.data.frame(cts_normalized)

#-------------------------------------------------------------------------------
#                        PREPARE OUTPUT TABLE
#-------------------------------------------------------------------------------
# SUBSET SAMPLES OF THE CURRENT CONTRAST FROM COUNTS
# factors to subset from conditions (using the contrast)
factors_subset = as.vector(strsplit(contrast_table$contrast, '-')[[1]])
# get a list of all the samples to use
use = rownames(conditions %>% filter(get(colname_factor) %in% factors_subset))
# subset counts
dataframe_save = dataframe_save %>% select(any_of(use))
dataframe_save = cbind(Geneid = rownames(cts_normalized), dataframe_save)

#-------------------------------------------------------------------------------
#                            SAVE OUTPUT
#-------------------------------------------------------------------------------
if ( save_output == TRUE ){
  ptm('Saving output files')
  
  # prepare the filename of the output
  output_filename = paste0('all-counts_',
                           software,
                           '_', 
                           contrast_table$label,
                           '_normalized.txt')
  
  path_output_file = paste0(path_output, '/', output_filename)
  write.table(dataframe_save, path_output_file,
              row.names = FALSE, col.names = TRUE, sep="\t")
}
