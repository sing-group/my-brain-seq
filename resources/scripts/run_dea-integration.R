## Script input parameters:
##  1.- path_deseq        : the file resulting from the DESeq2 analysis.
##  2.- path_edger        : the file resulting from the EdgeR analysis.
##  3.- input_contrast    : a line of the contrast_file.
##  4.- path_output       : output path to save the results.
##  5.- path_output_pipel : output path to save additional results.

# to print messages
ptm = function(text){
  header = '[PIPELINE -- dea-integration -- run_dea-integration.R]:'
  cat(paste(header, text), sep='\n')
}

#-------------------------------------------------------------------------------
#                                 INPUTS
#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
path_deseq = as.character(args[1])
path_edger = as.character(args[2])
input_contrast = as.character(args[3])
path_output = as.character(args[4])
path_output_pipel = as.character(args[5])

q_value_filter = 0.05
logFC = 0.5

#-------------------------------------------------------------------------------
#                                FUNCTIONS
#-------------------------------------------------------------------------------
ptm('Loading functions')
# to get the reference and condition factors from the input_contrast (pipeline)
get_contrast_factors = function(input_contrast) {
  contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=')
  contrast = strsplit(as.character(contrast_table$name), '-')[[1]]
  ref_factor = trimws(contrast[2])
  cond_factor = as.character(trimws(contrast[1]))
  result = list("reference" = ref_factor, "condition" = cond_factor)
}

# to print messages
ptm = function(text){
  header = '[PIPELINE -- dea-integration -- run_dea-integration.R]:'
  cat(paste(header, text), sep='\n')
}

# to print bare messages
pt = function(text){
  cat(text)
}

#-------------------------------------------------------------------------------
#                               READ FILES
#-------------------------------------------------------------------------------
ptm('Reading DEA files')
suppressMessages(library('data.table'))
suppressMessages(library('dplyr'))

# reading
deseq = read.delim(path_deseq)
edger = read.delim(path_edger)

#-------------------------------------------------------------------------------
#                           PREPARE DEA TABLES
#-------------------------------------------------------------------------------
ptm('Preparing DEA tables for the integration')
# keep only Feature, pvalue and qvalue columns and rename
pval_d = 'pvalue_deseq'; qval_d = 'qvalue_deseq'; fc_d = 'log2FC_deseq'
pval_e = 'pvalue_edger'; qval_e = 'qvalue_edger'; fc_e = 'log2FC_edger'

deseq = deseq %>%
  select(Feature, pvalue, qvalue, log2FC)
colnames(deseq) = c('Feature', pval_d, qval_d, fc_d)

edger = edger %>%
  select(Feature, pvalue, qvalue, log2FC)
colnames(edger) = c('Feature', pval_e, qval_e, fc_e)

#-------------------------------------------------------------------------------
#                   DESEQ-EDGER COINCIDENCES AFTER QVALUE FILTER
#-------------------------------------------------------------------------------
# First apply the qvalue filter to the dea files independently, then search for
# miRNA coincidences

ptm('Finding coincidences between DESeq2 and EdgeR after filtering by qvalue')
# filter the files by qvalue
deseq_filtered = deseq %>% filter(qvalue_deseq < q_value_filter)
edger_filtered = edger %>% filter(qvalue_edger < q_value_filter)

# filter the files by log2FC
deseq_filtered = deseq_filtered %>%
  filter(log2FC_deseq <= -logFC | log2FC_deseq >= logFC)
edger_filtered = edger_filtered %>%
  filter(log2FC_edger <= -logFC | log2FC_edger >= logFC)

# search for coincidences miRNA between the two tables, then averages p/qvalues
coincidences_full = inner_join(deseq_filtered, edger_filtered, by = 'Feature') %>%
  mutate(pvalue = (pvalue_deseq + pvalue_edger)/2, 
         qvalue = (qvalue_deseq + qvalue_edger)/2)

coincidences = coincidences_full %>% select(Feature, pvalue, qvalue)

# print a summary and save venn table
if (length(coincidences$Feature) > 0){
  # prepare the table for the venn diagram
  venn_d = deseq_filtered$Feature
  venn_e = edger_filtered$Feature
  max_length = max(c(length(venn_d), length(venn_e))) 
  venn_table = data.frame(
    deseq = c(venn_d,
              rep('', max_length - length(venn_d))),
    edger = c(venn_e,
              rep('', max_length - length(venn_e))))
  # summary
  pt('\n')
  pt('MiRNAs differentially expressed in DESeq and EdgeR results\n')
  cat(coincidences$Feature, sep = '\n')
  pt('\n')
} else {
  ptm(paste0('[WARNING]: No coincidences between EdgeR and DESeq2 results for a q-value < ', q_value_filter))
  ptm('[WARNING]: Venn plot skipped"')
}

#-------------------------------------------------------------------------------
#                            P/QVALUE AVERAGE
#-------------------------------------------------------------------------------
# Join the deseq and edger files by miRNA (keeping only coincidences) and then
# calculates the average qvalue and pvalue per miRNA, e.g.:
#    Average pvalue for the miRNA "a":
#        ave_pvalue(a) = (pvalue_deseq(a) + pvalue_edger(a))/2

ptm('Averaging miRNAs DESeq2 and EdgeR pvalues and q-values')
# perform the integration
integrated_full = deseq %>%
  inner_join(edger, by = 'Feature') %>%
  mutate(pvalue = (pvalue_deseq + pvalue_edger)/2, 
         qvalue = (qvalue_deseq + qvalue_edger)/2,
         log2FC = (log2FC_deseq + log2FC_edger)/2)

# table with just Features, pvalue, qvalue and log2FC columns
integrated = integrated_full %>%
  select(Feature, pvalue, qvalue, log2FC)

# table with just the DEmiRNAs
DEmiRNAs = integrated %>% 
  filter(qvalue < q_value_filter) %>%
  filter(abs(log2FC) >= logFC)

# print a summary
if (length(integrated$Feature) > 0){
  pt('\n')
  pt('Top integrated values\n')
  print(head(integrated %>% arrange(qvalue), n = 10))
  pt('\n')
} else {
  ptm('[WARNING]: No coincidences between EdgeR and DESeq2 results')
}

#-------------------------------------------------------------------------------
#                            SAVE RESULTS
#-------------------------------------------------------------------------------
ptm('Saving integrated results')
# build the label for the contrast
factors = get_contrast_factors(input_contrast)
contrast_name = paste0(factors$condition, '-', factors$reference)
# build the paths for the output files
path_output_integrated = paste0(path_output, '/', 'DEmiRNAs_', contrast_name, 
                                '_deseq-edger_integrated', '.tsv')
path_output_integrated_full = paste0(path_output_pipel, '/', 'DEmiRNAs_', 
                                     contrast_name, 
                                     '_deseq-edger_integrated_full', '.tsv')
path_output_DEmiRNAs = paste0(path_output, '/', 'differentially-expressed-miRNAs_',
                              contrast_name,
                              '_deseq-edger_integrated', '.tsv')
path_output_coincidences = paste0(path_output_pipel, '/', 'DEmiRNAs_', 
                                  contrast_name, 
                                  '_deseq-edger_coincidences', '.tsv')
path_output_venn = paste0(path_output_pipel, '/', 'DEmiRNAs_', 
                                  contrast_name, 
                                  '_deseq-edger_venn', '.tsv')
# write integrated
if (length(integrated$Feature) > 0){
write.table(integrated, path_output_integrated, 
            row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
write.table(integrated_full, path_output_integrated_full, 
              row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
}

# write integrated DE miRNAs
write.table(DEmiRNAs$Feature, path_output_DEmiRNAs, 
            row.names = FALSE, col.names = FALSE, sep="\t", quote = TRUE)

# write coincidences
if (length(coincidences$Feature) > 0){
  write.table(coincidences, path_output_coincidences, 
              row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
  write.table(venn_table, path_output_venn, 
              row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
}

ptm('Done')