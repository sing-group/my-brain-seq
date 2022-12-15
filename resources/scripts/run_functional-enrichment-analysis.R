# Performs a functional enrichment analysis using the method improved for miRNAs
# described in paper: 
#
#     10.1093/nar/gkv249, "Pathway analysis from lists of
#     microRNAs: common pitfalls and alternative strategy".
# 
#     and in the following RPubs article:
#       (https://rpubs.com/jrgonzalezISGlobal/enrichment)
#
# References for the Tarbase and Reactome files: 
#   - tarbase_db_file :
#       (https://dianalab.e-ce.uth.gr/html/diana/web/index.php?r=tarbasev8%2Fdownloaddataform)
#   - reactome_db_file:
#       (https://reactome.org/download/current/Ensembl2Reactome.txt)
#   - db_organism     :
#       ()

## Script input parameters:
##  1.- tarbase_db_file : a file downloaded from Tarbase.
##  2.- reactome_db_file: a file downloaded from Reactome database.
##  3.- dea_result_file : the file resulting from a DEA of the pipeline.
##  4.- input_contrast  : a line of the contrast_file.
##  5.- software        : the software that produced the "dea_result_file" (DESeq2/EdgeR/DESeq2-EdgeR).
##  6.- db_organism     : the name of the organism, e.g.: 'Homo sapiens'
##  7.- path_output     : path to save the output files.

# INPUTS
args = commandArgs(trailingOnly = TRUE)
tarbase_db_file  = as.character(args[1])
reactome_db_file = as.character(args[2])
dea_result_file  = as.character(args[3])
input_contrast   = as.character(args[4])
software         = as.character(args[5])
db_organism      = as.character(args[6])
path_output      = as.character(args[7])

q_value_filter = 0.05

# to print messages
ptm = function(text){
  header = paste0('[PIPELINE -- functional-enrichment -- run_functional-enrichment-analysis.R > ', software, ']: ')
  cat(paste(header, text), sep='\n')
}

# to print bare messages
pt = function(text){
  cat(text)
}

ptm("============================================================================")
ptm(' [PIPELINE -- functional-enrichment]: run_functional-enrichment-analysis.R  ')
ptm('............................................................................')
ptm(paste0("  tarbase_db_file :   ", as.character(args[1])))
ptm(paste0("  reactome_db_file:   ", as.character(args[2])))
ptm(paste0("  dea_result_file:    ", as.character(args[3])))
ptm(paste0("  input_contrast:     ", as.character(args[4])))
ptm(paste0("  software:           ", as.character(args[5])))
ptm(paste0("  db_organism:        ", as.character(args[6])))
ptm(paste0("  path_output:        ", as.character(args[7])))
ptm("============================================================================")

#-------------------------------------------------------------------------------
#                               LOAD TABLES
#-------------------------------------------------------------------------------
# install.packages("SnowballC")

suppressMessages(library('dplyr'))
suppressMessages(library('stringr'))
suppressMessages(library('tidytext'))
suppressMessages(library('SnowballC'))
suppressMessages(library('ggplot2'))

# loads Tarbase db file
tarbase_db = read.delim(tarbase_db_file)
# load Reactome db file
reactome_colnames = c('ensembl_id', 'reactome_id', 'url', 'description','evidence','organism')
reactome_db = read.delim(reactome_db_file, skip = 1, col.names = reactome_colnames)
# load dea result file
dea_result_table = read.delim(dea_result_file)

# TEST IF "db_organism" HAS A VALID VALUE
# get the valid options
options = tarbase_db %>%
  distinct(species) %>%
  filter(species %in% reactome_db$organism)

# make the test
if (db_organism %in% options$species){
  ptm(paste0('"', db_organism, '" is in TarBase and Reactome annotations'))
  ptm(paste0('filtering TarBase annotations by "', db_organism, '"'))
  tarbase_db = tarbase_db %>% filter(species == db_organism)
  ptm(paste0('filtering Reactome annotations by "', db_organism, '"'))
  reactome_db = reactome_db %>% filter(organism == db_organism)
}else{
  error_msg = paste0('[ERROR]: "', db_organism, '" is not a valid value for "db_organism", it must be one of the following: ', options$species)
  stop(error_msg)
}

#-------------------------------------------------------------------------------
#                               FUNCTIONS
#-------------------------------------------------------------------------------
ptm('Loading functions')
get_genes_in_pathway = function(reactomeId){
  result = pathways_table_DE %>% 
    filter(reactome_id == reactomeId) %>%
    select(ensembl_id)
  
  if (identical(result, character(0))){
    result = result$ensembl_id
    return('')
  }
  else {
    return(paste(result$ensembl_id, collapse = "|"))
  }
}

#Get the reference and condition factors from the input_contrast (pipeline)
get_contrast_factors = function(input_contrast) {
  #Get the current contrast
  contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=')
  #get the reference factor and the condition (removing trailing/leading white spaces)
  contrast = strsplit(as.character(contrast_table$name), '-')[[1]]
  ref_factor = trimws(contrast[2])
  cond_factor = as.character(trimws(contrast[1]))
  result = list("reference" = ref_factor, "condition" = cond_factor)
}

# Raises "message" if "column_test" is empty, then stops the script execution
# messages must be stored in a vector.
check_to_continue = function(column_test, message = NA) {
  if (length(column_test) == 0){
    if (!is.na(message[1])){
      ptm(paste0('No ', message[1], ' with a qvalue < ', q_value_filter, '.'))
      ptm(paste0(message[2], ' of "', contrast_label, '" results skipped.'))
      ptm('Done')
    }
    quit(save = 'no')
  }
}

#-------------------------------------------------------------------------------
#                        PREPARE USEFUL VARIABLES
#-------------------------------------------------------------------------------
# build the labels of the contrast
contrast = get_contrast_factors(input_contrast)
contrast_label = paste0(contrast$condition, '-', contrast$reference)

#-------------------------------------------------------------------------------
#                           TARGET PREDICTION
#-------------------------------------------------------------------------------
# NOTE: only targets with "direct_indirect" column equal to "DIRECT" are kept
ptm('Finding miRNA targets using Tarbase annotations')

# conversion database (to convert pathways-genes to pathways-miRNAs)
conversion_db = tarbase_db %>%
  select(geneId, mirna) %>%
  inner_join(reactome_db, by = c('geneId' = 'ensembl_id')) %>%
  select(mirna, reactome_id) %>%
  distinct(mirna, reactome_id)

# find differentially expressed miRNAs (DE miRNAs)
de_mirnas = dea_result_table %>% filter(qvalue < q_value_filter) %>% select(Feature)
colnames(de_mirnas) = 'mirna'

# stops the execution if no DE miRNAs
msg_no_de = c('miRNAs', 'Functional enrichment analysis')
check_to_continue(de_mirnas$mirna, msg_no_de)

# all miRNAs tested for differential expression analysis (DEA)
all_tested_mirna = dea_result_table %>% select(Feature)
colnames(all_tested_mirna) = 'mirna'

#-------------------------------------------------------------------------------
#                           PATHWAY ASSINGMENT
#-------------------------------------------------------------------------------
ptm('Assinging pathways to the DE miRNAs using Reactome annotations')

# pathways of the DE miRNAs 
pathways_table_DE = conversion_db %>%
  inner_join(de_mirnas, by='mirna')

# all pathways of the genes targeted by all the miRNAs tested for DE
pathways_table_ALL_TESTED = conversion_db %>%
  inner_join(all_tested_mirna, by='mirna')

# list with all the different pathways
pathway_list = reactome_db %>%
  distinct(reactome_id)

#-------------------------------------------------------------------------------
#                             ENRICHMENT
#-------------------------------------------------------------------------------
# Transcriptomic study, in which 12,671 genes have been tested for differential
# expression between two sample conditions (total_mirna_tested) and 529 genes
# were found DE (n_DE_mirnas). Among the DE genes, 28 are annotated to a specific
# functional gene set -or pathway- (n_DE_mirnas_in_pathway), which contains in
# total 170 genes (total_mirnas_in_pathway).
#
#                                https://rpubs.com/jrgonzalezISGlobal/enrichment

ptm("Performing the pathway enrichment using Fisher's exact test")

# constant values
total_mirna_tested = length(all_tested_mirna$mirna)
n_DE_mirnas = length(de_mirnas$mirna)

# variables dependent of the pathway
n_DE_mirnas_in_pathway = NA
total_mirnas_in_pathway = NA

# variable container for the enrichment values of the pathways
enrichment_table = pathway_list
enrichment_table['pvalue'] = 0

# fill the container with q-values resulting from a Fisher hypergeometric test
for (pth in 1:dim(enrichment_table)[1]){
  p = pathway_list$reactome_id[pth]
  
  # number of DE mirnas in the pathway
  DE_mirnas_in_pathway = pathways_table_DE %>%
    filter(reactome_id == p)
  n_DE_mirnas_in_pathway = length(DE_mirnas_in_pathway$mirna)
  
  # number of the total genes in the pathway
  total_mirnas_in_pathway = pathways_table_ALL_TESTED %>%
    filter(reactome_id == p)
  total_mirnas_in_pathway = length(total_mirnas_in_pathway$mirna)
  
  # derivated calcs from "variables dependent of the pathway"
  mirnas_not_in_pathway = total_mirnas_in_pathway - n_DE_mirnas_in_pathway
  DE_mirnas_not_in_pathway = n_DE_mirnas - n_DE_mirnas_in_pathway
  total_no_DE_mirnas = total_mirna_tested - n_DE_mirnas
  
  # build contingency table
  deTable <-  matrix(c(n_DE_mirnas_in_pathway, 
                       mirnas_not_in_pathway, 
                       DE_mirnas_not_in_pathway, 
                       total_no_DE_mirnas),
                     nrow = 2,
                     dimnames = list(DE=c("yes","no"),
                                     GeneSet=c("in","out")))
  
  # enrichment score for a pathway
  enrichment_score = fisher.test(deTable, alternative = "greater")
  
  # write the results to the empty container
  enrichment_table$pvalue[pth] = enrichment_score$p.value

  deTable <- NULL
  
}

#-------------------------------------------------------------------------------
#                            PREPARE RESULTS
#-------------------------------------------------------------------------------
ptm('Building the result table')

# remove pathways with p-value = 1 and order by q-value descending
enrichment_table = enrichment_table %>%
  filter(pvalue < 1)

# Benjamini & Hochberg FDR correction (equivalent to a qvalue with lambda = 0.5) https://www.biostars.org/p/128931/
enrichment_table$qvalue = p.adjust(enrichment_table$pvalue, method = 'BH')

# expand the results adding the explanation of each pathway
enrichment_table = reactome_db %>%
  select(reactome_id, description, url) %>%
  distinct(reactome_id, description, url) %>%
  filter(reactome_id %in% enrichment_table$reactome_id) %>%
  inner_join(enrichment_table, by = 'reactome_id') %>%
  arrange(qvalue)

# ADDING THE DE-MIRNAS PRESENT ON EACH PATHWAY
DE_mirnas = dea_result_table %>%
  filter(qvalue < q_value_filter) %>%
  arrange(qvalue)

# add a new column for the DE miRNAs on the enrichment result table
enrichment_table['DE_mirnas'] = ''

for (pth in 1:dim(enrichment_table)[1]){
  p = enrichment_table$reactome_id[pth]
  mirnas_id_in_DE_pathway = conversion_db %>% 
    filter(reactome_id == p & mirna %in% DE_mirnas$Feature) %>%
    distinct(mirna)
  
  mirnas_id_in_DE_pathway = str_c(mirnas_id_in_DE_pathway$mirna, collapse = '|')
  
  enrichment_table$DE_mirnas[pth] = mirnas_id_in_DE_pathway
}

# save the results of the enrichment analysis
output_table_filename = paste0('enrichment_table_', software, '_', contrast_label, '.tsv')
path_table = paste0(path_output, output_table_filename)
write.table(enrichment_table, path_table, row.names = FALSE, col.names = TRUE, sep = '\t')

# prepare the terminal output with top 10 enriched pathways
print_enriched = enrichment_table %>%
  filter(qvalue < q_value_filter) %>%
  select(description, qvalue) %>%
  arrange(qvalue)

# stops the execution if no DE pathways
msg_chart = c('enriched pathways', 'Lolipop chart')
check_to_continue(print_enriched$description, msg_chart)

# truncate the descriptions to 60 characters
print_enriched = print_enriched %>%
  mutate(pathway = str_trunc(description, 60)) %>%
  select(pathway, qvalue)

print_enriched = head(print_enriched, n = 10)
pt('\n')
pt('TOP 10 ENRICHED PATHWAYS\n')
print(print_enriched, right = F)
pt('\n')

#-------------------------------------------------------------------------------
#                             PLOTTING
#-------------------------------------------------------------------------------
ptm('Building the lolipop chart for the results')
  
# LOLIPOP CHART
# get all the text from the description of the enriched terms
text = enrichment_table %>% 
  filter(qvalue < q_value_filter) %>% 
  select(description)

# split the text in words and remove stop words using tidytext
text_table = text %>% 
  unnest_tokens(word, description) %>%
  anti_join(stop_words, by = 'word')

# removes words shorter than 2 characters
text_table = text_table %>% 
  filter(str_length(word) > 2)

# stem the words
stemmed_text = text_table %>%
  mutate(word = wordStem(word))

# computes frequencies
word_frequencies = stemmed_text %>%
  group_by(word) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

# remove words with only one repetition
word_frequencies = word_frequencies %>% 
  filter(n > 1) %>%
  arrange(desc(n))

# keep just the top x repeated words
n_top_words = 50
reduced_words_frequency = head(word_frequencies, n = n_top_words)

# build the title of the graph
graph_subtitle = paste0('Top ', n_top_words, ' repeated words in the description of Reactome pathways (qvalue < ', q_value_filter,')')
graph_caption = paste0("To calculate word repetitions, first Reactome pathways were filtered by a qvalue < ", q_value_filter,". Second, descriptions of those enriched pathway were joined and splitted into words. Third, the resulting text was cleaned by removing english stop words and performing lemmatisation. Finally, repetitions of the same words were calculated. The top ", n_top_words, " are presented in the graph.")
graph_caption = str_wrap(graph_caption, 100)

# Plot the lolipop chart
output_figure_filename = paste0('enriched_reactome_words_', software, '_', contrast_label, '.pdf')
path_figure = paste0(path_output, output_figure_filename)
pdf(file = path_figure)
reduced_words_frequency %>%
  arrange(n) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(word=factor(word, levels=word)) %>%   # This trick update the factor levels
  ggplot( aes(x=word, y=n)) +
  geom_segment( aes(xend=word, yend=0)) +
  geom_point( size=4, color="#FFBB32") +
  coord_flip() +
  xlab("") +
  ylab("repetitions") +
  labs(title = contrast_label,
       subtitle = graph_subtitle,
       caption = graph_caption) +
  theme(plot.caption.position = "plot",
        plot.caption = element_text(hjust = 1))

silence <- dev.off()

ptm('Done')