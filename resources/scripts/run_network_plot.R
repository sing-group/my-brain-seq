#-------------------------------------------------------------------------------
#                                INPUTS
#-------------------------------------------------------------------------------
## Script input parameters:
##  1.- path_tarbase_db : a file downloaded from Tarbase.
##  2.- path_reactome_db: a file downloaded from Reactome database.
##  3.- path_enrichment_table : the file resulting from a DEA of the pipeline.
##  4.- path_reactome_interaction : a file downloaded from Reactome database.
##  5.- input_contrast  : a line of the contrast_file.
##  6.- software        : the software that produced the "dea_result_file" (DESeq2/EdgeR/DESeq2-EdgeR).
##  7.- db_organism     : the name of the organism, e.g.: 'Homo sapiens'
##  8.- path_output     : path to save the output files.

args = commandArgs(trailingOnly = TRUE)
path_tarbase_db            = as.character(args[1])
path_reactome_db           = as.character(args[2])
path_enrichment_table      = as.character(args[3])
path_reactome_interaction  = as.character(args[4])
input_contrast             = as.character(args[5])
software                   = as.character(args[6])
db_organism                = as.character(args[7])
path_output                = as.character(args[8])

# the number of row from enrichment table from which the data should be taken
# this row number should be taken from a higer-low-pvalue-ordered table (top X).
nth_reactome_id = 1

#-------------------------------------------------------------------------------
#                     LIBRARIES, FUNCTIONS, HEADER
#-------------------------------------------------------------------------------
#libraries
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(igraph))
suppressMessages(library(networkD3))
suppressMessages(library(htmlwidgets))
suppressMessages(library(htmltools))

# function to print messages
ptm = function(text){
  header = paste0('[PIPELINE -- network -- run_network_plot.R > ', software, ']: ')
  cat(paste(header, text), sep='\n')
}

# header
ptm("============================================================================")
ptm('                [PIPELINE -- network -- run_network_plot.R                  ')
ptm('............................................................................')
ptm(paste0("  path_tarbase_db :             ", as.character(args[1])))
ptm(paste0("  path_reactome_db:             ", as.character(args[2])))
ptm(paste0("  path_enrichment_table:        ", as.character(args[3])))
ptm(paste0("  path_reactome_interaction:    ", as.character(args[4])))
ptm(paste0("  input_contrast:               ", as.character(args[5])))
ptm(paste0("  software:                     ", as.character(args[6])))
ptm(paste0("  db_organism:                  ", as.character(args[7])))
ptm(paste0("  path_output:                  ", as.character(args[8])))
ptm("============================================================================")

#-------------------------------------------------------------------------------
#                          READING INPUT FILES
#-------------------------------------------------------------------------------
ptm('Reading input files')
enrichment_table = read.csv(file = path_enrichment_table, 
                            sep = "\t",
                            header = TRUE)

reactome_colnames = c('ensembl_id', 'reactome_id', 'url', 'description','evidence','organism')
reactome_db = read.csv(file = path_reactome_db, 
                       sep = "\t",
                       col.names = reactome_colnames,
                       skip = 1)

tarbase_db = read.csv(file = path_tarbase_db,
                       sep = "\t", 
                       header = TRUE)

reactome_interaction = read.csv(file = path_reactome_interaction,               # https://reactome.org/download/current/interactors/reactome.all_species.interactions.tab-delimited.txt
                       sep = "\t", 
                       header = TRUE)

contrast_table = read.delim(text = str_replace_all(input_contrast, " ", ""),
                            sep = '=',
                            col.names = c('label', 'contrast'),
                            header = FALSE)

#-------------------------------------------------------------------------------
#                   TEST FOR VALID ORGANISM AND FILTER
#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------- TEST IF "db_organism" HAS A VALID VALUE
ptm('Checking if "db_organism" is valid')
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
  species_msg = paste0(options$species, collapse = ', ')
  error_msg = paste0('[ERROR]: "', db_organism, 
                     '" is not a valid value for "db_organism", it must be one of the following: ',
                     species_msg, '.')
  stop(error_msg)
}

#-------------------------------------------------------------------------------
#                        PREPARING THE TABLE
#-------------------------------------------------------------------------------
ptm('Preparing the interaction table')
#------------------------------------------------------------------------------- GET THE SELECTED PATHWAY ID
# arrange the enrichment table by pvalues (lower to higher)
enrichment_table = enrichment_table %>% arrange(pvalue)
# get the id of the pathway
sel_pathway_id = enrichment_table$reactome_id[nth_reactome_id]

#-- GET THE MIRNAS FROM THE SELECTED REACTOME ID (selected by "nth_reactome_id")
mirnas = strsplit(as.character(enrichment_table$DE_mirnas[nth_reactome_id]),
                  split = "|", fixed = TRUE)
mirnas = as.data.frame(mirnas, col.names = c('mirna'))

#------------------------------------------------------------------------------- GET THE GENES FROM THE SELECTED REACTOME ID
genes_in_pathway = reactome_db %>%
  filter(reactome_id == sel_pathway_id) %>%
  select(ensembl_id)

#------------------------------------------------------------------------------- GET THE MIRNA-GENE INTERACTIONS
# keep just geneID, miRNA target and direct or indirect columns on tarbas_db
tarbase_db_filtered = tarbase_db %>% select(geneId, mirna, direct_indirect)
# keep only miRNAs on the miRNA profile
table_mirna_gene = inner_join(genes_in_pathway, tarbase_db_filtered, 
                              by = c("ensembl_id" = "geneId")) %>%
  filter(mirna %in% mirnas$mirna)

#------------------------------------------------------------------------------- PREPARE THE REACTOME INTERACTIONS
# remove empty values from interactors
table_reactome_interactions = reactome_interaction %>% 
  filter(Interactor.1.Ensembl.gene.id != '-' & Interactor.2.Ensembl.gene.id != '-')

# select variables of interest and separate Interactor.1 and 2 columns in
# individual rows, one per ENSEMBL id.
table_reactome_interactions = table_reactome_interactions %>% 
  select(Interactor.1.Ensembl.gene.id,
         X..Interactor.1.uniprot.id,
         Interactor.2.Ensembl.gene.id,
         Interactor.2.uniprot.id,
         Interaction.type) %>%
  separate_rows(Interactor.1.Ensembl.gene.id, sep ="\\|") %>%
  separate_rows(Interactor.2.Ensembl.gene.id, sep ="\\|")

# rename columns and arrange the data by removing extra text
table_reactome_interactions = table_reactome_interactions %>%
  rename(interactor_1 = Interactor.1.Ensembl.gene.id,
         uniprot_1 = X..Interactor.1.uniprot.id,
         interactor_2 = Interactor.2.Ensembl.gene.id,
         uniprot_2 = Interactor.2.uniprot.id) %>%
  mutate(interactor_1 = str_remove(interactor_1, "ENSEMBL:")) %>%
  mutate(interactor_2 = str_remove(interactor_2, "ENSEMBL:")) %>%
  mutate(uniprot_1 = str_remove(uniprot_1, "uniprotkb:")) %>%
  mutate(uniprot_2 = str_remove(uniprot_2, "uniprotkb:"))

# select only the interactions present in the pathway id selected
table_reactome_interactions_filtered = table_reactome_interactions %>%
  filter(interactor_1 %in% genes_in_pathway$ensembl_id |
           interactor_2 %in% genes_in_pathway$ensembl_id)

# remove unwanted reactome references in uniprot colums         
table_reactome_interactions_filtered = table_reactome_interactions_filtered %>%
  filter(!grepl("reactome:", uniprot_1, fixed = TRUE),
           !grepl("reactome:", uniprot_2, fixed = TRUE))

#------------------------------------------------------------------------------- MAP ENSEMBL-MIRNA TO UNIPROT-MIRNA
# create the mapping table using the filtered table of reactome interactions
ensembl_uniprot = bind_rows(table_reactome_interactions_filtered[1:2],
                        table_reactome_interactions_filtered[3:4])
ensembl_uniprot = ensembl_uniprot %>% select(interactor_1, uniprot_1) %>%
  rename(ensembl = interactor_1, uniprot = uniprot_1) %>%
  distinct(ensembl, uniprot, .keep_all = FALSE)

# convert the gene-mirna table to uniprot-mirna
uniprot_mirna_filtered = table_mirna_gene %>% 
  inner_join(ensembl_uniprot, by = c("ensembl_id" = "ensembl")) %>%
  select(mirna, uniprot, direct_indirect) %>%
  rename(interactor_1 = mirna,
         interactor_2 = uniprot,
     Interaction.type = direct_indirect)

#------------------------------------------------------------------------------- BUILD THE CYTOSCAPE TABLE
# subset columns of table_reactome_interactions_filtered
table_reactome_interactions_filtered = table_reactome_interactions_filtered %>%
  select(uniprot_1, uniprot_2, Interaction.type) %>%
  rename(interactor_1 = uniprot_1,
         interactor_2 = uniprot_2)
cytoscape_table = bind_rows(table_reactome_interactions_filtered,
                        uniprot_mirna_filtered) %>%
  distinct()

#-------------------------------------------------------------------------------
#                              NETWORK FIGURE
#-------------------------------------------------------------------------------
ptm('Preparing the network figure')
output_file = paste0(path_output, 'cytoscape_table.tsv')
write.table(cytoscape_table, output_file, row.names = FALSE, 
            col.names = TRUE, sep = '\t', quote = FALSE)

# create for the graph
links <- data.frame(
  src=cytoscape_table$interactor_1,
  target=cytoscape_table$interactor_2,
  interaction = cytoscape_table$Interaction.type
)

#------------------------------------------------------------------------------- CREATES NODE INFO
nodes <- data.frame(
  name=unique(c(links$src, links$target)),
  stringsAsFactors = FALSE
)
nodes$id=0:(nrow(nodes) - 1)

# make a grouping variable that will match to colours
nodes = nodes %>% mutate(mirna = grepl('-miR-|-let-', name, fixed = FALSE))

#------------------------------------------------------------------------------- CREATES EDGE INFO
edges = links %>%
  select(-interaction) %>%
  left_join(nodes, by = c("src" = "name")) %>%
  select(-src, -mirna) %>%
  rename(source = id) %>%
  left_join(nodes, by = c("target" = "name")) %>%
  select(-target, -mirna) %>%
  rename(target = id)

edges$width <- 1

# create the network object
network <- graph_from_data_frame(d=links, vertices=nodes, directed=F)

ColourScale <- 'd3.scaleOrdinal()
            .domain(["TRUE", "FALSE"])
           .range(["#36B0AD", "#DF9F1F"]);'

# interactive network
p <- forceNetwork(Links = edges, Nodes = nodes,
                  Source = "source",
                  Target = "target",
                  NodeID ="name",
                  Group = "mirna",
                  zoom = T,
                  opacity = 0.9,
                  charge = -50,
                  colourScale = JS(ColourScale),
                  opacityNoHover = 0.8)

#------------------------------------------------------------------------------- ADD NETWORK CAPTIONS
#title and subtitle
title_label       = paste0(sel_pathway_id, ':')
description_label = reactome_db[reactome_db$reactome_id == sel_pathway_id,]$description[1]
p <- htmlwidgets::prependContent(p, htmltools::tags$p('', style="line-height:0;font-size:2px"))
p <- htmlwidgets::prependContent(p, htmltools::tags$h2(title_label, style='display:inline;color:#125251;font-family:helvetica'))
p <- htmlwidgets::prependContent(p, htmltools::tags$h3(description_label, style='display:inline;color:#348180;font-family:helvetica'))
p <- htmlwidgets::prependContent(p, htmltools::tags$hr(''))
# contrast
p <- htmlwidgets::prependContent(p, htmltools::tags$strong('Contrast: ', style='color:#1B6564;font:helvetica'))
p <- htmlwidgets::prependContent(p, htmltools::tags$em(contrast_table$label))
# software
p <- htmlwidgets::prependContent(p, htmltools::tags$p('', style="line-height:0;font-size:2px"))
p <- htmlwidgets::prependContent(p, htmltools::tags$strong('Software: ', style='color:#1B6564;font:helvetica'))
p <- htmlwidgets::prependContent(p, htmltools::tags$em(software))
# qvalue
qvalue_label = round(enrichment_table[nth_reactome_id,]$qvalue, 4)
p <- htmlwidgets::prependContent(p, htmltools::tags$p('', style="line-height:0;font-size:2px"))
p <- htmlwidgets::prependContent(p, htmltools::tags$em('q-value: ', style='color:#838383;font:helvetica;display:inline'))
p <- htmlwidgets::prependContent(p, htmltools::tags$em(qvalue_label, style='color:#838383'))

#-------------------------------------------------------------------------------
#                                 OUTPUT
#-------------------------------------------------------------------------------
ptm('Saving output files')

label = paste0('_', sel_pathway_id, '_', software, '_', contrast_table$label)

# save table compatible with Cytoscape
output_table_path = paste0(path_output, 'cytoscape_table', label, '.tsv')
write.table(cytoscape_table, output_table_path, row.names = FALSE, 
            col.names = TRUE, sep = '\t', quote = FALSE)

# save network widget
output_network_path = paste0(path_output, 'network_', label, '.html')
suppressWarnings(saveWidget(p, file=output_network_path))

# remove temporary files of the html widget
unlink(paste0(path_output, 'network_', label), recursive = TRUE)

ptm('Done')
