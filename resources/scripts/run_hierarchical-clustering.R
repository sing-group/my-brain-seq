#-------------------------------------------------------------------------------
#                                 INPUTS
#-------------------------------------------------------------------------------
## Script input parameters:
##  1.- path_hclust_file: path to the hclust.tsv file of a contrast.
##  2.- input_contrast: A line of the contrast file with one contrast.
##  3.- path_output: the output dir for the results.
##  4.- software: the software that produced the results (DESeq2/EdgeR/DESeq2-EdgeR).
##  5.- conditions_file: path to the conditions_file.

args = commandArgs(trailingOnly = TRUE)
path_hclust_file=as.character(args[1])
input_contrast=as.character(args[2])
path_output=as.character(args[3])
software=as.character(args[4])
conditions_file=as.character(args[5])

# CLUSTERING OPTIONS
distance_method = 'euclidean'
clustering_method = 'ward.D2'
cluster_number = 2
traspose_dendrogram = F

#-------------------------------------------------------------------------------
#                                 CODE
#-------------------------------------------------------------------------------
# to print bare text without "[1]"
pt = function(text){
  cat(text, sep='\n')
}
# to print messages
ptm = function(text, sft = software){
  header = '[PIPELINE -- hclust -- run_hierarchical-clustering.R > '
  header = paste0(header, sft, ']: ')
  cat(paste(header, text), sep='\n')
}

ptm("======================================================")
ptm(' [PIPELINE -- hclust]: run_hierarchical-clustering.R  ')
ptm('......................................................')
ptm(paste0("  Hclust file:     ", as.character(args[1])))
ptm(paste0("  Input contrast:  ", as.character(args[2])))
ptm(paste0("  Output dir:      ", as.character(args[3])))
ptm(paste0("  Software:        ", as.character(args[4])))
ptm(paste0("  conditions_file: ", as.character(args[5])))
ptm("======================================================")

ptm('Loading libraries')
suppressMessages(library('dplyr'))
suppressMessages(library('tibble'))
suppressMessages(library('RColorBrewer'))
suppressMessages(library('dendextend'))
suppressMessages(library("gplots"))

get_contrast_factors = function(input_contrast) {
  #Get the reference and condition factors from the input_contrast (pipeline)
  #Get the current contrast
  contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=')
  #get the reference factor and the condition (removing trailing/leading white spaces)
  contrast = strsplit(as.character(contrast_table$name), '-')[[1]]
  ref_factor = trimws(contrast[2])
  cond_factor = as.character(trimws(contrast[1]))
  result = list("reference" = ref_factor, "condition" = cond_factor)
}
#-------------------------------------------------------------------------------
#                         PREPARES THE COUNTS
#-------------------------------------------------------------------------------
ptm('Loading hclust file')
#Loads the counts data
cts = read.csv(file = path_hclust_file,
               sep = "\t",
               row.names = 1,
               header = TRUE,
               check.names=FALSE)

cts = as.data.frame(t(cts))

ptm('Loading conditions_file')
# loads conditions_file
conditions = read.table(
  file = conditions_file,
  sep = "\t",
  header = TRUE
)

# standarize the values of each row (mirnas, scales...)
cts_scaled = as.data.frame(scale(cts))

# traspose if needed
if (traspose_dendrogram == TRUE){
  cts_scaled = t(cts_scaled)
}

#-------------------------------------------------------------------------------
#                             DISTANCE MATRIX
#-------------------------------------------------------------------------------
ptm('Computing distance matrix')
# Calculated the distances using one of the following methods:
#   "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". 
dist_mat <- dist(cts_scaled, method = "euclidean")

#-------------------------------------------------------------------------------
#                         HIERARCHICAL CLUSTERING
#-------------------------------------------------------------------------------
ptm('Building hierarchical clustering')
# Build the hclust choosing a method, options below:
#   method: complete / average / single / ward.D2
hclust_avg <- hclust(dist_mat, method = clustering_method)
# set the threshold for the number of clusters
cut_avg <- cutree(hclust_avg, k = cluster_number)

#-------------------------------------------------------------------------------
#                              DENDROGRAM
#-------------------------------------------------------------------------------
ptm('Building dendrogram')
# plot the dendogram
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, k = cluster_number)

#-------------------------------------------------------------------------------
#                             HEAT MAP
#-------------------------------------------------------------------------------
ptm('Building heat map')
# color scale
color_scale <- c("black", "blue", "green", "yellow", "orange", "red")
pal <- colorRampPalette(color_scale)(100)
x=as.matrix(t(cts_scaled))

#-------------------------------------------------------------------------------
#                         PREPARE LABELS
#-------------------------------------------------------------------------------
ptm('Making Labels')
# make the label with the contrast, for the output files
contrast_factors = get_contrast_factors(input_contrast)
contrast_label = paste0(contrast_factors$condition, '-'
                        , contrast_factors$reference)

# make the title of the plots
title = paste0(software, ': ', contrast_label)

# find the longer label in the plot
a = as.vector(lapply(colnames(cts), nchar))
b = as.vector(lapply(rownames(cts), nchar))
c = unlist(c(a, b))
longer_label = max(c)

# uses the length of the longer label to calculate the canvas margins
margins = round(longer_label - (longer_label / 5))

#-------------------------------------------------------------------------------
#                         EXPORT DENDROGRAM
#-------------------------------------------------------------------------------
ptm('Exporting Dendrogram')
filename_dendrogram = paste0('dendrogram', '_', software, '_', contrast_label, '.pdf')
path_figure = paste(path_output, filename_dendrogram, sep = '')
pdf(file = path_figure)

# ASSIGN COLORS TO SAMPLE LABELS
# get the sample labels used in the dendrogram in order
labels_dend = as.data.frame(labels(avg_col_dend))
colnames(labels_dend) = c("sample_labels")
# create a variable with sample colors
colors_label_dend = merge(labels_dend, 
                          conditions, 
                          by.x = "sample_labels", 
                          by.y = "name",
                          sort = FALSE)
# change the conditions to colors
colors_label_dend$condition = replace(
  colors_label_dend$condition,
  colors_label_dend$condition == contrast_factors$reference,
  "#149E8B")
colors_label_dend$condition = replace(
  colors_label_dend$condition,
  colors_label_dend$condition == contrast_factors$condition,
  "#B53A5B")
  
# create a color vector
colors_label_dend = as.character(colors_label_dend$condition)

# assign the color vector to the dendrogram
labels_colors(avg_col_dend) <- colors_label_dend

# calculate margins
nlabels_dend = nrow(cts)
mar_C = nlabels_dend/25             # UP
mar_A = nlabels_dend/25             # DOWN
mar_B = nlabels_dend/25             # LEFT
mar_D = nlabels_dend/800            # RIGHT
mar_E = (nlabels_dend * 2)/100      # SCALE
mar_F = 100/((nlabels_dend*2)+5)    # LABEL SIZE

avg_col_dend = set(avg_col_dend, "labels_cex", mar_F)

par(mar=c(mar_A, mar_B, mar_C, mar_D) + mar_E)
plot(avg_col_dend,
     main = title,
     cex.main = 1.5, #title
     ylab="Euclidean distance",
     sub=NULL)

silence <- dev.off()

#-------------------------------------------------------------------------------
#                         EXPORT HEATMAP
#-------------------------------------------------------------------------------
# see https://www.rdocumentation.org/packages/gplots/versions/3.1.3/topics/heatmap.2
ptm('Building Heatmap')
filename_heatmap = paste0('heatmap', '_', software, '_', contrast_label, '.pdf')
path_figure = paste(path_output, filename_heatmap, sep = '')
pdf(file = path_figure)
heat_map = heatmap.2(x, 
                     distfun = function(x) dist(x, method="euclidean"),
                     hclustfun = function(x) hclust(x, method="ward.D2"),
                     col=pal,
                     main=title,
                     trace='none',
                     scale='none',
                     symm=F, #to avoid symmetrical x values
                     symkey=F,
                     symbreaks=F,
                     margins=c(margins, margins),
                     key.title = 'Expression',
                     key.xlab = 'Fold change',
                     key.ylab = 'Counts')

silence <- dev.off()