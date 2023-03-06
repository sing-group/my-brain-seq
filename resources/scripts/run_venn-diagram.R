## Script input parameters:
##  1.- venn_path: path to a file with the DE features returned by DESeq2 + EdgeR, one per line.
##  2.- venn_output: path to the output directory for the graph.
##  3.- venn_output_format: file format of the output image (png, svg or tiff).
##  4.- input_contrast : a line of the contrast_file.

suppressMessages(library('VennDiagram'))

# to print without "[1]"
pt = function(text){
  cat(text, sep='\n')
}
# to print messages
ptm = function(text){
  header = '[PIPELINE -- venn -- run_venn-diagram.R > deseq + edger]: '
  cat(paste(header, text), sep='\n')
}

#INPUTS
args = commandArgs(trailingOnly = TRUE)
venn_path = as.character(args[1])
venn_output = as.character(args[2])
venn_output_format = as.character(args[3])
input_contrast = as.character(args[4])

ptm("======================================================")
ptm('        [PIPELINE -- venn]: run_venn-diagram.R        ')
ptm('......................................................')
ptm(paste0("  Input venn file:   ", as.character(args[1])))
ptm(paste0("  Input output dir:  ", as.character(args[2])))
ptm(paste0("  Image format:      ", as.character(args[3])))
ptm("======================================================")

# FUNCTIONS
# to get the reference and condition factors from the input_contrast (pipeline)
get_contrast_factors = function(input_contrast) {
  contrast_table = read.delim(text = paste('name\n', input_contrast), sep = '=')
  contrast = strsplit(as.character(contrast_table$name), '-')[[1]]
ref_factor = trimws(contrast[2])
  cond_factor = as.character(trimws(contrast[1]))
  result = list("reference" = ref_factor, "condition" = cond_factor)
}

#LOADING FILE
venn_table = read.delim(venn_path)

#VENN DIAGRAM
#venn variables
deseq = deseq[venn_table$deseq != ""]
edger = edger[venn_table$edger != ""]

#colors of the venn diagram
myCol = c('#459991', '#FF972F')

# building the name of the output file
factors = get_contrast_factors(input_contrast)
venn_contrast_label = paste0(factors$condition, '-', factors$reference)
output_filename = paste0('DESeq2-EdgeR_', venn_contrast_label, '_results_venn.', venn_output_format)

#build venn chart
setwd(venn_output)
ptm('Saving venn diagram')
silence <- venn.diagram(
  x = list(deseq, edger),
  category.names = c('DESeq2', 'EdgeR'),
  imagetype=venn_output_format,
  filename = output_filename,
  output=TRUE,
  scaled=FALSE,
  #image
  height = 1050, 
  width = 1050, 
  resolution = 550,
  margin = 0.38,
  #title
  main = 'Differential expression analysis',
  main.fontface = "bold",
  main.fontfamily = "sans",
  main.cex = 0.6,
  main.pos = c(0.5, 0.95),
  main.col = '#002464',
  #subtitle
  sub = venn_contrast_label,
  sub.fontface = NULL,
  sub.fontfamily = "sans",
  sub.cex = 0.5,
  sub.pos = c(0.5, 0.7),
  sub.col = '#636363',
  #circles
  lty = 'blank',
  fill = myCol,
  #numbers
  cex = 0.6,
  fontfamily = "sans",
  #set names
  cat.cex = 0.6,
  cat.just = list(c(0.1,0) , c(0.9,0)),
  ext.text = TRUE,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
)
ptm('Done')