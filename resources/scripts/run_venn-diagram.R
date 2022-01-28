## Script input parameters:
##  1.- venn_path: path to a file with the DE features returned by DESeq2 + EdgeR, one per line.
##  2.- venn_output: path to the output directory for the graph.
##  3.- venn_output_format: file format of the output image (png, svg or tiff).

suppressMessages(library('VennDiagram'))

#INPUTS
args = commandArgs(trailingOnly = TRUE)
print("======================================================")
print('             [PIPELINE -- venn]: Rscript              ')
print('......................................................')
print(paste0("  Input venn file:   ", as.character(args[1])))
print(paste0("  Input output dir:  ", as.character(args[2])))
print(paste0("  Image format:      ", as.character(args[3])))
print('......................................................')
venn_path = as.character(args[1])
venn_output = as.character(args[2])
venn_output_format = as.character(args[3])

#LOADING FILE
venn_table = read.delim(venn_path, header = FALSE, col.names = c('DESeq2','EdgeR'))

#VENN DIAGRAM
#venn variables
deseq = venn_table$DESeq2
edger = venn_table$EdgeR

#colors of the venn diagram
myCol = c('#459991', '#FF972F')

#build venn chart
setwd(venn_output)
print('  Saving Venn diagram')
venn.diagram(
  x = list(deseq, edger),
  category.names = colnames(venn_table),
  imagetype=venn_output_format,
  filename = paste0('DESeq2-EdgeR_results_venn.', venn_output_format),
  output=TRUE,
  scaled=FALSE,
  #image
  height = 1050, 
  width = 1050, 
  resolution = 550,
  margin = 0.38,
  #title
  main = 'Common results',
  main.fontface = "bold",
  main.fontfamily = "sans",
  main.cex = 0.9,
  main.pos = c(0.5, 0.95),
  main.col = '#002464',
  #subtitle
  sub = 'Differential expression analysis',
  sub.fontface = NULL,
  sub.fontfamily = "sans",
  sub.cex = 0.6,
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
print('  Done')
print("======================================================")