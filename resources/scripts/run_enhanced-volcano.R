## Script input parameters:
##  1.- path_counts: path to the file in tsv with the deseq2/edger results
##  2.- path_output: the output path of the image.
##  3.- comparison_label: this text will be shown in the volcano title
##                        and in the output filename.
##  4.- software: the software that produced the tsv file (deseq/edger).

suppressMessages(library('EnhancedVolcano'))

#INPUT
args = commandArgs(trailingOnly = TRUE)
print("======================================================")
print('       [PIPELINE -- dea-integration]: Rscript         ')
print('......................................................')
print(paste0("  Input counts file: ", as.character(args[1])))
print(paste0("  Input output dir:  ", as.character(args[2])))
print(paste0("  Title:             ", as.character(args[3])))
print(paste0("  Software:          ", as.character(args[4])))
print('......................................................')

path_counts = as.character(args[1])
path_output = as.character(args[2])
comparison_label = as.character(args[3])
software = as.character(args[4])

#READING
res = read.csv(file = path_counts, 
               sep = "\t",
               row.names = 1,
               header = TRUE)

#CREATING VARIABLES
vol_x = 'log2FoldChange'
vol_y = 'padj'
sfw_name = 'DESeq2: '
if (software == 'edger') {
  vol_x = 'logFC'
  vol_y = 'PValue'
  sfw_name = 'EdgeR: '
} else if (software == 'both') {
  vol_x = 'log2FC'
  vol_y = 'q.value'
  sfw_name = 'DESeq2-EdgeR: '
}
titl = paste0(sfw_name, 'Differential expression')

#PLOT: Volcano plot
print('  Building Volcano plot')
volcano_name = paste0(comparison_label, '_volcano', '.pdf')
vol_path = paste0(path_output, volcano_name)
pdf(file=vol_path)
EnhancedVolcano(res,
                title = comparison_label,
                subtitle = titl,
                pCutoff = 10e-5,
                FCcutoff = 1,
                lab = rownames(res),
                x = vol_x,
                y = vol_y)
dev.off()
print("======================================================")