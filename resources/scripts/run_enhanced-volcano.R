## Script input parameters:
##  1.- path_counts: path to the file in tsv with the deseq2/edger results
##  2.- path_output: the output path of the image.
##  3.- comparison_label: this text will be shown in the volcano title
##                        and in the output filename.
##  4.- software: the software that produced the tsv file (deseq/edger).

suppressMessages(library('EnhancedVolcano'))

# to print bare text without "[1]"
pt = function(text){
  cat(text, sep='\n')
}
# to print messages
ptm = function(text){
  header = '[PIPELINE -- volcano -- run_enhanced-volcano.R]:'
  cat(paste(header, text), sep='\n')
}

#INPUT
args = commandArgs(trailingOnly = TRUE)

pt(''); pt('')
pt("======================================================")
pt('     [PIPELINE -- volcano]: run_enhanced-volcano.R    ')
pt('......................................................')
pt(paste0("  Input counts file: ", as.character(args[1])))
pt(paste0("  Input output dir:  ", as.character(args[2])))
pt(paste0("  Title:             ", as.character(args[3])))
pt(paste0("  Software:          ", as.character(args[4])))
pt("======================================================")
pt(''); pt('')

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
vol_x = 'log2FC'
vol_y = 'qvalue'
sfw_name = 'DESeq2: '
if (software == 'edger') {
  sfw_name = 'EdgeR: '
} else if (software == 'both') {
  sfw_name = 'DESeq2-EdgeR: '
}
titl = paste0(sfw_name, 'Differential expression')

#PLOT: Volcano plot
ptm('Building Volcano plot')
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
silence <- dev.off()
ptm('Done')