## Script input parameters:
##  1.- path_counts: path to the file in tsv with the deseq2/edger results
##  2.- path_output: the output path of the image.
##  3.- comparison_label: this text will be shown in the volcano title
##                        and in the output filename.
##  4.- software: the software that produced the tsv file (deseq/edger).

#-------------------------------------------------------------------------------
#                               INPUT
#-------------------------------------------------------------------------------
suppressMessages(library('EnhancedVolcano'))

args = commandArgs(trailingOnly = TRUE)
path_counts = as.character(args[1])
path_output = as.character(args[2])
comparison_label = as.character(args[3])
software = as.character(args[4])

# to print bare text without "[1]"
pt = function(text){
  cat(text, sep='\n')
}
# to print messages
ptm = function(text, sft = software){
  header = '[PIPELINE -- volcano -- run_enhanced-volcano.R > '
  header = paste0(header, sft, ']: ')
  cat(paste(header, text), sep='\n')
}

ptm("======================================================")
ptm('     [PIPELINE -- volcano]: run_enhanced-volcano.R    ')
ptm('......................................................')
ptm(paste0("  Input counts file: ", as.character(args[1])))
ptm(paste0("  Input output dir:  ", as.character(args[2])))
ptm(paste0("  Title:             ", as.character(args[3])))
ptm(paste0("  Software:          ", as.character(args[4])))
ptm("======================================================")



#READING
res = read.csv(file = path_counts, 
               sep = "\t",
               row.names = 1,
               header = TRUE)

#GET ONLY ROWS WITH QVALUE != NA
res_filtered = res[complete.cases(res$qvalue),]

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

#-------------------------------------------------------------------------------
#                             VOLCANO PLOT
#-------------------------------------------------------------------------------
ptm('Building Volcano plot')
volcano_name = paste0(comparison_label, '_volcano', '.pdf')
vol_path = paste0(path_output, volcano_name)
pdf(file=vol_path)

# AXIS LIMITS
# This section calculates the limits used to plot the volcano. There are minimum
# values for those limits (default values). These values correspond to the 
# thresholds used to plot the data (qvalue, fold-change). If no data points reach
# the minimum values, then the default values will be used.

# default values
qval = 0.05
logfc = 0.5
marg = 1.1

# default limits 
xlim_low_def = -logfc; xlim_high_def = logfc
ylim_high_def = -log10(qval)

# automatic limits
xlim_low  = min(res_filtered$log2FC) * marg
xlim_high = max(res_filtered$log2FC) * marg
ylim_high = max(-log10(res_filtered$qvalue)) * marg

# if default limits are not in automatic limits then use default limit
if (xlim_low > xlim_low_def){xlim_low = xlim_low_def}
if (xlim_high < xlim_high_def){xlim_high = xlim_high_def}
if (ylim_high < ylim_high_def){ylim_high = ylim_high_def}

# set xlim and ylim variables
xlim = c(xlim_low, xlim_high)
ylim = c(0, ylim_high)

EnhancedVolcano(res_filtered,
                x = vol_x, y = vol_y,
                xlim = xlim, ylim = ylim,
                title = comparison_label,
                subtitle = titl,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                col=c('#8D8E94', '#149E8B', '#106392', '#A62C7A'),
                colAlpha = 0.55,
                lab = rownames(res_filtered),
                labSize = 4)

silence <- dev.off()
ptm('Done')