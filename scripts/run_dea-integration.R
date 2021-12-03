## Script input parameters:
##  1.- path_results: path to the file with the results integrated with
##                    repetitions (DEmiRNAs_deseq-edger.tsv).
##  2.- path_output: path to the directory used to save the results.

suppressMessages(library('data.table'))

#INPUTS
args = commandArgs(trailingOnly = TRUE)
path_results = as.character(args[1])
path_output = as.character(args[2])

#READING
res = read.csv(file = path_results,
               sep = "\t",
               header = TRUE)

ordered = res[order(res$Feature),]

#ESTIMATION OF THE AVERAGE VALUE
#number of different features
n_feat = length(levels(as.factor(res$Feature)))
#empty matrix to build results
  #for the volcano
data = matrix(NA, n_feat, 3)
colnames(data) <- colnames(res)

  #for summarize qvals
ave_QVAL = matrix(NA, n_feat, 4)
colnames(ave_QVAL) = c('Feature', 'qval:Deseq', 'qval:EdgeR', 'Mean')
  #for summarize log2FC
ave_FC = matrix(NA, n_feat, 4)
colnames(ave_FC) = c('Feature', 'log2FC:Deseq', 'log2FC:EdgeR', 'Mean')

last = ''
counter=1
#parses res getting each row index
for (i in 1:nrow(ordered)) {
  feat = as.character(ordered$Feature[i])
  fc = as.numeric(ordered$log2FC[i])
  qval = ordered$q.value[i]
  # if current value is a duplicate averages values
  if (feat == last) {
    row1 = c(fc, qval) # row i values
    row2 = c(ordered$log2FC[i-1], ordered$q.value[i-1]) #row i-1 values
    #means between DESeq and EdgeR values
    row12_FC = mean(c(row1[1], row2[1]))
    row12_QVAL = mean(c(row1[2], row2[2]))
    rowM = c(feat, row12_FC, row12_QVAL)
    #builds the matrix each iteration
    data[counter, ] = rowM
    ave_QVAL[counter, ] = c(feat, row2[2], row1[2], row12_QVAL)
    ave_FC[counter, ] = c(feat, row2[1], row1[1], row12_FC)
    
    # EdgeR:
    #   log2FC ---> row1[1]
    #   qval -----> row1[2]
    #
    # DESeq2:
    #   log2FC ---> row2[1]
    #   qval -----> row2[2]
  
    counter=counter+1
  }
  last = feat
}

#BUILD THE RESULTS
#integrated results
integrated_RESULTS = as.data.frame(data)
# integrated_RESULTS = cbind(Feature = as.character(data[ , 1]), integrated_RESULTS)
# integrated_RESULTS$Feature = as.character(data[ , 1])
# integrated_RESULTS$log2FC = as.numeric(data[ , 2])
# integrated_RESULTS$qval = as.numeric(data[ , 3])
  
#summary table
summary_TABLE = data.table()
summary_TABLE$Feature = as.character(ave_QVAL[ , 1])
summary_TABLE$qval_Deseq = as.numeric(ave_QVAL[ , 2])
summary_TABLE$qval_EdgeR = as.numeric(ave_QVAL[ , 3])
summary_TABLE$qval_Average = as.numeric(ave_QVAL[ , 4])
summary_TABLE$log2FC_Deseq = as.numeric(ave_FC[ , 2])
summary_TABLE$log2FC_EdgeR = as.numeric(ave_FC[ , 3])
summary_TABLE$log2FC_Average = as.numeric(ave_FC[ , 4])
summary_TABLE = as.data.frame(summary_TABLE)

print('[PIPELINE -- dea-integration]: Average values calculated: ')
print(integrated_RESULTS)
#print(summary_TABLE)

#SAVE RESULTS
print('Saving integrated results')
path_output_file = paste0(path_output, '/', 'DEmiRNAs_deseq-edger_integrated', '.tsv')
write.table(integrated_RESULTS, path_output_file, row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
print('Saving summary table')
path_output_file = paste0(path_output, '/', 'DEmiRNAs_deseq-edger_summary', '.tsv')
write.table(summary_TABLE, path_output_file, row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)












