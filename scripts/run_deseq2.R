suppressMessages(library("DESeq2"))

#Command line arguments that will be used from pipeline.xml
args = commandArgs(trailingOnly = TRUE)

print(paste0("Input counts file: ", as.character(args[1])))
print(paste0("Input annotation file: ", as.character(args[2])))
print(paste0("Reference factor: ", as.character(args[3])))
print(paste0("Output directory: ", as.character(args[4])))

#The path to the featureCounts file
countFile=as.character(args[1])

#The path to the condition file (that produced by the pipeline)
annotationFile=as.character(args[2])

referenceFactor=paste(as.character(args[3]))

#The output directory to save the analysis results
outputDir=as.character(args[4])

#Loads the read counts of featureCounts
cts = read.csv(file = countFile,
                    sep = "\t",
                    skip = 1,
                    row.names = 1,
                    header = TRUE,
                    check.names=FALSE)

#Keeps only the gene_ids and the quantification results of each sample
cts = cts[-(1:5)]

#Loads the annotations of the samples
coldata = read.csv(file = annotationFile, 
                   sep = "\t",
                   row.names = 1,
                   header = TRUE)

#Set the factors for the comparison 
coldata = coldata[,c("condition","label")]
coldata$condition = factor(coldata$condition)

#Establish the main factor
coldata$condition = relevel(coldata$condition, referenceFactor)
coldata$label = factor(coldata$label)

# Checking that the cts and annotations are correct
test_col = all(rownames(coldata) %in% colnames(cts))
print(paste0("Rownames of 'conditions' are the same as colnames of 'counts': ", test_col))

#Test if colnames of cts are in the same order that rownames of coldata
if (all(rownames(coldata) == colnames(cts)) == FALSE) {
  cts = cts[, rownames(coldata)] #if not, arranges them
}

#Construct DESeq2 dataset
dds = DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

#Keep only rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Performs the differential expression analysis
dds <- DESeq(dds)

#Visualize the results
res <- results(dds)
summary(res)

#Order results by the smallest p-value
resOrdered <- res[order(res$pvalue),]

#Save results as tsv file
setwd(outputDir)
write.csv(as.data.frame(resOrdered),
          file=paste(coldata$condition, referenceFactor, sep="_"))
