library("DESeq2")

countFile="/home/dannyzimm/Documentos/mirna-pipeline/output/3_feature-counts/all-counts.txt"
annotationFile="/home/dannyzimm/Documentos/mirna-pipeline/data/conditions.txt"
referenceFactor="control" #usually control

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
all(rownames(coldata) %in% colnames(cts))

#Test if colnames of cts are in the same order that rownames of coldata
if (all(rownames(coldata) == colnames(cts)) == FALSE) {
  cts <- cts[, rownames(coldata)] #if not, arranges them
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
res

