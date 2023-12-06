
library("DESeq2")

# data <- read.delim('filtered_transcriptomic.txt', header = TRUE, sep = "\t")
data <- read.delim('proteome.txt', header = TRUE, sep = "\t")
metadata <- read.delim('metadata.txt', header = TRUE, sep = "\t")

# data1 <- data[,c(1,39:41,51:53)] # M, transcriptome
# data1 <- data[,c(1,45:47,51:53)] # L, transcriptome
# data1 <- data[,c(1,48:50,51:53)] # U, transcriptome
# data1 <- data[,c(2,41:43,53:55)] # M, proteome
# data1 <- data[,c(2,47:49,53:55)] # L, proteome
data1 <- data[,c(2,50:52,53:55)] # U, proteome
data1[,2:7] <- round(data1[,2:7],digits=0) # for proteome only
data2 <- metadata[c(40:42,46:48),]

dds <- DESeqDataSetFromMatrix(countData=data1, 
                              colData=data2, 
                              design=~id, tidy = TRUE)

dds <- DESeq(dds)

res <- results(dds)

result <- cbind(res@rownames,res$log2FoldChange,res$padj)

# head(results(dds, tidy=TRUE))

write.table(result, 'UvsP_proteome.txt', quote=FALSE, row.names=FALSE, sep="\t")

result_df <- as.data.frame(result)
result_df$V2 <- as.numeric(result_df$V2)
result_df$V3 <- as.numeric(result_df$V3)
result_df_up <- result_df[result_df[,2]>0.2,]
# result_sorted <- result_df[order(result_df$V2),]result_df
write.table(result_df_up, 'UvsP_upreg_proteome.txt', quote=FALSE, row.names=FALSE, sep="\t")

result_df_down <- result_df[result_df[,2]<0,]
result_df_down <- result_df_down[abs(result_df_down$V2)>0.2,]
write.table(result_df_down, 'UvsP_downreg_proteome.txt', quote=FALSE, row.names=FALSE, sep="\t")


