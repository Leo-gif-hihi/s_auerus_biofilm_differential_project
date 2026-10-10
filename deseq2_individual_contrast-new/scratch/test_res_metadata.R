library(DESeq2)
dds <- makeExampleDESeqDataSet(n=10, m=12)
dds$condition <- factor(rep(c("A", "B", "C"), each=4))
dds <- DESeq(dds, betaPrior=FALSE)

res <- results(dds, contrast=c("condition", "C", "B"))
print(str(res))
# Is the contrast vector stored in the metadata?
print(metadata(res))
