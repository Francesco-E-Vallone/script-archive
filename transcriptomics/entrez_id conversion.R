library(data.table)
library(dplyr)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

##loading data
tpm <- as.data.frame(read.delim("TPMs.txt", header = T))
rownames(tpm) <- tpm$Genes

meta <- read.delim("metadata.txt", header=T)

##ID conversion
hs <- org.Hs.eg.db
my_symbols <- tpm$Genes

#selecting IDs
entrez_id <- select(hs,
                    keys = my_symbols,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")

#merging IDs and gene symbols
new_tpm <- merge(entrez_id, tpm, by.x= "SYMBOL", by.y="Genes")
new_tpm <- new_tpm[!is.na(new_tpm$ENTREZID),]
new_tpm <- new_tpm[,-1]

##output
write.csv(new_tpm, "tpm.csv")
write.csv(meta, "metadata.csv")
