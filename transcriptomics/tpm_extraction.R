library(dplyr)
library(biomaRt)

##listing files
path <- "/Users/FV/Downloads/GSE282511_RAW/"
files <- list.files(path, pattern = "\\.sf\\.gz$", full.names = TRUE)

##extracting TPM data
#creating empty df
file_list <- list()

#looping for each file
for (file in files){
  f <- as.data.frame(read.table(file))
  df <- f[,c(1,4)] #extract ENST and TPM counts
  colnames(df) <- df[1,]
  df <- df[-1,] #delete colnames' row
  
  colnames(df)[2] <- gsub(pattern = "_quant.sf.gz$", replacement = "", basename(file))
  
  file_list[[basename(file)]] <- df #appending df to list
}

summary(file_list)

#merging all TPMs
tpm_data <- Reduce(function(x, y) merge(x, y, by = "Name"), file_list)

##convert gene ID
#connect to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#removing transcript version number
enst_ID <- tpm_data$Name %>%
  gsub(pattern = "\\.\\d+$", replacement = "")
tpm_data$ENST_ID <- enst_ID 

#getting gene names
anno <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_transcript_id",
  values = enst_ID,
  mart = ensembl
)

#merging gene annotations
tpm_data <- merge(tpm_data, anno, by.x = "ENST_ID", by.y = "ensembl_transcript_id") %>%
  relocate(c("ensembl_gene_id","external_gene_name"), .after = "Name")
head(tpm_data[,1:5])
colnames(tpm_data)[4] <- "Genes"

#trasforming in numeric
tpm_data[,-c(1:4)] <- lapply(tpm_data[,-c(1:4)], as.numeric)

#averaging transcript isoforms
tpm_data_fnl <- tpm_data %>%
  select(-ENST_ID, -Name) %>%
  group_by(Genes) %>%
  summarise(across(everything(), mean), .groups = "drop")
