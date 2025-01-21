library(dplyr)
library(readxl)
library(rtracklayer)
library(biomaRt)

##data import
#importing counts data
counts <- read_xlsx("raw counts.xlsx", sheet = 2)

#importing metadata
meta <- read_xlsx("raw counts.xlsx", sheet = 1)
meta <- as.data.frame(meta[,1:4])
meta$`Core ID` <- paste0("S", meta$`Core ID`) #pasting "S" to match counts colnames
head(meta)

#importing gene lengths data
hg38_data <- rtracklayer::import("gencode.v25.annotation.gtf.gz")

##Ensembl ID to gene symbol conversion
ensg <- counts$Geneid
ensg <- sub("\\.\\d+$", "", ensg) #removing version numb

#connecting to the ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#retrieving gene symbols
gene_symb <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",                    
                   values = ensg,                            
                   mart = mart)

#appending to dataframe
counts <- data.frame(ensg = ensg, counts)
colnames(gene_symb)[1] <- "ensg"
counts <- merge(gene_symb, counts, by = "ensg")
head(counts)

#maintaining relevant columns
counts_ok <- counts[,-c(1,4:7,9:32)]
colnames(counts_ok)

##TPM conversion
#function to calculate TPM for a given counts column
calculate_tpm <- function(counts_col, lengths_col) {
  rpk <- counts_col / (lengths_col / 1000) #calculate rpk
  scaling_factor <- sum(rpk) / 1e6 #scaling factor
  tpm <- rpk / scaling_factor #calculate tpm
  return(tpm)
}

#list of counts columns in the merged_data (excluding gene_symbol and gene_length)
count_columns <- setdiff(names(counts_ok), c("hgnc_symbol", "Geneid","Length"))

#applying TPM calculation for each counts column
tpm_data <- counts_ok[, c("hgnc_symbol", "Geneid", "Length")] 
for (col in count_columns) {
  tpm_col_name <- paste0("TPM_", col)
  tpm_data[[tpm_col_name]] <- calculate_tpm(counts_ok[[col]], counts_ok$Length)
}

#checking data
print(head(tpm_data))

#linking metadata and TPM data
list_all <- list("metadata" = as.data.frame(meta),
                 "TPM data" = as.data.frame(tpm_data))

#writing xlsx
writexl::write_xlsx(list_all, "TPM_data.xlsx")
