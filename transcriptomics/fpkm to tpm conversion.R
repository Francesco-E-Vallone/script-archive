library(tidyverse)

##loading data
oci_ly3_fpkm <- read.delim("OCI-LY3/dataset #1/fpkm_genenames_oci_ly3.txt", header = T)

##FPKM conversion
#function to convert FPKM to TPM
convert_fpkm_to_tpm <- function(fpkm_matrix) {
  #extract FPKM values
  fpkm_values <- fpkm_matrix[, -1] #remove gene column
  
  #calculate library sizes
  library_sizes <- colSums(fpkm_values)
  
  #calculate scaling factor
  scaling_factor <- library_sizes / 1e6
  
  #calculate TPM values
  tpm_values <- fpkm_values / scaling_factor
  
  #add Gene column to the TPM matrix
  tpm_matrix <- cbind(fpkm_matrix[, 1, drop = FALSE], tpm_values)
  
  return(tpm_matrix)
}

#converting FPKM to TPM
tpm_data <- convert_fpkm_to_tpm(oci_ly3_fpkm)

##saving file
write.csv(tpm_data, file = "tpm-converted data_oci-ly3.csv")
