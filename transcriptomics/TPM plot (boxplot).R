library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(gridExtra)

#load and process data
data <- as.data.frame(read_xlsx("TPMs_RS PDX.xlsx"))
rownames(data) <- data$Genes
data <- data[,-c(1,2)]

#create metadata
meta <- data.frame(Samples = colnames(rxra),
                   type = c("RS1316","RS1316","RS1316","RS1316",
                            "RS1316", "RS1316", "RS9737", "RS9737",
                            "RS9737","RS9737", "IP867","IP867",
                            "IP867", "IP867", "RS1050", "RS1050",
                            "RS1050"))

#plot function
plot_TPM <- function(gene){
  #select gene
  gene_data <- data[rownames(data) == gene, ] 
  
  #build data.frame
  gene_data_df <- data.frame(TPM = as.numeric(gene_data), 
                             type = meta$type[match(colnames(gene_data), meta$Samples)])
  
  #plot
  plot <- ggplot(gene_data_df, aes(x = type, y = TPM, fill = type)) +
    geom_boxplot(outlier.shape = NA, color = "black", lwd = 0.5) +  # Frame the boxplots
    geom_jitter(position = position_jitter(width = 0.2), size = 2, color = "black") +  # Add jittered dots for individual samples
    labs(x = "Sample Type", y = "TPM", title = paste0(gene, " Expression")) +
    theme_minimal() +
    theme(
      legend.position = "none",  # Remove legend (since fill already shows the groups)
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Frame the panel
    )
  
  #return the plot
  print(plot)
}

genes_of_interest <- c("RXRA", "NR1H3", "NR1H2", "PPARG", "RARA", "RARB", 
                       "NCOA2", "NCOA1", "PPARA", "NR1H4", "NR1I2")
genes_of_interest2 <- c("RARA", "RARB", "RARG", "PPARG", "NCOA1", "NCOA2", "NCOA3", 
                        "NCOR1", "NCOR2", "HDAC1", "HDAC2", "HDAC3","HDAC4")
#plot recursively
for (gene in genes_of_interest2){
  plot_2 <- plot_TPM(gene)
}
