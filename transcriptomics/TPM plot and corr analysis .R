library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(gridExtra)
library(TranscripTools)
library(writexl)

#load expression data
data <- as.data.frame(
  read_xlsx("TPMs_RS_CLL.xlsx", sheet = 2)
)

#remove first non-expression column
data <- data[, -1, drop = F]

#use gene symbols as row names
rownames(data) <- make.names(data$symbol, unique = T)

#remove gene-symbol column so expression matrix contains numeric values only
data <- data[, colnames(data) != "symbol", drop = F]

#load metadata
meta <- as.data.frame(
  read_xlsx("TPMs_RS_CLL.xlsx", sheet = 1)
)

#convert Class to factor
meta$Class <- as.factor(meta$Class)

#remove samples where Class is the literal string "NA"
meta <- subset(meta, Class != "NA")

#remove Nadeu's samples
meta <- subset(
  meta,
  !Class %in% c(
    "Richter (Nadeu): fase CLL",
    "Richter (Nadeu): fase RS"
  )
)

#drop unused factor level
meta$Class <- droplevels(meta$Class)

#check which metadata samples are absent from the expression matrix
missing_in_data <- setdiff(meta$ID, colnames(data))
print(missing_in_data)

#identify samples present in both metadata and expression matrix
common_samples <- intersect(colnames(data), meta$ID)

#keep only expression samples that are also present in metadata
data <- data[, common_samples, drop = F]

#keep only metadata samples that actually have expression data
meta <- meta[meta$ID %in% common_samples, , drop = F]

#put metadata in exactly the same order as expression columns
meta <- meta[match(colnames(data), meta$ID), , drop = F]

#sanity checks
stopifnot(
  ncol(data) == nrow(meta),
  !anyNA(meta$ID),
  identical(meta$ID, colnames(data))
)

#plot
genes <- c("ATR",
           "PRKCA",
           "PRKACA",
           "CAMK4",
           "PRKCB",
           "LMNB1",
           "LMNB2",
           "LMNA",
           "DENND4A",
           "TCAF2",
           "LTBP2",
           "AMPD3",
           "WDR44",
           "SH3PXD2A",
           "MYH3",
           "TRPM8",
           "TRPM2",
           "TRPM7") #old one remains for reference
genes <- c("CD274","EPHA2","EPHA4","EPHA10") #new one
whiskyplot(
  data,
  genes,
  meta,
  sample_col = "ID",
  group_col = "Class"
)

##"raw" visualization
#select only genes of interest before converting to long format
tpm <- data[genes[genes %in% rownames(data)], , drop = F] %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "ID",
    values_to = "TPM"
  ) %>%
  left_join(
    meta[, c("ID", "Class")],
    by = "ID"
  )

#plot raw TPM values without statistical tests
ggplot(
  tpm,
  aes(x = Class, y = TPM, fill = Class)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.2,
    size = 1,
    alpha = 0.6
  ) +
  facet_wrap(
    ~Gene,
    scales = "free_y"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Gene expression",
    x = "Class",
    y = "TPM",
    fill = NULL
  )

##save excel file
#select genes of interest from the expression matrix
selected_data <- data[genes[genes %in% rownames(data)], , drop = F]

#restore gene symbols as a regular column for the Excel output
selected_data <- data.frame(
  symbol = rownames(selected_data),
  selected_data,
  check.names = F
)

#save metadata and raw TPM values for selected genes in the same Excel workbook
write_xlsx(
  list(
    metadata = meta,
    expression_TPM = selected_data
  ),
  path = "selected_genes_TPM.xlsx"
)

##PD-L1 vs EPHA correlation analysis
#keep only CLL samples
cll_meta <- meta %>%
    filter(Class == "CLL")

#check again how many CLL samples
nrow(cll_meta)

#extract CD274 and EPHA expression only from CLL samples
cor_data <- data[genes, cll_meta$ID, drop = F] %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "ID", values_to = "TPM") %>%
  pivot_wider(names_from = Gene, values_from = TPM)

#inspect the data used for the correlations
head(cor_data)

#calculate Spearman correlations between CD274 and each EPHA gene
cor_results <- lapply(
  c("EPHA2", "EPHA4", "EPHA10"),
  
  function(gene) {
    test <- cor.test(cor_data$CD274, cor_data[[gene]], method = "spearman", exact = F)
    data.frame(
      Gene = gene,
      n = nrow(cor_data),
      rho = unname(test$estimate),
      p_value = test$p.value)
  }
) %>%
  bind_rows()

#adjust the three p-values for multiple testing
cor_results <- cor_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

#print correlation results
cor_results

##visualize correlations
#convert the correlation data to long format for plotting
cor_plot <- cor_data %>% 
  select(ID,CD274,EPHA2,EPHA4,EPHA10) %>%
  pivot_longer(
    cols = c(EPHA2, EPHA4, EPHA10),
    names_to = "Gene",
    values_to = "EPHA_TPM") %>%
  mutate(
    CD274_logTPM = log2(CD274 + 1),
    EPHA_logTPM = log2(EPHA_TPM + 1))

#prepare correlation labels for each panel
cor_labels <- cor_results %>%
  mutate(
    label = paste0(
      "rho = ", round(rho, 2),
      "\np = ", format.pval(p_value, digits = 2, eps = 0.001),
      "\nFDR = ", format.pval(p_adj, digits = 2, eps = 0.001),
      "\nn = ", n))

#plot PD-L1 expression against each EPHA gene
ggplot(
  cor_plot,
  aes(x = CD274_logTPM,
    y = EPHA_logTPM)) +
  geom_point(
    size = 2,
    alpha = 0.7) +
  geom_smooth(
    method = "lm",
    se = T,
    linewidth = 0.7) +
  geom_text(
    data = cor_labels,
    aes(x = -Inf,y = Inf,
      label = label),
    inherit.aes = F,
    hjust = -0.1,
    vjust = 1.1,
    size = 3.5) +
  facet_wrap(
    ~Gene,
    scales = "free_y") +
  theme_bw() +
  labs(title = "PD-L1 and EPHA expression in CLL",
    x = expression(CD274~log[2](TPM + 1)),
    y = expression(EPHA~log[2](TPM + 1)))

##save excel file
#select genes of interest from the expression matrix
selected_data <- data[genes[genes %in% rownames(data)], , drop = F]

#restore gene symbols as a regular column for the Excel output
selected_data <- data.frame(
  symbol = rownames(selected_data),
  selected_data,
  check.names = F)

#save metadata, expression values and correlation results
write_xlsx(
  list(
    metadata = meta,
    expression_TPM = selected_data,
    CLL_correlations = cor_results
    ),
  path = "PDL1_EPHA_CLL_correlations.xlsx"
)
