library(tidyverse)
library(readr)
library(cmapR)

#setting wd
wd <- setwd("~/Scienze Mediche Dropbox/Pers.immuno/fe.vallone/collabs/chemotaxis & cell adhesion/RNA-seq/data/gsea/") #set in the same folder, otherwise returns an error

##loading files 
files <- list.files(path = wd, 
                    pattern = "*.gmt")

#creating empty list
gmt_list <- list()

#looping for each gmt file
for (file in files){
  #read each gmt file
  gmt <- readLines(file) #read.delim is not ideal
  gmt <- strsplit(gmt, "\t")
  
  #appending to nested list (same name as the file)
  gmt_list[[file]] <- gmt[[1]][-c(1,2)]
}

head(gmt_list)

##writing file
output_file <- "adhesion_and_chemotaxis_sig.gmt"

#open a connection (write line by line)
connect <- file(output_file, open = "w") 

for (file in names(gmt_list)){
  #extracting dataset name (without .gmt)
  db_name <- sub("\\.gmt", "", file)
  
  #extracting genes for each dataset
  genes <- as.character(unlist(gmt_list[[file]], use.names = FALSE))

  #creating a single gmt line
  gmt_line <- paste(c(db_name, genes), collapse = "\t")
  
  #append the lines to the file
  writeLines(gmt_line, con = connect)
}

#closing file connection (end of writing)
close(connect)

readLines("adhesion_and_chemotaxis_sig.gmt")
