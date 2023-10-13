#!/usr/local/bin/Rscript
options(gsubfn.engine = "R")

library(tidyverse)
library(data.table)
library(tidyverse)
library(httr)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(R453Plus1Toolbox)
library(sqldf)
library(jsonlite)
source("/storage1/fs1/bolton/Active/Users/DucTran/UkbbVar/Script/supportFunction.R", local = TRUE)

args = commandArgs(trailingOnly=TRUE)

message("Reading in the CSV file...")
df = fread(args[1], sep = ",", header = T, data.table = FALSE)
message("Finished reading file.")

# Protected files, rarely changes
{
  TRUNCATING <- "/storage1/fs1/bolton/Active/Protected/Annotation_Files/BB.truncating.more.than.1.tsv"
  BOLTON_BICK_VARS <- "/storage1/fs1/bolton/Active/Protected/Annotation_Files/bick.bolton.vars3.txt"
  MUT2_BICK <- "/storage1/fs1/bolton/Active/Protected/Annotation_Files/topmed.n2.mutation.c.p.txt"
  MUT2_KELLY <- "/storage1/fs1/bolton/Active/Protected/Annotation_Files/kelly.n2.mutation.c.p.txt"
  MATCHES2 <- "/storage1/fs1/bolton/Active/Protected/Annotation_Files/matches.2.c.p.txt"
  GENE_LIST <- "/storage1/fs1/bolton/Active/Protected/Annotation_Files/oncoKB_CGC_pd_table_disparity_KB_BW.csv"
}
ONCO_KB_AVAILABLE_GENE <- "/storage1/fs1/bolton/Active/Users/DucTran/UkbbVar/Data/oncoKbCancerGeneList.tsv"

# Load all support data
{
  supportData <- list()
  
  # Cosmic
  chrList <- paste0("chr", c(1:22, "X", "Y"))
  supportData[["cosmic"]] <- lapply(chrList, function(x) {
    cosmic <- fread(paste0("/storage1/fs1/bolton/Active/Protected/Annotation_Files/cosmic/CosmicMutantExport.final.minimal.", x, ".tsv"), header = T, sep = "\t", quote = "", drop = c(2, 4))
    colnames(cosmic) <- c("COSMIC_ID", "var_key", "CosmicCount", "heme_cosmic_count", "myeloid_cosmic_count", "Gene_HGVSp_VEP")
    cosmic
  })
  names(supportData[["cosmic"]]) <- chrList
  supportData[["fullCosmic"]] <- supportData[["cosmic"]] %>% do.call(what = rbind)
  supportData[["cosmic"]] <- NULL
  
  supportData[["topmed.mutation.2"]] <- fread(MUT2_BICK, sep = "\t", header = T, data.table = FALSE)
  supportData[["kelly.mutation.2"]] <- fread(MUT2_KELLY, sep = "\t", header = T, data.table = FALSE)
  supportData[["matches.2.c.p"]] <- fread(MATCHES2, sep = "\t", header = T, data.table = FALSE)
  supportData[["vars"]] <- fread(BOLTON_BICK_VARS, sep = "\t", header = T, data.table = FALSE)
  
  supportData[["vars.truncating"]] <- fread(TRUNCATING, sep = "\t", header = T, data.table = FALSE)
  supportData[["gene_list"]] <- fread(GENE_LIST, sep = ",", header = T, data.table = FALSE)
  supportData[["gene_list"]] <- supportData[["gene_list"]] %>% dplyr::filter(remove != 1)
  
  # Read updated list of genes available in oncoKB, including aliases
  # The file can be directly read from the API link but download once at the beginning for consistency
  tmp <- fread(ONCO_KB_AVAILABLE_GENE, data.table = FALSE) %>%
    dplyr::filter(`OncoKB Annotated` == "Yes")
  tmp <- c(tmp$`Hugo Symbol`, str_squish(unlist(str_split(tmp$`Gene Aliases`, ",")))) %>%
    unique() %>%
    .[. != ""]
  supportData[["oncoKbAvailGene"]] <- tmp
  
  supportData[["redisHost"]] <- "9999"
}
  
message("Calculating Complexity")
df <- annotateComplexity(df)
message("Done Complexity")

df$key = paste0(df$CHROM, ":", df$POS, ":", df$REF, ":", df$ALT) 

message("Run AnnotatePD")
annotation <- annotatePD(df, supportData)
df <- left_join(df, annotation, by = c("CHROM", "POS", "REF", "ALT"))
message("Finished AnnotatePD")

write.csv(df, args[2], row.names = FALSE)