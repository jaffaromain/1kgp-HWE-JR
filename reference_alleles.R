
# data:
# https://www.well.ox.ac.uk/~wrayner/strand/RefAlt.html
library(readr)
library(readxl)
library(tidyverse)
library(data.table)
# humanomni <- read_csv("~/Downloads/HumanOmni25-8v1-2_A1.csv", 
#                       skip = 7)
# homniM <- as.data.frame(read.table("~/Downloads/HumanOmni25M-8v1-1_B-b37.txt", header = F))
# 
# # head(humanomni)
# humanomni <- humanomni %>% select(Chr, MapInfo, Name)
# homniM <- homniM %>% select("Chr" = V2, "Name" = V1, "MapInfo" = V3)
# homniM$Chr[homniM$Chr == "X"] <- 23
# homniM$Chr[homniM$Chr == "Y"] <- 24
# homniM$Chr[homniM$Chr == "MT"] <- 26
# 
# 
# humanomni <- rbind(humanomni, homniM)
# snp_table <- as.data.frame(
#   read.table("~/Desktop/CHL5207 - Lab Work/snpTable.txt",   header = T, 
#              sep = "",   fill = TRUE,  quote = "", check.names = F)) 
# humanomni$GRCh37 <- humanomni$MapInfo
# humanomni$Chr <- as.numeric(humanomni$Chr)
# 
# full <- full_join(snp_table, humanomni, by = c("Chr", "GRCh37"))
# 
# full <- full %>% select(-c(CallRate, CallErr, MendErr, `P(HWE)`, FailedQC, cM,
#                            A1, A2, MAF, Genotypes)) %>% 
#   filter(!is.na(Chr))
################################################################################
omni1 <- as.data.frame(
  fread("~/Desktop/CHL5207 - Lab Work/HumanOmni25-8v1-2_A1.b37.RefAlt"))
colnames(omni1) <- c("SNP", "A1")


omni2 <- as.data.frame(
  fread("~/Desktop/CHL5207 - Lab Work/HumanOmni25M-8v1-1_B-b37.strand.RefAlt")) 
colnames(omni2) <- c("SNP", "A1")
  
full_omni <- full_join(omni1, omni2) %>% drop_na()
# fulldata <-  full_join(full_omni, full) %>% filter(!is.na(SNP))
fulldata <- distinct(full_omni,SNP, .keep_all= TRUE)

fwrite(fulldata, "reference_alleles.txt", col.names = F, quote = F, sep = "\t", row.names = F)


