###the data used here have been downloaded by https://hgserver1.amc.nl/cgi-bin/r2/main.cgi
###studies considered:Pfister mas5.0, Pfister fpkm, Gilbertson mas5.0
###downloaded through the data grabber function of R2 for each subtype (group_3, group_4, shh, wnt) using the HugoOnce option
###setwd("~/ambra_cmyc_medulloblastomaR2/")

### libraries required

library(gtools)
library(dplyr)
library(ggpubr)
library("CEMiTool")

###

#### Pfister fpkm databasets

#### read the data and check group sizes ####


Pg3.fpkm <- read.table("downloaded_data/Pfister_fpkm/pfister_g3_fpkm_all_log2.txt", sep = '\t', header = FALSE)
View(Pg3.fpkm)
###check the group size - Pshh.fpkm should contain 41 samples
dim(Pg3.fpkm)
# 43601 43
Pg3.fpkm <- as.data.frame(Pg3.fpkm)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pg3.fpkm <- select(Pg3.fpkm, V2:V43)
Pg3.fpkm <- t(Pg3.fpkm)
dim(Pg3.fpkm)

### coexpression analyses using CEMiTools  ###
rownames(Pg3.fpkm) <- Pg3.fpkm[,1]
Pg3.fpkm <- select(Pg3.fpkm, V3:V43)
#CEMiTool with pearson correlation
cem.Pg3.fpkm <- cemitool(Pg3.fpkm, diss_thresh = 0.6)
cem.Pg3.fpkm
filter_expr(cem.Pg3.fpkm)
#12 modules
#4245 genes
modules.Pg3.fpkm<- module_genes(cem.Pg3.fpkm)
#MYC in module 4 - AMBRA not found

#### Gilbertson mas5 databasets

#### read the data and check group sizes ####

Gg3.mas5 <- read.table("downloaded_data/Gilbertson_mas5.0/gilbertson_g3_mas5_all_log2.txt", sep = '\t', header = FALSE)
View(Gg3.mas5)
###check the group size - Pshh.fpkm should contain 16 samples
dim(Gg3.mas5)
# 2 18
Gg3.mas5 <- as.data.frame(Gg3.mas5)
rownames(Gg3.mas5) <- Gg3.mas5[,1]
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Gg3.mas5 <- select(Gg3.mas5, V3:V18)
Gg3.mas5 <- t(Gg3.mas5)
dim(Gg3.mas5)


### coexpression analyses using CEMiTools  ###

#CEMiTool with pearson correlation
cem.Gg3.mas5 <- cemitool(Gg3.mas5, diss_thresh = 0.6)
cem.Gg3.mas5
filter_expr(cem.Gg3.mas5)
#9 modules
#1542 genes
modules.Gg3.mas5<- module_genes(cem.Gg3.mas5)
#MYC in module 5 - AMBRA not found
