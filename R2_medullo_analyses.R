###the data used here have been downloaded by https://hgserver1.amc.nl/cgi-bin/r2/main.cgi
###studies considered:Pfister mas5.0, Pfister fpkm, Gilbertson mas5.0
###probes used:  52731_at and 202431_s_at probes for AMBRA1 and c-MYC, respectively
###downloaded through the data grabber function of R2 for each subtype (group_3, group_4, shh, wnt)
###setwd("~/medulloblastoma/")

### libraries required

library(gtools)
library("tibble")
library(dplyr)

###

###read the data

Pshh.fpkm <- read.table("pfister_shh_fpkm_no.txt", sep = '\t', header = FALSE)
Pg4.fpkm <- read.table("pfister_g4_fpkm_no.txt", sep = '\t', header = FALSE)
Pg4.fpkm.df <-  as_data_frame(Pg4.fpkm)
P4g.fpkm.df.ok <-   select(my_data, V:Petal.Length)
Pg3.fpkm <- read.table("pfister_g3_fpkm_no.txt", sep = '\t', header = FALSE)
Pwnt.fpkm <- read.table("pfister_wnt_fpkm_no.txt", sep = '\t', header = FALSE)

Pshh.mas5 <- read.table("pfister_shh_mas5_no.txt", sep = '\t', header = FALSE)
Pg4.mas5 <- read.table("pfister_g4_mas5_no.txt", sep = '\t', header = FALSE)
Pg3.mas5 <- read.table("pfister_g3_mas5_no.txt", sep = '\t', header = FALSE)
Pwnt.mas5 <- read.table("pfister_wnt_mas5_no.txt", sep = '\t', header = FALSE)

Gshh.mas5 <- read.table("gilbertson_shh_mas5_no.txt", sep = '\t', header = FALSE)
Gg3.mas5 <- read.table("gilbertson_g3_mas5_no.txt", sep = '\t', header = FALSE)
Gg4.mas5 <- read.table("gilbertson_g4_mas5_no.txt", sep = '\t', header = FALSE)
Gwnt.mas5 <- read.table("gilbertson_wnt_mas5_no.txt", sep = '\t', header = FALSE)

#check group size

dim(Pshh.fpkm)
#2 48
dim(Pg4.fpkm)
#2 66
dim(Pwnt.fpkm)
#2 18    
dim(Pshh.mas5)
#2 61
dim(Pg4.mas5)
# 2 93
dim(Pg3.mas5)
#2 58
dim(Pwnt.mas5)
#2 19
dim(Gshh.mas5) 
# 
Gg3.mas5 <- read.table("gilbertson_g3_mas5_no.txt", sep = '\t', header = FALSE)
    Gg4.mas5 <- read.table("gilbertson_g4_mas5_no.txt", sep = '\t', header = FALSE)
    Gwnt.mas5 <- read.table("gilbertson_wnt_mas5_no.txt", sep = '\t', header = FALSE)


###vectors of data###
Pfpkm<-list()
Pmas5<-list()
Gmas5<-list()

Pfpkm[['shh']]<-Pshh.fpkm[,3:48]
Pfpkm[['g4']]<-Pg4.fpkm[,3:66]
Pfpkm[['g3']]<-Pg3.fpkm[,3:43]
Pfpkm[['wnt']]<-Pwnt.fpkm[,3:18]

Pmas5[['shh']]<-Pshh.mas5[,3:61]
Pmas5[['g4']]<-Pg4.mas5[,3:93]
Pmas5[['g3']]<-Pg3.mas5[,3:58]
Pmas5[['wnt']]<-Pwnt.mas5[,3:19]

Gmas5[['shh']]<-Gshh.mas5[,3:12]
Gmas5[['g4']]<-Gg4.mas5[,3:41]
Gmas5[['g3']]<-Gg3.mas5[,3:18]
Gmas5[['wnt']]<-Gwnt.mas5[,3:10]


###Pariwise Statistical Tests between the 4 subtypes####

Ttests.Pfpkm<-list()
Ttests.Pmas5<-list()
Ttests.Gmas5<-list()

subtypes<-c('shh', 'g4', 'g3', 'wnt')
d1<-combn(subtypes, 2)

pairwise.combs<-combinations(v=subtypes, n=4, r=2)

for(i in 1:nrow(pairwise.combs)){

  type1<-pairwise.combs[i, 1]
  type2<-pairwise.combs[i, 2]

  el.name<-paste0(type1,"VS",type2)
  print(el.name)
  
  Ttests.Pfpkm[[el.name]]<-t.test(as.vector(Pfpkm[[type1]][1,]),as.vector(Pfpkm[[type2]][1,]))
  Ttests.Pmas5[[el.name]]<-t.test(as.vector(Pmas5[[type1]][1,]),as.vector(Pmas5[[type2]][1,]))
  Ttests.Gmas5[[el.name]]<-t.test(as.vector(Gmas5[[type1]][1,]),as.vector(Gmas5[[type2]][1,]))

}



Ttests.Gmas5$shhVSwnt



Pg3.fpkm[1,3:43]

dim(Pg3.fpkm)
dim(Psonic.fpkm)
