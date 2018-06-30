###the data used here have been downloaded by https://hgserver1.amc.nl/cgi-bin/r2/main.cgi
###studies considered:Pfister mas5.0, Pfister fpkm, Gilbertson mas5.0
###probes used:  52731_at and 202431_s_at probes for AMBRA1 and c-MYC, respectively
###downloaded through the data grabber function of R2 for each subtype (group_3, group_4, shh, wnt)
###setwd("~/ambra_cmyc_medulloblastomaR2/")

### libraries required

library(gtools)
library(dplyr)
library(ggpubr)

###

#### Pfister fpkm databasets

#### read the data and check group sizes ####

Pshh.fpkm <- read.table("downloaded_data/Pfister_fpkm/pfister_shh_fpkm_log2.txt", sep = '\t', header = FALSE)
View(Pshh.fpkm)
###check the group size - Pshh.fpkm should contain 46 samples
dim(Pshh.fpkm)
# 2 48
Pshh.fpkm <- as.data.frame(Pshh.fpkm)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pshh.fpkm <- select(Pshh.fpkm, V3:V48)
Pshh.fpkm <- t(Pshh.fpkm)
dim(Pshh.fpkm)

Pg3.fpkm <- read.table("downloaded_data/Pfister_fpkm/pfister_g3_fpkm_log2.txt", sep = '\t', header = FALSE)
View(Pg3.fpkm)
###check the group size - Pshh.fpkm should contain 41 samples
dim(Pg3.fpkm)
# 2 43
Pg3.fpkm <- as.data.frame(Pg3.fpkm)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pg3.fpkm <- select(Pg3.fpkm, V3:V43)
Pg3.fpkm <- t(Pg3.fpkm)
dim(Pg3.fpkm)


Pg4.fpkm <- read.table("downloaded_data/Pfister_fpkm/pfister_g4_fpkm_log2.txt", sep = '\t', header = FALSE)
View(Pg4.fpkm)
###check the group size - Pshh.fpkm should contain 64 samples
dim(Pg4.fpkm)
# 2 66
Pg4.fpkm <- as.data.frame(Pg4.fpkm)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pg4.fpkm <- select(Pg4.fpkm, V3:V66)
Pg4.fpkm <- t(Pg4.fpkm)
dim(Pg4.fpkm)

Pwnt.fpkm <- read.table("downloaded_data/Pfister_fpkm/pfister_wnt_fpkm_log2.txt", sep = '\t', header = FALSE)
View(Pwnt.fpkm)
###check the group size - Pwnt.fpkm should contain 16 samples
dim(Pwnt.fpkm)
# 2 18
Pwnt.fpkm <- as.data.frame(Pwnt.fpkm)
#to remove the duplicate column with gene names
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pwnt.fpkm <- select(Pwnt.fpkm, V3:V18)
Pwnt.fpkm <- t(Pwnt.fpkm)
dim(Pwnt.fpkm)
#par(mar=c(1,1,1,1))

### correlation between AMBRA and MYC ###


cor.test (Pshh.fpkm[,1], Pshh.fpkm[,2], method= "pearson")
cor.test (Pshh.fpkm[,1], Pshh.fpkm[,2], method= "kendall")
cor.test (Pshh.fpkm[,1], Pshh.fpkm[,2], method= "spearman")

plot(Pshh.fpkm[,1], Pshh.fpkm[,2])

colnames(Pshh.fpkm) <- c("AMBRA","MYC")
pdf("plots/Pshh_fpkm_correlations_heatmap.pdf")
heatmap(t(Pshh.fpkm[1:20,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pshh.fpkm[21:46,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

cor.test (Pg3.fpkm[,1], Pg3.fpkm[,2], method= "pearson")
cor.test (Pg3.fpkm[,1], Pg3.fpkm[,2], method= "kendall")
cor.test (Pg3.fpkm[,1], Pg3.fpkm[,2], method= "spearman")

plot(Pg3.fpkm[,1], Pg3.fpkm[,2])

colnames(Pg3.fpkm) <- c("AMBRA","MYC")
pdf("plots/Pg3_fpkm_correlations_heatmap.pdf")
heatmap(t(Pg3.fpkm[1:20,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pg3.fpkm[21:41,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

cor.test (Pg4.fpkm[,1], Pg4.fpkm[,2], method= "pearson") # 0.42, p-value = 0.000489
cor.test (Pg4.fpkm[,1], Pg4.fpkm[,2], method= "kendall") #0.30, p-value = 0.0003276
cor.test (Pg4.fpkm[,1], Pg4.fpkm[,2], method= "spearman") #ties

plot(Pg4.fpkm[,1], Pg4.fpkm[,2])

colnames(Pg4.fpkm) <- c("AMBRA","MYC")
pdf("plots/Pg4_fpkm_correlations_heatmap.pdf")
heatmap(t(Pg4.fpkm[1:20,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pg4.fpkm[21:41,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pg4.fpkm[42:64,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

#weak correlation with kendall and pearson

cor.test (Pwnt.fpkm[,1], Pwnt.fpkm[,2], method= "pearson") 
cor.test (Pwnt.fpkm[,1], Pwnt.fpkm[,2], method= "kendall") 
cor.test (Pwnt.fpkm[,1], Pwnt.fpkm[,2], method= "spearman") 

plot(Pwnt.fpkm[,1], Pwnt.fpkm[,2])

colnames(Pwnt.fpkm) <- c("AMBRA","MYC")
pdf("plots/Pwnt_fpkm_correlations_heatmap.pdf")
heatmap(t(Pwnt.fpkm[1:20,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

#weak correlation with kendall and pearson

#plot data

pdf("plots/P.fpkm_four_AMBRA.pdf")
boxplot(Pshh.fpkm[,1],Pg4.fpkm[,1], Pg3.fpkm[,1],Pwnt.fpkm[,1])
dev.off()
pdf("plots/P.fpkm_four_MYC.pdf")
boxplot(Pshh.fpkm[,2],Pg4.fpkm[,2], Pg3.fpkm[,2],Pwnt.fpkm[,2])
dev.off()

#### verify if the data are normally distributed ####

#evaluate data distribution
d<- density(Pshh.fpkm[,1])
plot(d)
#skewness on the left
#Q-Q plot
ggqqplot(Pshh.fpkm[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pshh.fpkm[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.015 - data not normally distributed!
d<- density(Pshh.fpkm[,2])
plot(d)
#skewness to the right 
#Q-Q plot
ggqqplot(Pshh.fpkm[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pshh.fpkm[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value=0.14 


#evaluate data distribution
d<- density(Pg4.fpkm[,1])
plot(d)

#Q-Q plot
ggqqplot(Pg4.fpkm[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pg4.fpkm[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.27
d<- density(Pg4.fpkm[,2])
plot(d)
#skewness on the right 
#Q-Q plot
ggqqplot(Pg4.fpkm[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pg4.fpkm[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value=1.762e-05 not normally distributed!


#evaluate data distribution
d<- density(Pg3.fpkm[,1])
plot(d)

#Q-Q plot
ggqqplot(Pg3.fpkm[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pg3.fpkm[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.02 not normally distributed
d<- density(Pg3.fpkm[,2])
plot(d)

#Q-Q plot
ggqqplot(Pg3.fpkm[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pg3.fpkm[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value=0.55 

#evaluate data distribution
d<- density(Pwnt.fpkm[,1])
plot(d)
#skewness to the left
#Q-Q plot
ggqqplot(Pwnt.fpkm[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pwnt.fpkm[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.91
d<- density(Pwnt.fpkm[,2])
plot(d)
#skewness to the left
#Q-Q plot
ggqqplot(Pwnt.fpkm[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pwnt.fpkm[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value=0.04 not normally distributed!

#t-test pairwise
t.test(as.vector(Pshh.fpkm[,1]), as.vector(Pg3.fpkm[,1]))  # p-value = 0.00254
t.test(as.vector(Pg4.fpkm[,1]), as.vector(Pg3.fpkm[,1])) # p-value = 0.008657
t.test(as.vector(Pwnt.fpkm[,1]), as.vector(Pg3.fpkm[,1])) # p-value = 1.937e-06
t.test(as.vector(Pwnt.fpkm[,1]), as.vector(Pshh.fpkm[,1])) # p-value = 7.862e-08
t.test(as.vector(Pwnt.fpkm[,1]), as.vector(Pg4.fpkm[,1])) #  p-value = 1.332e-07
t.test(as.vector(Pg4.fpkm[,1]), as.vector(Pshh.fpkm[,1])) #  p-value = 0.543   NOT SIGNIFICANT

#wilcox test - better when there are deviation from normality
wilcox.test(as.vector(Pshh.fpkm[,1]), as.vector(Pg3.fpkm[,1]))  #  p-value = 0.0224
wilcox.test(as.vector(Pg4.fpkm[,1]), as.vector(Pg3.fpkm[,1])) # p-value = 0.03793
wilcox.test(as.vector(Pwnt.fpkm[,1]), as.vector(Pg3.fpkm[,1])) # ties
wilcox.test(as.vector(Pwnt.fpkm[,1]), as.vector(Pshh.fpkm[,1])) #p-value = 6.665e-11
wilcox.test(as.vector(Pwnt.fpkm[,1]), as.vector(Pg4.fpkm[,1]))  #p-value = 1.417e-08
wilcox.test(as.vector(Pg4.fpkm[,1]), as.vector(Pshh.fpkm[,1])) # p-value = 0.5325  NOT SIGNIFICANT


#### Pfister mas5 databasets

#### read the data and check group sizes ####

Pshh.mas5 <- read.table("downloaded_data/Pfister_mas5.0/pfister_shh_mas5_log2.txt", sep = '\t', header = FALSE)
View(Pshh.mas5)
###check the group size - Pshh.mas5 should contain 59 samples
dim(Pshh.mas5)
# 2 61
Pshh.mas5 <- as.data.frame(Pshh.mas5)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pshh.mas5 <- select(Pshh.mas5, V3:V61)
Pshh.mas5 <- t(Pshh.mas5)
dim(Pshh.mas5)

Pg3.mas5 <- read.table("downloaded_data/Pfister_mas5.0/pfister_g3_mas5_log2.txt", sep = '\t', header = FALSE)
View(Pg3.mas5)
###check the group size - Pshh.fpkm should contain 56 samples
dim(Pg3.mas5)
# 2 58
Pg3.mas5 <- as.data.frame(Pg3.mas5)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pg3.mas5 <- select(Pg3.mas5, V3:V58)
Pg3.mas5 <- t(Pg3.mas5)
dim(Pg3.mas5)

Pg4.mas5 <- read.table("downloaded_data/Pfister_mas5.0/pfister_g4_mas5_log2.txt", sep = '\t', header = FALSE)
View(Pg4.mas5)
###check the group size - Pg4.mas5 should contain 91 samples
dim(Pg4.mas5)
# 2 93
Pg4.mas5 <- as.data.frame(Pg4.mas5)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pg4.mas5 <- select(Pg4.mas5, V3:V93)
Pg4.mas5 <- t(Pg4.mas5)
dim(Pg4.mas5)

Pwnt.mas5 <- read.table("downloaded_data/Pfister_mas5.0/pfister_wnt_mas5_log2.txt", sep = '\t', header = FALSE)
View(Pwnt.mas5)
###check the group size - Pwnt.mas5 should contain 17 samples
dim(Pwnt.mas5)
# 2 19
Pwnt.mas5 <- as.data.frame(Pwnt.mas5)
#to remove the duplicate column with gene names
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Pwnt.mas5 <- select(Pwnt.mas5, V3:V19)
Pwnt.mas5 <- t(Pwnt.mas5)
dim(Pwnt.mas5)
#par(mar=c(1,1,1,1))


### correlation between AMBRA and MYC ###

cor.test (Pshh.mas5[,1], Pshh.mas5[,2], method= "pearson")
cor.test (Pshh.mas5[,1], Pshh.mas5[,2], method= "kendall")
cor.test (Pshh.mas5[,1], Pshh.mas5[,2], method= "spearman")

plot(Pshh.mas5[,1], Pshh.mas5[,2])

colnames(Pshh.mas5) <- c("AMBRA","MYC")
pdf("plots/Pshh_mas5_correlations_heatmap.pdf")
heatmap(t(Pshh.mas5[1:20,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pshh.mas5[21:41,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pshh.mas5[42:59,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()


cor.test (Pg3.mas5[,1], Pg3.mas5[,2], method= "pearson")
cor.test (Pg3.mas5[,1], Pg3.mas5[,2], method= "kendall")
cor.test (Pg3.mas5[,1], Pg3.mas5[,2], method= "spearman")

plot(Pg3.mas5[,1], Pg3.mas5[,2])

colnames(Pg3.mas5) <- c("AMBRA","MYC")
pdf("plots/Pg3_mas5_correlations_heatmap.pdf")
heatmap(t(Pg3.mas5[1:20,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pg3.mas5[21:56,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

cor.test (Pg4.mas5[,1], Pg4.mas5[,2], method= "pearson")
cor.test (Pg4.mas5[,1], Pg4.mas5[,2], method= "kendall")
cor.test (Pg4.mas5[,1], Pg4.mas5[,2], method= "spearman") 

plot(Pg4.mas5[,1], Pg4.mas5[,2])

colnames(Pg4.mas5) <- c("AMBRA","MYC")
pdf("plots/Pg4_mas5_correlations_heatmap.pdf")
heatmap(t(Pg4.mas5[1:20,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pg4.mas5[21:41,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pg4.mas5[42:62,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Pg4.mas5[63:91,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

cor.test (Pwnt.mas5[,1], Pwnt.mas5[,2], method= "pearson")  #correlation 0.46 p-value = 0.06319
cor.test (Pwnt.mas5[,1], Pwnt.mas5[,2], method= "kendall") #cor. 0.27
cor.test (Pwnt.mas5[,1], Pwnt.mas5[,2], method= "spearman") #cor. 0.39

plot(Pwnt.fpkm[,1], Pwnt.fpkm[,2])

colnames(Pwnt.mas5) <- c("AMBRA","MYC")
pdf("plots/Pwnt_mas5_correlations_heatmap.pdf")
heatmap(t(Pwnt.mas5[1:17,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

# weak correlation

#plot all four
pdf("plots/P_mas5_four_AMBRA.pdf")
boxplot(Pshh.mas5[,1],Pg4.mas5[,1], Pg3.mas5[,1],Pwnt.mas5[,1])
dev.off()
pdf("plots/P.mas5_four_MYC.pdf")
boxplot(Pshh.mas5[,2],Pg4.mas5[,2], Pg3.mas5[,2],Pwnt.mas5[,2])
dev.off()

#### verify if the data are normally distributed ####

#evaluate data distribution
d<- density(Pshh.mas5[,1])
plot(d)
#skewness to the right
#Q-Q plot
ggqqplot(Pshh.mas5[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pshh.mas5[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained  0.3903 
d<- density(Pshh.mas5[,2])
plot(d)
#skewness to the right 
#Q-Q plot
ggqqplot(Pshh.mas5[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pshh.mas5[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value= 0.5896 

#evaluate data distribution
d<- density(Pg4.mas5[,1])
plot(d)

#Q-Q plot
ggqqplot(Pg4.mas5[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pg4.mas5[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.02632 - not normally distributed
d<- density(Pg4.mas5[,2])
plot(d)
#skewness on the right 
#Q-Q plot
ggqqplot(Pg4.mas5[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pg4.mas5[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value=0.006855 not normally distributed!


#evaluate data distribution
d<- density(Pg3.mas5[,1])
plot(d)

#Q-Q plot
ggqqplot(Pg3.mas5[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pg3.mas5[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.9295
d<- density(Pg3.mas5[,2])
plot(d)

#Q-Q plot
ggqqplot(Pg3.mas5[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pg3.mas5[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value= 0.2623 

#evaluate data distribution
d<- density(Pwnt.mas5[,1])
plot(d)

#Q-Q plot
ggqqplot(Pwnt.mas5[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pwnt.mas5[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.804
d<- density(Pwnt.mas5[,2])
plot(d)

#Q-Q plot
ggqqplot(Pwnt.mas5[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Pwnt.mas5[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value= 0.03069 not normally distributed!

#t-test pairwise
t.test(as.vector(Pshh.mas5[,1]), as.vector(Pg3.mas5[,1]))  # p-value = 0.03510302
t.test(as.vector(Pg4.mas5[,1]), as.vector(Pg3.mas5[,1])) # p-value =  0.03301
t.test(as.vector(Pwnt.mas5[,1]), as.vector(Pg3.mas5[,1])) # p-value =  6.361e-06
t.test(as.vector(Pwnt.mas5[,1]), as.vector(Pshh.mas5[,1])) # p-value = 2.084e-06
t.test(as.vector(Pwnt.mas5[,1]), as.vector(Pg4.mas5[,1])) #  p-value = 1.025e-06
t.test(as.vector(Pg4.mas5[,1]), as.vector(Pshh.mas5[,1])) #  p-value = 0.3759   NOT SIGNIFICANT

#wilcox test - better when there are deviation from normality
wilcox.test(as.vector(Pshh.mas5[,1]), as.vector(Pg3.mas5[,1]))  #  p-value = 0.1191 NOT SIGNIFICANT
wilcox.test(as.vector(Pg4.mas5[,1]), as.vector(Pg3.mas5[,1])) # p-value = 0.05708
wilcox.test(as.vector(Pwnt.mas5[,1]), as.vector(Pg3.mas5[,1])) # p-value = 1.495e-07
wilcox.test(as.vector(Pwnt.mas5[,1]), as.vector(Pshh.mas5[,1])) #p-value =  1.521e-08
wilcox.test(as.vector(Pwnt.mas5[,1]), as.vector(Pg4.mas5[,1]))  #p-value = 5.293e-09
wilcox.test(as.vector(Pg4.mas5[,1]), as.vector(Pshh.mas5[,1])) # p-value = 0.7539  NOT SIGNIFICANT


#### Gilbertson mas5 databasets

#### read the data and check group sizes ####

Gshh.mas5 <- read.table("downloaded_data/Gilbertson_mas5.0/gilbertson_shh_mas5_log2.txt", sep = '\t', header = FALSE)
View(Gshh.mas5)
###check the group size - Gshh.mas5 should contain 10 samples
dim(Gshh.mas5)
# 2 12
Gshh.mas5 <- as.data.frame(Gshh.mas5)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Gshh.mas5 <- select(Gshh.mas5, V3:V12)
Gshh.mas5 <- t(Gshh.mas5)
dim(Gshh.mas5)

Gg3.mas5 <- read.table("downloaded_data/Gilbertson_mas5.0/gilbertson_g3_mas5_log2.txt", sep = '\t', header = FALSE)
View(Gg3.mas5)
###check the group size - Pshh.fpkm should contain 16 samples
dim(Gg3.mas5)
# 2 18
Gg3.mas5 <- as.data.frame(Gg3.mas5)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Gg3.mas5 <- select(Gg3.mas5, V3:V18)
Gg3.mas5 <- t(Gg3.mas5)
dim(Gg3.mas5)


Gg4.mas5 <- read.table("downloaded_data/Gilbertson_mas5.0/gilbertson_g4_mas5_log2.txt", sep = '\t', header = FALSE)
View(Gg4.mas5)
###check the group size - Gg4.mas5 should contain 39 samples
dim(Gg4.mas5)
# 2 41
Gg4.mas5 <- as.data.frame(Gg4.mas5)
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Gg4.mas5 <- select(Gg4.mas5, V3:V41)
Gg4.mas5 <- t(Gg4.mas5)
dim(Gg4.mas5)

Gwnt.mas5 <- read.table("downloaded_data/Gilbertson_mas5.0/gilbertson_wnt_mas5_log2.txt", sep = '\t', header = FALSE)
View(Gwnt.mas5)
###check the group size - Gwnt.mas5 should contain 8 samples
dim(Gwnt.mas5)
# 2 10
Gwnt.mas5 <- as.data.frame(Gwnt.mas5)
#to remove the duplicate column with gene names
#to remove the duplicate column with gene names and transpose the matrix to have genes on the columns and samples on the rows
#first column is AMBRA, second column is MYC
Gwnt.mas5 <- select(Gwnt.mas5, V3:V10)
Gwnt.mas5 <- t(Gwnt.mas5)
dim(Gwnt.mas5)
#par(mar=c(1,1,1,1))


### correlation between AMBRA and MYC ###

cor.test (Gshh.mas5[,1], Gshh.mas5[,2], method= "pearson")
cor.test (Gshh.mas5[,1], Gshh.mas5[,2], method= "kendall") # cor. 0.33
cor.test (Gshh.mas5[,1], Gshh.mas5[,2], method= "spearman") #cor. 0.38

plot(Gshh.mas5[,1], Gshh.mas5[,2])

colnames(Gshh.mas5) <- c("AMBRA","MYC")
pdf("plots/Gshh_mas5_correlations_heatmap.pdf")
heatmap(t(Gshh.mas5[1:10,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

# weak correlation

cor.test (Gg3.mas5[,1], Gg3.mas5[,2], method= "pearson") #cor. 0.39
cor.test (Gg3.mas5[,1], Gg3.mas5[,2], method= "kendall") #cor. 0.25
cor.test (Gg3.mas5[,1], Gg3.mas5[,2], method= "spearman") #cor. 0.40

plot(Gg3.mas5[,1], Gg3.mas5[,2])

colnames(Gg3.mas5) <- c("AMBRA","MYC")
pdf("plots/Gg3_mas5_correlations_heatmap.pdf")
heatmap(t(Gg3.mas5[1:16,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

# weak correlation

cor.test (Gg4.mas5[,1], Gg4.mas5[,2], method= "pearson") 
cor.test (Gg4.mas5[,1], Gg4.mas5[,2], method= "kendall") 
cor.test (Gg4.mas5[,1], Gg4.mas5[,2], method= "spearman") 

plot(Gg4.mas5[,1], Gg4.mas5[,2])

colnames(Gg4.mas5) <- c("AMBRA","MYC")
pdf("plots/Gg4_mas5_correlations_heatmap.pdf")
heatmap(t(Gg4.mas5[1:20,]), cexRow = 0.9, Colv=NA, Rowv=NA)
heatmap(t(Gg4.mas5[21:39,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

cor.test (Gwnt.mas5[,1], Gwnt.mas5[,2], method= "pearson")  
cor.test (Gwnt.mas5[,1], Gwnt.mas5[,2], method= "kendall")
cor.test (Gwnt.mas5[,1], Gwnt.mas5[,2], method= "spearman") 

plot(Gwnt.mas5[,1], Gwnt.mas5[,2])

colnames(Gwnt.mas5) <- c("AMBRA","MYC")
pdf("plots/Gwnt_mas5_correlations_heatmap.pdf")
heatmap(t(Gwnt.mas5[1:8,]), cexRow = 0.9, Colv=NA, Rowv=NA)
dev.off()

#plot all four
pdf("plots/G_mas5_four_AMBRA.pdf")
boxplot(Gshh.mas5[,1],Gg4.mas5[,1], Gg3.mas5[,1],Gwnt.mas5[,1])
dev.off()
pdf("plots/G.mas5_four_MYC.pdf")
boxplot(Gshh.mas5[,2],Gg4.mas5[,2], Gg3.mas5[,2],Gwnt.mas5[,2])
dev.off()

#### verify if the data are normally distributed ####

#evaluate data distribution
d<- density(Gshh.mas5[,1])
plot(d)

#Q-Q plot
ggqqplot(Gshh.mas5[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Gshh.mas5[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained  0.01976 not normally distributed 
d<- density(Gshh.mas5[,2])
plot(d)

#Q-Q plot
ggqqplot(Gshh.mas5[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Gshh.mas5[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value= 0.8622 

#evaluate data distribution
d<- density(Gg4.mas5[,1])
plot(d)

#Q-Q plot
ggqqplot(Gg4.mas5[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Gg4.mas5[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained  0.002488 - not normally distributed
d<- density(Gg4.mas5[,2])
plot(d)

#Q-Q plot
ggqqplot(Gg4.mas5[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Gg4.mas5[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value= 0.02317 not normally distributed!


#evaluate data distribution
d<- density(Gg3.mas5[,1])
plot(d)

#Q-Q plot
ggqqplot(Gg3.mas5[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Gg3.mas5[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.764
d<- density(Gg3.mas5[,2])
plot(d)

#Q-Q plot
ggqqplot(Gg3.mas5[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Gg3.mas5[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value= 0.648 

#evaluate data distribution
d<- density(Gwnt.mas5[,1])
plot(d)

#Q-Q plot
ggqqplot(Gwnt.mas5[,1])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Gwnt.mas5[,1])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained 0.542
d<- density(Gwnt.mas5[,2])
plot(d)

#Q-Q plot
ggqqplot(Gwnt.mas5[,2])
#A 45-degree reference line is also plotted.As all the points fall approximately along this reference line, we can assume normality.
#normality test with Shapiro-Wilk???s method, based on the correlation between the data and the corresponding normal scores.
shapiro.test(Gwnt.mas5[,2])
#if p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution,we can assume the normality.
#obtained p-value= 0.2506

#t-test pairwise
t.test(as.vector(Gshh.mas5[,1]), as.vector(Gg3.mas5[,1]))  # p-value = 0.905 NOT SIGNIFICANT
t.test(as.vector(Gg4.mas5[,1]), as.vector(Gg3.mas5[,1])) # p-value =  0.9656 NOT SIGNIFICANT
t.test(as.vector(Gwnt.mas5[,1]), as.vector(Gg3.mas5[,1])) # p-value =  0.01454
t.test(as.vector(Gwnt.mas5[,1]), as.vector(Gshh.mas5[,1])) # p-value = 0.0154
t.test(as.vector(Gwnt.mas5[,1]), as.vector(Gg4.mas5[,1])) #  p-value = 0.01349
t.test(as.vector(Gg4.mas5[,1]), as.vector(Gshh.mas5[,1])) #  p-value = 0.9178   NOT SIGNIFICANT

#wilcox test - better when there are deviation from normality
wilcox.test(as.vector(Gshh.mas5[,1]), as.vector(Gg3.mas5[,1]))  #  p-value = 0.1191 NOT SIGNIFICANT
wilcox.test(as.vector(Gg4.mas5[,1]), as.vector(Gg3.mas5[,1])) # p-value = 0.05708
wilcox.test(as.vector(Gwnt.mas5[,1]), as.vector(Gg3.mas5[,1])) # p-value = 1.495e-07
wilcox.test(as.vector(Gwnt.mas5[,1]), as.vector(Gshh.mas5[,1])) #p-value =  1.521e-08
wilcox.test(as.vector(Gwnt.mas5[,1]), as.vector(Gg4.mas5[,1]))  #p-value = 5.293e-09
wilcox.test(as.vector(Gg4.mas5[,1]), as.vector(Gshh.mas5[,1])) # p-value = 0.7539  NOT SIGNIFICANT
