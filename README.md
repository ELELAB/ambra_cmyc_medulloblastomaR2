Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark

Repository associated to the publication:

Targeting cancer stem cells in medulloblastoma by inhibiting AMBRA1 dual function in autophagy and STAT3 signalling.
Nazio F, Po A, Abballe L, Ballabio C, Diomedi Camassei F, Bordi M, Camera A, Caruso S, Caruana I, Pezzullo M, Ferraina C, Milletti G, Gianesello M, Reddel S, De Luca CD, Ceglie D, Marinelli S, Campello S, Papaleo E, Miele E, Cacchione A, Carai A, Vinci M, Velardi E, De Angelis B, Tiberi L, Quintarelli C, Mastronuzzi A, Ferretti E, Locatelli F, Cecconi F.
Acta Neuropathol. 2021 Sep;142(3):537-564. doi: 10.1007/s00401-021-02347-7.


author of the repository: Elena Papaleo, elenap@cancer.dk

This repository contains normalized and log2 transformed data from the R2.mc.nl database with AMBRA1 and MYC expression levels in medulloblastoma.
The repository  was made with the intent of openly sharing both the raw input data used at the time of the analyses and the R-scripts employed to carry out the study.

Start using the rscripts/R2_medullo_analyses.R

Requirements:

R version 3.3.1 or higher Rstudio version 1.1.383 or higher
Bioconductor version 3.6 or higher

Other packages required:

CRAN:

ggplot
dplyr
ggpubr

BIOCONDUCTOR:
cemitool

NOTES:

a) We suggest to use Rstudio to run the scripts of interest so that you can follow the analyses one line at the time and digest the results.


INFORMATION ON DATA DOWNLOAD:
We used the R2 genomics platform (https://hgserver1.amc.nl/cgi-bin/r2/main.cgi) to identify suitable transcriptomics datasets with molecular and histological subtypes available for medulloblastoma. We selected the following datasets: Gilbertson (76 samples, MAS5.0 normalized, u133p2 array), Pfister (223 samples, MAS5.0 normalized, u133ps array), Pfister (167 samples, fpkm normalized, mb500rs1 array) and Cavalli (763 samples, rma_sketch normalized, hugene11t chip).
For each of them we download through the “data grabber” function the data for each subgroup separately and for the 52731_at and 202431_s_at probes for AMBRA1 and c-MYC, respectively, which are the probes with the largest coverage for the two genes. We selected log2 transformed data.

