
################################################################################

# packages 

################################################################################

# data organisation

library(dplyr)
library(tidyverse)
library(ggthemes)

# RNA-seq specific packages  

library(edgeR) # import, organise, filter and normalise the data, 
library(limma) # voom method, linear modelling and empirical Bayes moderation for DEG and GSE
library(Glimma) #nteractive exploration of the results to query individual samples 
library(org.Hs.eg.db) # human annotations 

# visualisation

library(RColorBrewer) # colour palette 
library(gplots) # plots 
library(ggfortify) # pca
library(viridis) # for one palette with 242 col


################################################################################

# folder 

################################################################################

setwd("C:/Users/sabrinai/OneDrive - The University of Melbourne/PHD/Chapter2/3.AllenBrain")


################################################################################

# about the data and the experiment 

################################################################################

# soruce: https://human.brain-map.org/static/download (raw counts )
# tissue: cortex and subcortex of human
# cell types :  
## 29 substrucutres
## 10 main strucutures 
# replicates: 8  
# libraries : 29x8 + (10 extra) = 242
# sequencing:  Illumina HiSeq 2000 to obtain 50 bp single-end
# depth: 30 MM reads/sample
# aligned to the human genome using RNA-Seq by Expectation-Maximization 
# transcripts defined using the knownGene table from the UCSC Genome Browser
# workflow based on https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

################################################################################

# Objective of this code 

################################################################################

# I will be performing quality-control of the RNA seq counts

# I will go through quality-control checkpoints, such as 

  # transforming the data to a scale that allows good exploration 
  
  # filtering genes with low expression
  
  # normlising the libraries to be more comparable 
  
  # unsupervised clustering to see how the samples group together 

################################################################################

# section 1: import the data 

################################################################################

# step 1:  loading RAW gene counts 

# sample 1

ab_count_1 <- read.csv("Data/RNAseqCounts.csv",
                       header = F,
                       check.names = F,
                       row.names = 1)

# samples annotation 

ab_samples_1 <- read.csv("Data/SampleAnnot.csv",
                         header = TRUE,
                         check.names = F,
                         row.names = 1)

ab_tissues_1 <- data.frame(ab_samples_1[ , c(6,7)]) # extract tissue info 

ab_tissues_1$SAMPID <- rownames(ab_tissues_1) # add sample id 

ab_tissues_1$sub_structure <- as.factor(ab_tissues_1$sub_structure) # make tissue infor factor 

colnames(ab_count_1) <-  rownames(ab_tissues_1) # fit colnames and rownames 

# sample 2

ab_count_2 <- read.csv("Data2/RNAseqCounts.csv",
                       header = F,
                       check.names = F,
                       row.names = 1)


# samples annotation 

ab_samples_2 <- read.csv("Data2/SampleAnnot.csv",
                         header = TRUE,
                         check.names = F,
                         row.names = 1)

ab_tissues_2 <- data.frame(ab_samples_2[ , c(6,7)]) # extract tissue info 

ab_tissues_2$SAMPID <- rownames(ab_tissues_2) # ass sample id 

ab_tissues_2$sub_structure <- as.factor(ab_tissues_2$sub_structure) # make tissue info factor 

colnames(ab_count_2) <-  rownames(ab_tissues_2)  # fit colnames and rownames 


# combine

ncol(ab_count_1) # 121 

ncol(ab_count_2) # 121 

ab_count_whole <- cbind(ab_count_1, ab_count_2)

ncol(ab_count_whole) # 242 

################################################################################

# step 2: organising the sample information 

# to examine biological (cell type) and technical (lane) variables 

# cell type information 

ab_tissues <- rbind(ab_tissues_1, ab_tissues_2)

nrow(ab_tissues) == ncol(ab_count_whole) # same number of rows and columns 

tissues <- ab_tissues[, c(1:3)] # making another df for more manipulation 

head(tissues) # has main and sub structures and sample id 

# sequencing information

ab_samples <- rbind(ab_samples_1, ab_samples_2)

nrow(ab_samples) == ncol(ab_count_whole) # same 

samples <- ab_samples [, c(2,15)] # sampid and sequencing id extract

head(samples) # worked 

###############################################################################

# step 3: organising gene-level info

ab_annot_1 <- data.frame(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = rownames(ab_count_1),
                                               columns = c ("GENENAME"),
                                               keytype="GENENAME"))

length(unique(ab_annot_1$GENENAME)) # no duplicates 

nrow(ab_annot_1) == nrow(ab_count_1) # T


ab_annot_2 <- data.frame(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = rownames(ab_count_2),
                                               columns = c ("GENENAME"),
                                               keytype="GENENAME"))


length(unique(ab_annot_2$GENENAME)) # no duplicates again

nrow(ab_annot_2) == nrow(ab_count_2) # T

summary(ab_annot_1$GENENAME == ab_annot_2$GENENAME) #   same genes in both data  

ab_annot_whole <- ab_annot_1 # for DGE

# reading as a DGE list object 

x <- DGEList(counts = ab_count_whole, genes = ab_annot_whole)

class(x) # is DGE list 

x$samples # sample information 

# library sizes are automatically calculated for each sample and normalisation factors are set to 1

x$counts # count data 

dim(x)

# 22318 rows: Entrez gene identifiers (IDs); 242 columns

# barplot of library sizes

# (not a step from the paper)

par(mfrow=c(1,1))

plot_col <- brewer.pal(n =12, name = "Paired")

barplot(colSums(x$counts[, 21:40]) / 1e+6, 
        las = 2, 
        col = plot_col,
        main = "library size in millions") 

# library size is variable 

##############################################################################

# section 2: preprocessing 

#############################################################################

# step 1: transformation

# library SIZE difference: libraries sequenced at a greater depth will result in higher counts.

# that's why we transform raw counts onto a scale that accounts for such library size differences

# popular options; cpm, logcpm, rpkm and fpkm

# cpm and logcpm can be calculated with just count matrix. gene length is assumed to be constant to compare between conditions.

# rpkm and fpkm accounts for gene length. 

# we won't be comparing expression across multiple genes or drawing conclusions on absolute levels of expression

cpm <- cpm(x) # cpm transformation 

lcpm <- log2(cpm(x$counts) + 1) # log2 + cpm transformation 

# used for exploratory plots 

# when log == T, cpm adds an offset

# the offset is  2/L where 2 is the "prior count" and L is the average library size in millions,

# so it is actually log2(CPM + 2/L). (I used log2(CPM+ 1), voom used log2(cpm+ 0.5))

# I used 1 as offest to firstly, remove negative valus for ease on interpretation 

# secondly, I got the idea from Jarny 

# they all serve the same purpose of avoiding taking the logarithm of 0

# thus avoiding - tive log2 values for genes with very low count 

L <- mean(x$samples$lib.size) * 1e-6 # mean library size is 23.35

M <- median(x$samples$lib.size) * 1e-6 # median library size is 22.65

c(L, M)

###############################################################################

# step 2: filtering

table(rowSums(x$counts==0)==242) 

# 1396 have 0 counts (6.5%)

# 1396/22318*100 is 6.5 %

# 6.5 % of genes in this dataset have zero counts across all nine samples.

# why do we remove genes with very low reads/counts

# genes that not expressed at a biologically meaningful level in any condition are not of interest

# to reliably estimate the mean-variance trend

# to reduce the number of tests looking for DEG

# The filterByExpr function in the edgeR package provides an automatic way to filter genes, 

# while keeping as many genes as possible with worthwhile counts.

# keeps genes with about 10 read counts or more in a minimum number of samples, 

# where the number of samples is chosen according to the minimum group sample size. 

# The actual filtering uses CPM values rather than counts in order 

# to avoid giving preference to samples with large library sizes.

# cutoff = 10/ median library size

# if median library size is bigger cutoff becomes smaller

# cz because larger library sizes == better resolution to explore more genes at lower expression levels

# but smaller library size == less resolution into marginal genes 

# For this dataset, 

# the median library size is about 22.65 million 

# so the filterByExpr function keeps genes that have a CPM of 10/51 = 0.44 cpm

# in at least 8 samples, because there are three replicates for each group

# my detection threshold was 1.3 

keep.exprs <- filterByExpr(x, group = tissues$main_structure)  

fil_x <- x[keep.exprs,, keep.lib.sizes=FALSE]

dim(fil_x)

(22318-17340)/22318 # app 22% got removed 

# filtering counts below 1.3 in 8 samples 

thresh <- lcpm > 1.3

head(thresh) # logical 

keep <- (rowSums(thresh) > 8) # at least 1/10 samples TRUES i

table(keep)

summary(keep) 

lcpm_keep <- lcpm [keep, ]

(22318-15820)/22318 # app 29% got removed 

# effect of filtering visualise 

lcpm.cutoff <- log2(10/M + 2/L) # cutoff as used by edgeR  

col <- brewer.pal(11, "Paired")

# this lcpm is x from before filtering 

# lcpm post filtering 

lcpm_fil <- log2(cpm(fil_x$counts) + 1)

nsamples_i <- 11 # number of samples 

# library 232 to 242 from here 

samplenames_i <- as.vector(as.character(ab_tissues$SAMPID[232:242]))

data_u <- as.data.frame(lcpm[, 232:242]) 

data_f <- as.data.frame(lcpm_fil[, 232:242])

data_k <- as.data.frame(lcpm_keep[, 232:242])

par(mfrow = c(1,3))

# no filtering  

plot(density(data_u[,1]),
     col= col[1],
     lwd=2,
     xlab="log2(cpm+1)",
     ylab = "density",
     ylim=c(0,0.4),
     las=2,
     main= "Samples 232:242")

for (j in 2:nsamples_i)
{
  den_j <- density(data_u[,j])
  lines(den_j$x,
        den_j$y,
        col= col[j],
        lwd=2)
}
legend("topright", samplenames_i, text.col= col, bty="n")

# filtering with edgeR 

plot(density(data_f[,1]),
     col= col[1],
     lwd=2,
     xlab="log2(cpm+1)",
     ylab = "density",
     ylim=c(0,0.4),
     las=2,
     main= "filterbyExp")

for (j in 2:nsamples_i)
{
  den_j <- density(data_f[,j])
  lines(den_j$x,
        den_j$y,
        col= col[j],
        lwd=2)
}
legend("topright", samplenames_i, text.col= col, bty="n")

# filtering with the elbow  

plot(density(data_k[,1]),
     col= col[1],
     lwd=2,
     xlab="log2(cpm+1)",
     ylab = "density",
     ylim=c(0,0.4),
     las=2,
     main= "1.3")

for (j in 2:nsamples_i)
{
  den_j <- density(data_k[,j])
  lines(den_j$x,
        den_j$y,
        col= col[j],
        lwd=2)
}
legend("topright", samplenames_i, text.col= col, bty="n")

################################################################################

# step 3: normalising 

# normalisation addresses the expression distribution across libraries

# which should be similar 

# however, experimentally one batch can have more/less expressed genes than the other 

# distribution of log-cpm can be visulaised as density (above) or boxplot

# calcNormFactors function in edgeR used here to do TMM (trimmed mean of M-values) normalisation

# this method uses the norm.factors stored in the DGE$samples to scale the libraries

# in this experiment, all norm.factors are close to 1 so lil effect 

keep_x <- x[keep, ] # filtering the dge with my elbow method threshold 

dim(keep_x) # 15820

dim(fil_x) # 17340

# again, 2000 more genes removed 

norm_x <- calcNormFactors(keep_x, method = c("TMM")) # TMM normlisation 

norm_x$samples$norm.factors # norm.factors are mostly close to 1

# lets plot a histogram to see how many normalizing factors deviate too far from 1

par(mfrow = c(1,1))

hist(keep_x$samples$norm.factors, main="norm.factors of libraries post calcNormFactors()",xlab="norm.factors")

# in the unnormalised dge, all normalising factors = 1

hist(norm_x$samples$norm.factors, main="norm.factors of libraries post calcNormFactors()",xlab="norm.factors")

# post-normalising, a few normalising factors deviate from 1

# which ones are those?

length(rownames(norm_x$samples[which(norm_x$samples$norm.factors < 0.95), ])) # 12 samples 

length(rownames(norm_x$samples[which(norm_x$samples$norm.factors > 1.05), ])) # 10 samples 

100 - ((22/242)*100) # 91% normlising factors are within 0.95 to 1.05, but 9% libraries differ 

# logcpm tranformation of the normalised library 

lcpm_norm <- log2(cpm(norm_x$counts) + 1)

# the effect of normalising on those 22 outlier samples 

outlier_samp <- c("S010157_L7.LB15",
                   "S010217_L5.LB6", 
                   "S010363_L27.LB4",
                   "S010416_L6.LB8",
                   "S010417_L9.LB5",
                   "S010419_L3.LB16", 
                   "S010420_L8.LB14" , 
                   "S010536_L1.LB20",  
                   "S010536_L8.LB20",
                   "S010543_L2.LB4",  
                   "S010556_L2.LB1",   
                   "S010557_L35.LB23",
                   "S010066_L6.LB12",
                   "S020032_L7.LB13", 
                   "S020055_L3.LB12",
                   "S020064_L4.LB25", 
                   "S020205_L2.LB10", 
                   "S020206_L6.LB7",
                   "S020208_L2.LB7",
                   "S020210_L7.LB15", 
                   "S020215_L4.LB23", 
                   "S020697_L1.LB3")

length(outlier_samp)

par(mfrow=c(1,2))

# unnormalised 

boxplot(lcpm_keep[, colnames(lcpm_keep) %in% outlier_samp], las=2, col=col, main="")
title(main="A. Unnormalised Samples",ylab="log(cpm+1)")

# normalised 
boxplot(lcpm_norm[, colnames(lcpm_norm) %in% outlier_samp], las=2, col=col, main="")
title(main="B. Normalised Samples",ylab="log(cpm+1)")

# visually there is little change 

# but if we examine the count values 

df_un <- lcpm_keep[, colnames(lcpm_keep) %in% outlier_samp]

df_norm <- lcpm_norm[, colnames(lcpm_norm) %in% outlier_samp]

head(df_un[ ,1:4])

head(df_norm[ ,1:4])

# we see the count values have changed

# and if we look at total count for each sample 

table(colSums(df_un) == colSums(df_norm)) # values are different 

################################################################################

#  section 3: unsupervised clustering 

###############################################################################

# to look at how the samples group together without any supervision

# in an ideal world

# samples would cluster well within the primary condition of interest, 

# replicates should cluster together 

# in real world 

# any sample straying far from its group could be identified and followed up for sources of error or extra variation.


# the first dimension represents the leading-fold-change that best separates samples 

# and explains the largest proportion of variation in the data

# and so on and so forth for each dimension


# step 1: MDS 

par(mfrow=c(1,2))

# setting colour 

tissues$main_structure <- as.factor(tissues$main_structure)

levels(tissues$main_structure) # 10 

tissues$donor <- c("H0351.2001")

tissues$donor[122:242] <- c("H0351.2002")

tissues$donor <- as.factor(tissues$donor)

levels(tissues$donor)

samples$rnaseq_run_id <- as.factor(samples$rnaseq_run_id)

levels(samples$rnaseq_run_id) # 10 

col.group <- c("#9e0142", # cbcx
               "#9970ab", # cgg
               "#f46d43", # Fl
               "#fdae61", # GP
               "#c51b7d", # Ins
               "#e6f598", # OL
               "#d1e5f0" , # PHG
               "#66c2a5", # PL
               "#3288bd", # Str
               "#5e4fa2") [tissues$main_structure]

col.lane <- c("#9e0142", # cbcx
              "#9970ab", # cgg
              "#f46d43", # Fl
              "#fdae61", # GP
              "#c51b7d", # Ins
              "#e6f598", # OL
              "#d1e5f0" , # PHG
              "#66c2a5", # PL
              "#3288bd", # Str
              "#5e4fa2")  [samples$rnaseq_run_id]

# MDS by cell type 

par(mfrow = c(1,2))

plotMDS(lcpm_norm, labels = tissues$main_structure, col=col.group)
legend("bottomleft", 
       legend =levels(tissues$main_structure), 
       cex = 0.8)
title(main="MDS of `H0351.2001`, `H0351.2002` by brain main structures")

# MDS by sequencing 

plotMDS(lcpm_norm, labels = samples$rnaseq_run_id, col=col.lane)
legend("bottomleft", 
       legend =levels(samples$rnaseq_run_id), 
       cex = 0.8)
title(main="B. Sample lanes")


# we see clear grouping of cortical, cerebellar, GP and Str samples 

# in my judgment, seq run id doesn't impact clsutering much 

# but donor may 

###############################################################################

# step 2: PCA

## Perform pca

pca <- prcomp(t(lcpm_norm))

## color scheme 

plot_main <-  c("#9e0142", # cbcx
                "#9970ab", # cgg
                "#f46d43", # Fl
                "#fdae61", # GP
                "#c51b7d", # Ins
                "#e6f598", # OL
                "#d1e5f0" , # PHG
                "#66c2a5", # PL
                "#3288bd", # Str
                "#5e4fa2")   # TL # brain tissue 

plot_donor <- c("#0072B2", "#D55E00") # donor id 

## Plot PCA

par(mfrow = c(1,1))


## by brain structure 

autoplot(pca, 
         x = 1,
         y = 2,
         # frame = TRUE,
         # frame.type = 'norm',
         data = tissues, 
         colour = "main_structure", 
         size = 4, 
         main = "PCA of `H0351.2001`, `H0351.2002` by brain main structures") +
  scale_fill_manual(values = plot_main) + 
  scale_color_manual(values =  plot_main) +
  theme_classic() +
  theme(text = element_text(family = "Arial"))

## by donor 

autoplot(pca, 
         x = 1,
         y = 2,
         data = tissues, 
         colour = "donor", 
         size = 4, 
         main = "PCA coloured by donor") +
  scale_fill_manual(values = plot_donor) + 
  scale_color_manual(values =  plot_donor) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# by run id 

autoplot(pca, 
         x = 1,
         y = 2,
         data = samples, 
         colour = "rnaseq_run_id", 
         size = 4, 
         main = "PCA coloured by sequencing batch") +
  scale_fill_manual(values = plot_col) + 
  scale_color_manual(values =  plot_col) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# cells do cluster by cortex, cerebellum, gp and striatum

# rna seq id has some effect along the 3rd dimension

# donor variation has an effect 

###############################################################################@

# step 5: modelling for all libraries

################################################################################

par(mfrow = c(1,1))

boxplot(lcpm_keep, las=2, col=col, main="")
title(main="Not Normalised",ylab="log2(cpm+1)")


boxplot(lcpm_norm, las=2, col=col, 
        main="")
title(main="Logcpm transformed normlised libraries",ylab="log2(cpm+1)")

# 
gg_data <- as.data.frame(lcpm_norm) %>% gather(key = "Sample", value = "Count")

head(gg_data)

gg_data$Sample <- as.factor(gg_data$Sample)

pal <- viridis_pal(option = "D")(242)

gg_data %>%
  ggplot() + 
  geom_boxplot(aes(x = Sample, y = Count), 
               fill = pal) +
  labs(title = "Samples from the Allen Brain Atlas data post quality control and normalisation",
       y = "Count (log2(cpm+1))")+
  scale_y_continuous(breaks =seq(0, 15, 1), limit = c(0, 15)) +
  theme_tufte() +
  theme(text = element_text(family = "Arial"),
    axis.text.x = element_blank())


# exporting

write.csv(norm_x$counts, "Data/NormalisedLibraries.csv")

write.csv(norm_x$counts, "Data/NormalisedGenes.csv")
