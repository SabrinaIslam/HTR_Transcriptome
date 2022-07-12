################################################################################

# packages 

################################################################################

# RNA-seq specific packages  

library(edgeR) # import, organise, filter and normalise the data, 
library(limma) # voom method, linear modelling and empirical Bayes moderation for DEG and GSE
library(Glimma) #nteractive exploration of the results to query individual samples 
library(org.Hs.eg.db) # human annotations 

# visualisation

library(RColorBrewer) # colour palette 
library(gplots) # plots 
library(ggfortify) # pca

# heatmap 

library(ComplexHeatmap)
library(GetoptLong)
library(circlize)
library(dendextend)
library(corrplot)
library(factoextra)
library(paletteer)


# font issue 

windowsFonts("Arial" = windowsFont("Arial"))


################################################################################

# folder 

################################################################################

setwd("C:/Users/sabrinai/OneDrive - The University of Melbourne/PHD/Chapter2/3.AllenBrain")


################################################################################

# about the data and the experiment 

#########################################################################

# soruce: https://human.brain-map.org/static/download
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


################################################################################

# Objective of this code 

################################################################################

# I will filter and normalise the count data using previously established workflow

# then I will design a set of comparisons (the tissues)

# I will fit the data to the linear model. I will look for variability acorss tissues

# and account for donor variation

################################################################################

# section 1:import the data 

################################################################################


# step 1: reading count the data 

# loading RAW gene counts 
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

ab_tissues_1 <- data.frame(ab_samples_1[ , c(6,7)])

ab_tissues_1$SAMPID <- rownames(ab_tissues_1)

ab_tissues_1$sub_structure <- as.factor(ab_tissues_1$sub_structure)

colnames(ab_count_1) <-  rownames(ab_tissues_1) 

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

ab_tissues_2 <- data.frame(ab_samples_2[ , c(6,7)])

ab_tissues_2$SAMPID <- rownames(ab_tissues_2)

ab_tissues_2$sub_structure <- as.factor(ab_tissues_2$sub_structure)

colnames(ab_count_2) <-  rownames(ab_tissues_2) 


# combine

ncol(ab_count_1)

ncol(ab_count_2)

ab_count_whole <- cbind(ab_count_1, ab_count_2)

ncol(ab_count_whole)

##################################################

# step 2: organising the sample information 

# to examine biological (cell type) and technical (lane) variables 

# cell type information 

ab_tissues <- rbind(ab_tissues_1, ab_tissues_2)

nrow(ab_tissues) == ncol(ab_count_whole)

tissues <- ab_tissues[, c(1:3)]

head(tissues)

# sequencing information

ab_samples <- rbind(ab_samples_1, ab_samples_2)

nrow(ab_samples) == ncol(ab_count_whole)

samples <- ab_samples [, c(15)]

head(samples)

# join 

sampleinfo <- cbind(tissues, samples)

head(sampleinfo)

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

summary(ab_annot_1$GENENAME == ab_annot_2$GENENAME) # same 

ab_annot_whole <- ab_annot_1 # for DGE

# reading as a DGE list object 

x <- DGEList(counts = ab_count_whole, genes = ab_annot_whole)

class(x) # is DGE list 

dim(x)

################################################################################

# section 2: QC 

################################################################################

# log cpm transnformation 

lcpm <- log2(cpm(x$counts) + 1)


# filtering 

# filtering counts below 1.3 in 8 samples 

thresh <- lcpm > 1.3

head(thresh) # logical 

keep <- (rowSums(thresh) > 8) # at least 1 tissues (8 replicates) expressed 

table(keep)

summary(keep) 

lcpm_keep <- lcpm [keep, ]

(22318-15820)/22318 # app 29% got removed 

fil_k <- x[keep, ]

dim(fil_k) # 15820

###############################################################################

# normalisation 

x_norm <- calcNormFactors(fil_k, method = c("TMM")) # TMM normlisation 

x_norm$samples$norm.factors # norm.factors are close to 1

#################################################################################

### section 3: voom

################################################################################
# setting levels

sampleinfo$main_structure <- as.factor(sampleinfo$main_structure)

levels(sampleinfo$main_structure) # 10 

sampleinfo$donor <- c("H0351.2001")

sampleinfo$donor[122:242] <- c("H0351.2002")

sampleinfo$donor <- as.factor(sampleinfo$donor)

length(levels(sampleinfo$donor))

sampleinfo$samples <- as.factor(sampleinfo$samples)

levels(sampleinfo$samples)

# step 1: matrix 

# create design matrix 

# design matrix will accommodate the information from unsupervised clustering

# samples cluster together based on brain main structure, donor and lanes. 

# I am comparing the  brain structures are the variables 

# and I am blocking the covariates donor and lanes 

# thus, the brain structures are compared only within each blocks

# source: https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

# section 9.4.2 page 45 

design <- model.matrix(~0 + sampleinfo$main_structure + sampleinfo$donor + sampleinfo$samples) 

# ~0 removes intercept / the intercept is no longer relative to the first group 

# setting up model contrasts is more straight forward in the absence of an intercept

# and our target comparison is group, so we remove intercept from that one 

# changing column names 

colnames(design)[1:10] <- levels(sampleinfo$main_structure)

colnames(design)[11] <- levels(sampleinfo$donor)[2]

colnames(design)[12:20] <- levels(sampleinfo$samples)[2:10]

design # readable 

dim(design) # rows are samples, columns are the parameters 

# step 2: contrast matrix

# represents combinations of coefficients that I am interested in

# 9 contrasts for each tissue, each a coef

# Rows are the model parameters

# Columns are the contrast of interest 

# cbcx

contr.matrix.cbcx <- makeContrasts(
  CbCxVsCgG = CbCx - CgG, # coef 1
  CbCxVsFL = CbCx - FL, # coef 2
  CbCxVsGP = CbCx - GP, # coef 3
  CbCxVsIns = CbCx - Ins, # coef 4
  CbCxVsOL = CbCx - OL, # coef 5
  CbCxVsPHG = CbCx - PHG,# coef 6
  CbCxVsPL = CbCx - PL, # coef 7
  CbCxVsStr = CbCx - Str, # coef 8
  CbCxVsTL = CbCx - TL, # coef 9
  levels = colnames(design))

(contr.matrix.cbcx) 

# CgG

contr.matrix.cgg <- makeContrasts(
  CgGVsCbCx = CgG - CbCx, # coef 1
  CgGVsFL = CgG - FL, # coef 2
  CgGVsGP = CgG - GP, # coef 3
  CgGVsIns = CgG - Ins, # coef 4
  CgGVsOL = CgG - OL, # coef 5
  CgGVsPHG = CgG - PHG,# coef 6
  CgGVsPL = CgG - PL, # coef 7
  CgGVsStr = CgG - Str, # coef 8
  CgGVsTL = CgG - TL, # coef 9
  levels = colnames(design))

(contr.matrix.cgg) 

# FL

contr.matrix.fl <- makeContrasts(
  FLVsCbCx = FL - CbCx, # coef 1
  FLVsCgG = FL - CgG, # coef 2
  FLVsGP = FL - GP, # coef 3
  FLVsIns = FL - Ins, # coef 4
  FLVsOL = FL - OL, # coef 5
  FLVsPHG = FL - PHG,# coef 6
  FLVsPL = FL - PL, # coef 7
  FLVsStr = FL - Str, # coef 8
  FLVsTL = FL - TL, # coef 9
  levels = colnames(design))

(contr.matrix.fl) 

# GP

contr.matrix.gp <- makeContrasts(
  GPVsCbCx = GP - CbCx, # coef 1
  GPVsCgG = GP - CgG, # coef 2
  GPVsFL = GP - FL, # coef 3
  GPVsIns = GP - Ins, # coef 4
  GPVsOL = GP - OL, # coef 5
  GPVsPHG = GP - PHG,# coef 6
  GPVsPL = GP - PL, # coef 7
  GPVsStr = GP - Str, # coef 8
  GPVsTL = GP - TL, # coef 9
  levels = colnames(design))

(contr.matrix.gp) 


# INS

contr.matrix.ins <- makeContrasts(
  InsVsCbCx = Ins - CbCx, # coef 1
  InsVsCgG = Ins - CgG, # coef 2
  InsVsFL = Ins - FL, # coef 3
  InsVsGP = Ins - GP, # coef 4
  InsVsOL = Ins - OL, # coef 5
  InsVsPHG = Ins - PHG,# coef 6
  InsVsPL = Ins - PL, # coef 7
  InsVsStr = Ins - Str, # coef 8
  InsVsTL = Ins - TL, # coef 9
  levels = colnames(design))

(contr.matrix.ins) 

# OL

contr.matrix.ol <- makeContrasts(
  OLVsCbCx = OL - CbCx, # coef 1
  OLVsCgG = OL - CgG, # coef 2
  OLVsFL = OL - FL, # coef 3
  OLVsGP = OL - GP, # coef 4
  OLVsIns = OL - Ins, # coef 5
  OLVsPHG = OL - PHG,# coef 6
  OLVsPL = OL - PL, # coef 7
  OLVsStr = OL - Str, # coef 8
  OLVsTL = OL - TL, # coef 9
  levels = colnames(design))

(contr.matrix.ol) 

# PHG

contr.matrix.phg <- makeContrasts(
  PHGVsCbCx = PHG - CbCx, # coef 1
  PHGVsCgG = PHG - CgG, # coef 2
  PHGVsFL = PHG - FL, # coef 3
  PHGVsGP = PHG - GP, # coef 4
  PHGVsIns = PHG - Ins, # coef 5
  PHGVsOL = PHG - OL,# coef 6
  PHGVsPL = PHG - PL, # coef 7
  PHGVsStr = PHG - Str, # coef 8
  PHGVsTL = PHG - TL, # coef 9
  levels = colnames(design))

(contr.matrix.phg) 

# PL

contr.matrix.pl <- makeContrasts(
  PLVsCbCx = PL - CbCx, # coef 1
  PLVsCgG = PL - CgG, # coef 2
  PLVsFL = PL - FL, # coef 3
  PLVsGP = PL - GP, # coef 4
  PLVsIns = PL - Ins, # coef 5
  PLVsOL = PL - OL,# coef 6
  PLVsPHG = PL - PHG, # coef 7
  PLVsStr = PL - Str, # coef 8
  PLVsTL = PL - TL, # coef 9
  levels = colnames(design))

(contr.matrix.pl) 


# str 

contr.matrix.str <- makeContrasts(
  StrVsCbCx = Str - CbCx, # coef 1
  StrVsCgG = Str - CgG, # coef 2
  StrVsFL = Str - FL, # coef 3
  StrVsGP = Str - GP, # coef 4
  StrVsIns = Str - Ins, # coef 5
  StrVsOL = Str - OL, # coef 6
  StrVsPHG = Str - PHG,# coef 7
  StrVsPL = Str - PL, # coef 8
  StrVsTL = Str - TL, # coef 9
  levels = colnames(design))

(contr.matrix.str) 


# TL 

contr.matrix.tl <- makeContrasts(
  TLVsCbCx = TL - CbCx, # coef 1
  TLVsCgG = TL - CgG, # coef 2
  TLVsFL = TL - FL, # coef 3
  TLVsGP = TL - GP, # coef 4
  TLVsIns = TL - Ins, # coef 5
  TLVsOL = TL - OL,# coef 6
  TLVsPHG = TL - PHG, # coef 7
  TLVsPL = TL - PL, # coef 8
  TLVsStr = TL - Str, # coef 9
  levels = colnames(design))

(contr.matrix.tl) 


##############################################################################

# section 4: fitting the model 

############################################################################

# we estimate effect for different groups in a linear model y = mx +c where

# Where y is the individual row (genes) from the expression matrix
# m is the average gene expression changes that we try to estimate from our data across the replicate sample
# x is the design matrix that links the data to the coefficients we want to estimate
# C is the standard error of the least square estimates. We assume its equally distributed   

# step 1: mean variance difference: removing heteroscedascity from count data

# assumptions of linear modeling the logcpm counts:

# edgeR and deSEQ2 models the data in its natural, negative binomial distribution of raw counts 

# Negative bionomail count

# In binomial distribution, unlike normal distribution, the median is not equal to mean 

# Because the probability of finding a transcript is small but number of samples is huge.

# One way to describe is Poisson distribution but the assumption is homoscedasticity 

# What we do get in count data is Heteroscedasticity

# For a given expression level in the low range we observe a lot of variability in the variance values. 

# For the genes with low mean expression we see quite a bit of scatter.

# Therefore, NB fit is better. 

# It allows for overdispersion due to the biological variability and also a lot of low counts/ zeroes 

# but limma models the linear modeling on log-cpm counts, which has a normal distribution


# what does voom do? 

# homoscedastic distribution of the data is needed for DGE   

# voom calculates the trend by calculating the weight matrix in the regression of each and every gene (precision weights) 

par(mfrow=c(1,1))

v <- voom(x_norm, design, plot=TRUE) # voom object with a model being plotted 

v # is the voom 

# x axis is the average expression

# y axis is variation 

# the red line is the smooth lowes curve 

# quadratic mean-variance relationship (	 graph's curve forms a parabola y= ax2+bx+c)

# information from the voom plot

# 1. typically, the voom-plot shows a decreasing trend between the means and variances

# resulting from a combination of technical variation in the sequencing experiment 

# and biological variation among the replicate samples from different cell populations. 

# Experiments with high biological variation usually result in flatter trends, 

# where variance values plateau at high expression values. 

# Experiments with low biological variation tend to result in sharp decreasing trends.

# 2. visual check on the upstream filtering

# insufficient filtering == very small counts at the low end of the scale == a drop in the variance

# which means we should filter more 

# voom transforms the counts in the v$E

# v$genes == dge$genes

# v$targets == dge$samples with the norm.factors and lib.size 

# v#designs == design matrix

# v$weight == precision weights 

# voom$E == log(cpm + 0.05)  


# step 2: covariates 

# accounting for donor effect

donorcor <- duplicateCorrelation(v, design, block = sampleinfo$donor)

donorcor$consensus.correlation # 0.91

# is there a sequencing lane effect 

lanecor <- duplicateCorrelation(v, design, block = sampleinfo$samples)
  
lanecor$consensus.correlation # 0.93

# correlation within donor and sequencing lanes is strong and + tive  

# step 3: fitting linear model

# we are assuming the underlying data is normally distributed

# but we are not blocking here with duplicateCorrelation() because the blocks in design matrix

# already adjusts for difference in patient

# duplicateCorrelation() accounts for variation between what we perceive to be duplicate

# as I am comfortable with the variability within levels donor / lane 

# I will not penalize genes that show variation within donor / lane

# and move on after adding these to my model.matrix

# source: https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

# https://support.bioconductor.org/p/125489/

# https://support.bioconductor.org/p/94280/#94290

# https://support.bioconductor.org/p/86166/  

vfit <- lmFit(v, design) 

################################################################################

# section 5: looking the genes 

###############################################################################

# lmfit 

# assumes the population SD is same for all groups  

# will fit the contrast matrix 

# and fits the weight calculated in voom

# Lmfit$coefficients are the values of coefficients defined in model matrix 

# efit

# Empirical Bayes moderation is modified t test: a t test that is moderated for genomic context 

# which returns us v fit + 

# moderated t statistic, P vale and F scores

# by empirical Bayes moderation of the standard errors towards a common value.

# significance is determined by genes where p.adj < 0.05

# in the topTreat()

# logFC == log fold change between the two tissues

# AveExpr: the average log2-expression level for that gene across all samples. partifular to each gene

# t : moderated t statistic.it tells us about the difference between the two groups  

# Recall that t-statistic is the ratio 

# of difference(between the two groups) 

# and standard error of the difference. 

# SE for a random variable is sd/sqrt(length)

# SE of difference is sqrt(var(group1)/length(group1) + var(group2)/length(group2))

# P. Value : associated P value 

# adj.P.Val : P value adjusted for multiple comparison. Default method is BH by default.

# B : (a.k.a lods) is the log-odds that the gene is deferentially expressed.

# i don't understand this lods

#################################################################################


# we will fit different contrasts to same code 

cfit <- contrasts.fit(vfit, contrasts=contr.matrix.cbcx) # fitting linear model to contrasts matrix 

cfit$coefficients # values of coefs from the desing matrix 

efit <- eBayes(cfit) # Empirical Bayes moderation gets info all the genes 

summary(decideTests(efit)) # quick summary of the DEG

summary_DE <- decideTests(efit) # 0 == NotSig, 1 == up-regulated, -1 == down-regulated.

# P value plot 

hist(efit$p.value)
abline(v = 0.05, col = "red")

# most P values are significant

# t test 

htrgenes <- c("HTR1A",
              "HTR1B",
              "HTR1D",
              "HTR1E",
              "HTR1F",
              "HTR2A",
              "HTR2B",
              "HTR2C",
              "HTR4",
              "HTR5A",
              "HTR6",
              "HTR7") # my genes are the receptor genes and G alpha 

htr_dge <- list() # empty index

# loop for 9 comparisons 

for (i in 1:9) {
  treat_efit <- topTreat(efit, coef= i, n=Inf, sort.by = "none") # n = inf is all results  %>%
  
  dge_i <- treat_efit[which(rownames(treat_efit) %in% htrgenes),] # topTreat arranges genes from smallest to largest 
  
  # adjusted p-value with associated gene information, log-FC, average log-CPM, moderated t-statistic, 
  
  # raw and adjusted p-value for each gene
  
  print(colnames(efit)[i]) # checking 
  
  dge_i$diff <- as.character(colnames(efit)[i]) # the name of the comparison 
  
  htr_dge[[i]] <- dge_i # joining
  
}


htr_dge_df <- do.call(rbind, htr_dge) # combining all 9 comparison

htr_dge_df

dim(htr_dge_df)

for (i in 1:9) 
  {
  # for a specicic coef in the linear mmodel 
  
  o <- which(rownames(efit$p.value) %in% htrgenes) # the index
  
  x <- efit$Amean # x axis in MD plot 
  
  y_i <- efit$coefficients[,i] # y axis in MD plot
  
  G <- efit$genes$GENENAME # Label 
  
  # locating the genes in MD
  
  # X axis: average gene expression (Means)
  
  # Y axis: estimates coefficient for my model (log 2 fold change) (Difference) 
  
  plotMD(efit, 
         hl.col=c("turquoise", "salmon"),
         column = i, 
         coef = i,
         main = colnames(efit)[i], 
         status = summary_DE[,i], # statistical significance 
         xlim=c(-8,13))
  text(x[o], y_i[o], labels=G[o], col = "black")
  
  # volcano plot 
  
  # can locate genes with large lfc that are also significant 
  
  # X axis: (log 2 fold change) (Difference)
  
  # the more left the more downregulated 
  
  # the more right the more upregulated 
  
  # Y axis: P Value
  
  # more top is mode significant 
  
  # volcanoplot(efit,
  #             coef= i,
  #             main = colnames(efit)[i],
  #             names = efit$genes$GENENAME)
  # text(x[o], y_i[o], labels=G[o], col = "red")
  
}


################################################################################

# writing results 

# write.csv(htr_dge_df, "limma_output/1.limma_output/htrde_tl.csv")


###############################################################################
