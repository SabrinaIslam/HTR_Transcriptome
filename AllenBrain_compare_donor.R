################################################################################

# packages 

################################################################################

# working with data 

library(dplyr) # for wrangling data frames 
library(tidyverse) # tidy data 
library(ggpubr) # putting images together 

# visualisation 

library(ggplot2) # plotting 
library(gplots) # plotting data 
library(RColorBrewer) # build color-pallates for plots 
library(ggthemes) # themes 

# anova and related 

library(rstatix)
library(reshape)
library(plyr)
library(datarium)
library(paletteer)


# statistics 

library(spgs) # tuning point test 
library(matrixStats) # calculating matrix statistics 

# RNA-seq specific packages  

library(limma) # for expression data 
library(edgeR) # for RNA- seq data 
library(affy) # plotDensity polynomial fitted plots  
library(org.Hs.eg.db) # human annotations 


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


################################################################################

# Objective of this code 

################################################################################

# I will compare the two donors to evaulate if I can treat them as replicates 

# I will visulaise libraires from both donors side by side

# I will check which tissue regions tend to be more variable amoong samples

# I will perform correlation on all counts from both donors 

################################################################################

# section 1: import the data 

################################################################################

# step 1: import Raw counts of genes

#################################################################################

# librarY 1

ab_count_1 <- read.csv("Data/RNAseqCounts.csv",
                       header = F,
                       check.names = F,
                       row.names = 1)

# samples annotation 

ab_samples_1 <- read.csv("Data/SampleAnnot.csv",
                         header = TRUE,
                         check.names = F,
                         row.names = 1)

ab_tissues_1 <- data.frame(ab_samples_1[ , c(6,7)]) # extracting the brain strucutres 

ab_tissues_1$SAMPID <- rownames(ab_tissues_1) # adding the sample id

str(ab_tissues_1) # all columns are character vectors  

# changing them to factor 

ab_tissues_1$sub_structure <- as.factor(ab_tissues_1$sub_structure)

ab_tissues_1$main_structure <- as.factor(ab_tissues_1$main_structure)

# make columns match

colnames(ab_count_1) <-  rownames(ab_tissues_1) 

# preparing annotation

ab_annot_1 <- data.frame(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = rownames(ab_count_1),
                                               columns = c ("GENENAME"),
                                               keytype="GENENAME"))

length(unique(ab_annot_1$GENENAME)) # no duplicates 

nrow(ab_annot_1) == nrow(ab_count_1) # T

# librarY 2

ab_count_2 <- read.csv("Data2/RNAseqCounts.csv",
                       header = F,
                       check.names = F,
                       row.names = 1)

# samples annotation 

ab_samples_2 <- read.csv("Data2/SampleAnnot.csv",
                         header = TRUE,
                         check.names = F,
                         row.names = 1)

ab_tissues_2 <- data.frame(ab_samples_2[ , c(6,7)]) # extracting brain strucutre 

ab_tissues_2$SAMPID <- rownames(ab_tissues_2) # adding sample information 

str(ab_tissues_2) # again, character 

# making them factors 

ab_tissues_2$sub_structure <- as.factor(ab_tissues_2$sub_structure)

ab_tissues_2$main_structure <- as.factor(ab_tissues_2$main_structure)

# make columns match

colnames(ab_count_2) <-  rownames(ab_tissues_2) 

# preparing annotation

ab_annot_2 <- data.frame(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = rownames(ab_count_2),
                                               columns = c ("GENENAME"),
                                               keytype="GENENAME"))


length(unique(ab_annot_2$GENENAME)) # no duplicates 

nrow(ab_annot_2) == nrow(ab_count_2) # T 

#################################################################################

# step 2: dge object 

# library 1

ab_dge_1 <- DGEList(counts = ab_count_1, genes = ab_annot_1)

names(ab_dge_1) # sample, counts, and genes 

log2_ab_dge_1 <- log2(cpm(ab_dge_1$counts) + 1) # logcpm transformation

# library 2

ab_dge_2 <- DGEList(counts = ab_count_2, genes = ab_annot_2)

names(ab_dge_2) # sample, counts, and genes 

log2_ab_dge_2 <- log2(cpm(ab_dge_2$counts) + 1) # logcpm transformation 

#################################################################################

# step 3: preparing the tissue data 

## means of each tissue 

### sample 1

tissue_list <- levels(as.factor(ab_tissues_1$sub_structure)) # tissue list 

mean_list_1 <- list() # empty list 

for (i in 1:length (tissue_list)) # remember to change to 1: length (tissue_list)
{
  # step 1 selecting the samples from each tissue 
  
  tissue_i <- ab_tissues_1 %>%
    filter (ab_tissues_1$sub_structure == tissue_list[i]) %>%
    dplyr::select(SAMPID) %>% unlist()
  
  print(tissue_i)
  
  # step 2 dge
  
  i_dge <- (log2_ab_dge_1[, colnames(log2_ab_dge_1) %in%  tissue_i])
   
  print(head(i_dge))
   
  # step 3 making a mean of of the HTR tissues being present in the samples of each tissue
  
  mean_tissue_i <- data.frame(rowMeans(as.matrix(i_dge)))

  names(mean_tissue_i)[names(mean_tissue_i) == 'rowMeans.as.matrix.i_dge..'] <-
    sprintf( "%s", tissue_list[i])

  mean_list_1[[i]] <- mean_tissue_i # adding to the empty list
}


tissue_mean_1 <- do.call(cbind, mean_list_1) # binding all values 


### sample 2

tissue_list <- levels(as.factor(ab_tissues_2$sub_structure)) # tissue list 

mean_list_2 <- list() # empty list 

for (i in 1:length (tissue_list)) # remember to change to 1: length (tissue_list)
{
  # step 1 selecting the samples from each tissue 
  
  tissue_i <- ab_tissues_2 %>%
    filter (ab_tissues_2$sub_structure == tissue_list[i]) %>%
    dplyr::select(SAMPID) %>% unlist()
  
  print(tissue_list[i])
  
# step 2 dge
  
  i_dge <- data.frame(log2_ab_dge_2[, colnames(log2_ab_dge_2) %in%  tissue_i])
  
  # step 3 making a mean of of the HTR tissues being present in the samples of each tissue
  
  mean_tissue_i <- data.frame(rowMeans(as.matrix(i_dge)))
  
  names(mean_tissue_i)[names(mean_tissue_i) == 'rowMeans.as.matrix.i_dge..'] <-
    sprintf( "%s", tissue_list[i])
  
  mean_list_2[[i]] <- mean_tissue_i # adding to the empty list
}


tissue_mean_2 <- do.call(cbind, mean_list_2) # binding all values 

tissue_mean_2

# do columns from means from each donor (tissue) match 

colnames(tissue_mean_1) %in% colnames(tissue_mean_2) # T

#################################################################################

# section 2: compare two donors  

#################################################################################

# step 1:  making long data 

### sample 1

long_mean_1 <- tissue_mean_1 %>%
  gather(key = "Tissue", value = "Mean") %>%
  convert_as_factor(Tissue) 

long_mean_1 <- cbind(long_mean_1, Sample = "H0351.2001") # add tag for sample

### sample 2

long_mean_2 <- tissue_mean_2 %>%
  gather(key = "Tissue", value = "Mean") %>%
  convert_as_factor(Tissue) 

long_mean_2 <- cbind(long_mean_2, Sample = "H0351.2002") # add tag for sample

# combine 

all_mean <- rbind(long_mean_1, long_mean_2)

head(all_mean, 3)

#################################################################################

# step 2: summarise 

summary <- all_mean %>%
  group_by(Sample) %>%
  get_summary_stats(Mean, type = "mean_sd")

data.frame(summary)

# exact same n, mean and sd are very close 

summary <- all_mean %>%
  group_by(Sample) %>%
  get_summary_stats(Mean, type = "full")

data.frame(summary)

# min, max, iqr, mean, sd are very close and se and ci are the same 

#################################################################################

# step 3: simple statistics  

# means and medians of the donor 

median(log2_ab_dge_1)  # 3.11

median(log2_ab_dge_2) # 3.10

mean(log2_ab_dge_1) # 3.11

mean(log2_ab_dge_2) # 3.15

# quantiles    

donor_1_qq <-quantile(log2_ab_dge_1, probs = c(0.05,0.25, 0.5, 0.75, 0.95)) 

donor_2_qq <- quantile(log2_ab_dge_2, probs = c(0.05,0.25, 0.5, 0.75, 0.95))  

donor_stats <- rbind(donor_1_qq, donor_2_qq)

donor_stats

#################################################################################

# step 4: visualising 

all_mean %>%
ggplot() +
  geom_boxplot(aes(x = Tissue, y = Mean, fill = Sample)) + 
  labs(title = "Count of the two donors",
           x = "Samples", # by donor 
           y = "log2 cpm") +
      scale_fill_manual(values = c("#0072B2", "#D55E00")) +  # color by donor
      theme(axis.line = element_line(),  
            axis.text.x = element_text(angle = 90, vjust=0.6, hjust=0),
            panel.background = element_blank())

# visually, Caudate, Putamen, Str_v1, GP are tissues where the donor vary  

################################################################################

# section 3: what are the variable tissues 

################################################################################

# step 1: compute variation in eahc tissue  

log2_ab_dge_whole <- cbind (log2_ab_dge_1, log2_ab_dge_2) # combine the 2 dataframe 

ab_tissues <- rbind(ab_tissues_1, ab_tissues_2) # and tissues 

tissue_list <- levels(as.factor(ab_tissues$sub_structure))

tissue_list

var_list <- list() # empty list 

for (i in 1:length (tissue_list)) # remember to change to 1: length (tissue_list)
{
  # step 1 selecting the samples from each tissue 
  
  tissue_i <- ab_tissues %>%
    filter (ab_tissues$sub_structure %in% tissue_list[i]) %>%
    dplyr::select(SAMPID) %>% unlist()
  
  print(tissue_i)
  
  # step 2 dge

  i_dge <- (log2_ab_dge_whole[, colnames(log2_ab_dge_whole) %in%  tissue_i])

  print(head(i_dge))

  # step 3 variance for each gene 

  var_tissue_i <- data.frame(rowVars(as.matrix(i_dge)))

  names(var_tissue_i)[names(var_tissue_i) == 'rowVars.as.matrix.i_dge..'] <-
    sprintf( "%s", tissue_list[i])

  var_list[[i]] <- var_tissue_i # adding to the empty list
}

tissue_var <- do.call(cbind, var_list) # binding all values 

vars <- data.frame(apply(tissue_var, 2, mean)) # mean variance of all genes for each column

vars <- rename(vars, c("apply.tissue_var..2..mean." = "variance"))

#################################################################################

# step 2: density 

sample_density <- density(vars$variance) # density 

sample_delta <-  diff(sample_density$y) # differences 

sample_turns <- which(sample_delta[-1] * sample_delta[-length(sample_delta)] < 0) + 1 # turning points 

plot(sample_density,
     lty=1,
     col= "#0072B2",
     xlab="log2(cpm+1)",
     ylab = "density",
     lwd=2,
     ylim=c(0, 30),
     main=sprintf( "Variance"))
points(sample_density$x[sample_turns], 
       sample_density$y[sample_turns], 
       pch=16, col="red")
abline(v = sample_density$x[sample_turns][2], col = "red")
grid()

sample_density$x[sample_turns][2] # is 0.095

tissues_higher_var <- vars %>%
  filter(vars$variance > 0.095)

tissues_higher_var # caudate, GP, Pest_V2, PHG, Putamen, Str_V1

###############################################################################

# section 4: compare between the two donors

###############################################################################

par(mfrow = c(1,1))

# step 1: both donor density plot 

plot(density(log2_ab_dge_1),
     lty=1,
     col= "#0072B2",
     xlab="log2(cpm+1)",
     ylab = "density",
     lwd=2,
     ylim=c(0,0.8),
     main = "Density of the two donors")

lines(density(log2_ab_dge_2),
      lty=1,
      col= "#D55E00")
abline(v = 1.7)
legend("topright",    
       c("H0351.2001", "H0351.2002"),
       col = c("#0072B2", "#D55E00"),
       lty = 1,
       lwd = 2)
axis(side =1, 
     at = seq(0,15,1), 
     labels = seq(0,15,1))
grid()

#################################################################################

# step 2: correlation

# donor 1 long data 

donor_1 <- as.data.frame(log2_ab_dge_1)

donor_1$Genes <- rownames(donor_1)

donor_1 <- donor_1 %>%
  gather(key = "sample", value = "Count", -"Genes")

donor_1$Donor <- c("H0351.2001")

# donor_2 long data 

donor_2 <- as.data.frame(log2_ab_dge_2)

donor_2$Genes <- rownames(donor_2)

donor_2 <- donor_2 %>%
  gather(key = "sample", value = "Count", -"Genes")

donor_2$Donor <- c("H0351.2002")

comb <- rbind(donor_1, donor_2)

summary <- comb %>%
  group_by(Donor) %>%
  get_summary_stats(Count, type = "full")

data.frame(summary)

summary <- data.frame(t(summary))

names(summary) <- c("H0351.2001", "H0351.2002")

summary <- summary[-c(1,2), ]


##################################################


donor_1 <- donor_1[ , c(3)] 

donor_2 <- donor_2[ , c(3)] 


# pearson 

cor_donor_p <- cor(donor_1, donor_2, method = "pearson") 

cor_donor_p # 0.97

# spearman 

cor_donor_s <- cor(donor_1, donor_2, method = "spearman") 

cor_donor_s # 0.97

#############################################################################

# outcome:

# the key sample statistics for both donors are very close

# boxplots of library sizes are very similar 

# the correlation between them is 97%

# density plots align 

# conclusions: we can combine them in a single data frame 
