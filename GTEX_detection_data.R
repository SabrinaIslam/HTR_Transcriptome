################################################################################

# packages 

################################################################################

# data manipulation 

library(dplyr) # for wrangling data frames 
library(tidyverse) # tidy data 

# visualisation 

library(ggplot2) # plotting 
library(gplots) # plotting data 
library(RColorBrewer) # build color-pallates for plots 
library(ggthemes) # themes 
library(hrbrthemes) # themes 
library(viridis) # colour 
library(matrixStats) # row statisitcs 


## RNA-seq specific  

library(limma) # for expression data 
library(edgeR) # for RNA- seq data 
library(stringr) # to create a matrix using split


# font issue 

windowsFonts("Arial" = windowsFont("Arial"))


################################################################################

# folder 

################################################################################


setwd("C:/Users/sabrinai/OneDrive - The University of Melbourne/PHD/Chapter2/2.GTEXRnaSeq")


################################################################################

# about the data and the experiment 

################################################################################

# source: https://gtexportal.org/home/datasets 

# details on https://gtexportal.org/home/tissueSummaryPage

# tissue: 55 human tissues 
# replicates: none   
# libraries: 17382
# sequencing:  Illumina TrueSeq HiSeq 2000 and 2500 RNA sequencing to obtain 76bp paired-end reads
# median coverage was ~82M total reads
# aligned to the human genome hg38/GRCh38 human genome using STAR v2.5.3a
# transcripts defined using the the GENCODE 26 transcriptome


################################################################################

# Objective of this code 

################################################################################

# WHOLE GTEX LIBRARY PROB 
        
################################################################################

# DATA READYING 

################################################################################

# step 1: filetring the names of the tissues

meta <- read.csv("Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                 header = T, 
                 sep = "\t",     
                 quote = "")

tissue_list <- levels(as.factor(meta$SMTSD))

tissue_list # 55 tissue 

tissue_list <- tissue_list[-c(24)]  

tissue_list # 54 after remvong cancer cells

# step 2: preparing the trasnformed data 

htr_count <- read.csv("Data/GTEX_HTR.csv", 
                      row.names = 1,
                      header = T)

rownames(htr_count) # checks 

htr_exp <- (htr_count[, c(3:17384)]) # just the count

htr_annot <- (htr_count[, c(1:2)]) # annotation 

head(htr_annot) # checks 

htr_dge <- DGEList(counts = htr_exp, genes = htr_annot) # DGE 

log2_htr_dge <- log2((htr_dge$counts) + 1) # log

# step 3: names of the samples 

gtex_tissue <- read.csv("Data/sampleTissues.csv",
                header = TRUE,
                check.names = F,
                row.names = 1) # sample annotations 


colnames(gtex_tissue)

gtex_tissue$SAMPID <- gsub("\\-", ".", (gtex_tissue$SAMPID)) # changing "." to "_"

head(gtex_tissue$SAMPID) == head(colnames(log2_htr_dge)) # T 

gtex_tissue %>%
  filter(SMTSD %in% tissue_list) %>%
  select(SAMPID) %>% unlist # testing 

rownames(htr_dge$counts) # checks 

rownames(log2_htr_dge) # checks 


################################################################################

# MEDIAN, MEAN AND PERCENATGE OF HTR EXPRESSION ACROSS TISSUES 

################################################################################

# empty index 

exp_prob_list <- list()

exp_median_list <- list()

exp_mean_list <- list()

# the for loop   

# remember to change to 1: length (tissue_list)

for (i in 1:length (tissue_list)) {
        
        # step 1 selecting the samples from each tissue 
  
        tissue_i <- gtex_tissue %>%
        filter (SMTSD %in% tissue_list[i]) %>%
        dplyr::select(SAMPID) %>% unlist 
        
        print(tissue_i)

        # step 2 divide

        htr_exp_i <- data.frame(log2_htr_dge[, colnames(log2_htr_dge) %in%  tissue_i])
        
        # step 3: making a percetange of of the HTR tissues being present in the samples of each tissue

        sample_num_i <- ncol(htr_exp_i) # for division

        bin_htr_tissue_i <- data.frame(ifelse(htr_exp_i[,1:sample_num_i] > 1.15, 1, 0)) # binarising matrix

        htr_probability_i <- (data.frame(rowSums(bin_htr_tissue_i)/sample_num_i)*100) # number of 1s

        names(htr_probability_i)[names(htr_probability_i) == 'rowSums.bin_htr_tissue_i..sample_num_i'] <-
                sprintf( "%s", tissue_list[i])

        exp_prob_list[[i]] <- htr_probability_i # adding to the empty list
        
        # step 4 making a median of of the HTR tissues being present in the samples of each tissue

        median_htr_tissue_i <- data.frame(rowMedians(as.matrix(htr_exp_i)))

        names(median_htr_tissue_i)[names(median_htr_tissue_i) == 'rowMedians.as.matrix.htr_exp_i..'] <-
                sprintf( "%s", tissue_list[i])

        exp_median_list[[i]] <- median_htr_tissue_i # adding to the empty list

        # step 5 making a mean of of the HTR tissues being present in the samples of each tissue

        mean_htr_tissue_i <- data.frame(apply((as.matrix(htr_exp_i)), 1, mean, trim = 0.25))

        names(mean_htr_tissue_i)[names(mean_htr_tissue_i) == 'apply..as.matrix.htr_exp_i....1..mean..trim...0.25.'] <-sprintf( "%s", tissue_list[i])

        exp_mean_list[[i]] <- mean_htr_tissue_i # adding to the empty list
}

################################################################################

# mean data 

htr_mean <- do.call(cbind, exp_mean_list) # binding all mean values 

colnames(htr_mean)

rownames(htr_mean) <- c("ENSG00000158748.3" = "HTR6",
                     "ENSG00000179546.4" = "HTR1D",
                     "ENSG00000135914.5" = "HTR2B",
                     "ENSG00000179097.5" = "HTR1F",
                     "ENSG00000178394.4" = "HTR1A",
                     "ENSG00000164270.17" = "HTR4",
                     "ENSG00000135312.6" = "HTR1B",
                     "ENSG00000168830.7" = "HTR1E",
                     "ENSG00000157219.3" = "HTR5A",
                     "ENSG00000148680.15" = "HTR7",
                     "ENSG00000102468.10" = "HTR2A",
                     "ENSG00000147246.9" = "HTR2C")


htr_mean$Receptor <- rownames(htr_mean) # renaming by receptor 

long_htr_mean <- htr_mean %>% gather(key = "Tissue", value = "Mean", -"Receptor") 

long_htr_mean$Mean <- as.numeric(long_htr_mean$Mean) # "Makes value numeric"

long_htr_mean$Tissue <- as.factor(long_htr_mean$Tissue) # should be factored 

long_htr_mean$Receptor <- as.factor(long_htr_mean$Receptor)

#################################################################################

# median data 

htr_m <- do.call(cbind, exp_median_list) # binding all prob values 

rownames(htr_m) <- c("ENSG00000158748.3" = "HTR6",
                     "ENSG00000179546.4" = "HTR1D",
                     "ENSG00000135914.5" = "HTR2B",
                     "ENSG00000179097.5" = "HTR1F",
                     "ENSG00000178394.4" = "HTR1A",
                     "ENSG00000164270.17" = "HTR4",
                     "ENSG00000135312.6" = "HTR1B",
                     "ENSG00000168830.7" = "HTR1E",
                     "ENSG00000157219.3" = "HTR5A",
                     "ENSG00000148680.15" = "HTR7",
                     "ENSG00000102468.10" = "HTR2A",
                     "ENSG00000147246.9" = "HTR2C")


htr_m$Receptor <- rownames(htr_m) # renaming by receptor 

long_htr_m <- htr_m %>% gather(key = "Tissue", value = "Median", -"Receptor") 

long_htr_m$Median <- as.numeric(long_htr_m$Median) # "Makes value numeric"

long_htr_m$Tissue <- as.factor(long_htr_m$Tissue) # should be factored 

long_htr_m$Receptor <- as.factor(long_htr_m$Receptor)

#################################################################################

# probability data  

htr_pro <- do.call(cbind, exp_prob_list) # binding all prob values 

rownames(htr_pro) <- c("ENSG00000158748.3" = "HTR6",
                       "ENSG00000179546.4" = "HTR1D",
                       "ENSG00000135914.5" = "HTR2B",
                       "ENSG00000179097.5" = "HTR1F",
                       "ENSG00000178394.4" = "HTR1A",
                       "ENSG00000164270.17" = "HTR4",
                       "ENSG00000135312.6" = "HTR1B",
                       "ENSG00000168830.7" = "HTR1E",
                       "ENSG00000157219.3" = "HTR5A",
                       "ENSG00000148680.15" = "HTR7",
                       "ENSG00000102468.10" = "HTR2A",
                       "ENSG00000147246.9" = "HTR2C")


htr_pro$Receptor <- rownames(htr_pro) # renaming by receptor 

long_htr_p <- htr_pro %>% gather(key = "Tissue", value = "Probability", -"Receptor") 

long_htr_p$Probability <- as.numeric(long_htr_p$Probability) # "Makes value numeric"

long_htr_p$Tissue <- as.factor(long_htr_p$Tissue) # should be factored 

long_htr_p$Receptor <- as.factor(long_htr_p$Receptor) # factored 


################################################################################

# visualisatin 

################################################################################

probplot <- ggplot(long_htr_p,
                   aes(x= Receptor, 
                       y=Tissue, 
                       fill = Probability)) +
  geom_tile(aes(fill = Probability),
            colour = "white",
            size = 0.5) + # making the tile size nicer 
  scale_y_discrete(limits = rev) +
  scale_fill_gradient(low= "white", high="#0e7f92", guide="colorbar") + # changing color gradient
  ggtitle("Detection % of HTR in 17382 GTEx samples") + # adding title 
  theme_tufte() +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10))  # "TPM" = size 8

probplot


meanplot <- ggplot(long_htr_mean,
                aes(x= Receptor, 
                    y=Tissue, 
                    fill = Mean)) +
  geom_tile(aes(fill = Mean),
            colour = "white",
            size = 0.5) + # making the tile size nicer 
  scale_y_discrete(limits = rev) +
  scale_fill_gradient(low= "white", high="#bf0016", guide="colorbar") + # changing color gradient
  ggtitle("Mean Expression of HTR in Tissues") + # adding title 
  theme_tufte() +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10))  # "TPM" = size 8

meanplot


mplot <- ggplot(long_htr_m,
                aes(x= Receptor, 
                    y=Tissue, 
                    fill = Median)) +
  geom_tile(aes(fill = Median),
            colour = "white",
            size = 0.5) + # making the tile size nicer 
  scale_y_discrete(limits = rev) +
  scale_fill_gradient(low= "white", high="#bf0016", guide="colorbar") + # changing color gradient
  ggtitle("Median Expression of HTR in Tissues") + # adding title 
  theme_tufte() +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10))  # "TPM" = size 8

mplot

    
        
#################################################################################
        
# writing csv for hclust 

#################################################################################

write.csv(htr_pro, "Data/HTR_Prob.csv", sep = "\t")
        
write.csv(htr_m, "Data/HTR_Median.csv", sep = "\t")
        
        write.csv(htr_mean, "Data/HTR_Mean.csv", sep = "\t")
        