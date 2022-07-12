################################################################################

# packages 

################################################################################

# working with data 

library(dplyr) # for wrangling data frames 
library(tidyverse) # tidy data 


# statistics

library(matrixStats) # calculating matrix statistics 
library(rstatix) # for converting to factor 

# visualisation 

library(ggplot2) # plotting 
library(RColorBrewer) # build color-palates for plots 
library(ggthemes) # themes 


# RNA-seq specific packages

library(limma) # for expression data 
library(edgeR) # for RNA- seq data 
library(org.Hs.eg.db) # annotation 

# font issue 

windowsFonts("Arial" = windowsFont("Arial"))

# pca

library(gplots) # plots 
library(ggfortify) # pca

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

# I will look at the relative abundance of the HTR in brain tissues

# I will show that HTR2B and HTR1F can not be reliably detected in the CNS 

################################################################################

# section 1: import the data 

################################################################################

# importing raw counts 

ab_norm_counts <- read.csv("Data/NormalisedLibraries.csv",
                           header = T,
                           check.names = F,
                           row.names = 1)


# preparing annotation

ab_norm_annot <- read.csv("Data/NormalisedGenes.csv",
                          header = T,
                          check.names = F,
                          row.names = 1)

# dge object 

ab_dge_whole <- DGEList(counts = ab_norm_counts, genes = ab_norm_annot)

dim(ab_dge_whole$counts)

# logcpm transformation 

log2_ab_dge_whole <- log2(cpm(ab_dge_whole$counts) + 1)


# combine

# tissues from donor 1

ab_samples_1 <- read.csv("Data/SampleAnnot.csv",
                         header = TRUE,
                         check.names = F,
                         row.names = 1)

ab_tissues_1 <- data.frame(ab_samples_1[ , c(6,7)]) # main and sub brain structures

ab_tissues_1$SAMPID <- rownames(ab_tissues_1) # add sampid

# factorize 

ab_tissues_1$sub_structure <- as.factor(ab_tissues_1$sub_structure) 

ab_tissues_1$main_structure <- as.factor(ab_tissues_1$main_structure) 

# tissues from donor 2

ab_samples_2 <- read.csv("Data2/SampleAnnot.csv",
                         header = TRUE,
                         check.names = F,
                         row.names = 1)

ab_tissues_2 <- data.frame(ab_samples_2[ , c(6,7)]) # main and sub brain structures

ab_tissues_2$SAMPID <- rownames(ab_tissues_2) # add sampid

# factorize 

ab_tissues_2$sub_structure <- as.factor(ab_tissues_2$sub_structure)

ab_tissues_2$main_structure <- as.factor(ab_tissues_2$main_structure)

ab_tissues <- rbind(ab_tissues_1, ab_tissues_2)

################################################################################

# section 2: selecting just HTR 

################################################################################

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
              "HTR7") # are my genes 

table(rownames(log2_ab_dge_whole) %in% htrgenes) # all 11 genes are found except HTR2B 

htr_ind <- rownames(log2_ab_dge_whole) %in% htrgenes # indexing the genes

colnames(log2_ab_dge_whole) # samples 

htr_dge <- as.data.frame(log2_ab_dge_whole %>%
                           subset(rownames(log2_ab_dge_whole) %in%  htrgenes)) # finding from df

# comparison 

htr_dge

main_str <- levels(ab_tissues$main_structure)

main_str # 10 main brain structures 


sub_str <- levels(ab_tissues$sub_structure)

sub_str # 29 main brain structures 


tissue_list <- levels(as.factor(ab_tissues$sub_structure))

tissue_list # a list of tissues 

tissue_cns_i <- ab_tissues %>%
  filter(sub_structure %in%  tissue_list) %>%
  dplyr::select(SAMPID) %>% unlist()

tissue_cns_i # a list of samples 


################################################################################

# section 4: comparing main_structures 

################################################################################

main_col <- c("CgG" = "#9970ab",
              "Ins" = "#c51b7d",
              "FL" = "#f46d43",
              "PL" = "#66c2a5" ,
              "TL" = "#5e4fa2",
              "PHG" = "#41ab5d",
              "OL" = "#e6f598",
              "Str" = "#6baed6",
              "GP" = "#fdae61",
              "CbCx" = "#9e0142")

for (i in 1:11)
{
  long_i <- as.data.frame(htr_dge[i, colnames(htr_dge) %in%  tissue_cns_i]) %>%
    gather(key = "Sample", value = "Value" ) 
  
  long_i <- cbind(long_i,
                  ab_tissues %>%
                  filter(ab_tissues$SAMPID %in% long_i$Sample) %>%
                  dplyr::select(main_structure)) 
  
  
  long_i$main_structure <- as.factor(long_i$main_structure)
  
  htr_i <- rownames (htr_dge[i,] )
  
  structure_order <- factor(long_i$main_structure, levels = c("CgG", "Ins", "FL", "PL", "TL", 
                                                              "PHG", "OL", 
                                                              "Str", "GP", "CbCx"))
  
  box_i <- long_i %>%
    ggplot() +
    geom_violin(aes(x = structure_order, y = Value, fill = main_structure,col = main_structure)) +
    geom_boxplot(aes(x = structure_order, y = Value, col = main_structure), width = 0.05) +
    scale_y_continuous(breaks =seq(0, 6, .5), limit = c(0, 6)) +
    labs(title = sprintf( "%s", htr_i),
         x = "Tissue",
         y = "log2(cpm+1)") +
    scale_fill_manual(values = main_col) +
    scale_color_manual(values = main_col) +
    theme_tufte()+
    theme(text = element_text(family = "Arial"),
          axis.line = element_line(),
          axis.text.x = element_text(size = 10, angle = 90),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position ="none",
          plot.title = element_text(hjust = 0.5))+ 
    geom_hline(yintercept = 1.3, size = 1, linetype ="dashed", color = "grey") + 
    # geom_hline(yintercept = quantile(log2_ab_dge_whole, probs = 0.25 ), size = 1, linetype ="dashed", color = "green") +
    # geom_hline(yintercept = quantile(log2_ab_dge_whole, probs = 0.75 ), size = 1, linetype ="dashed", color = "green") +
    geom_hline(yintercept = median(log2_ab_dge_whole), size = 1, linetype ="dashed", color = "grey") 
    
  print(box_i) 
  }


################################################################################

# summary of main structures 

#################################################################################

rownames(htr_dge)

# empty frame 

htr_stat <- list()

for (i in 1:11) {
  
  long_i <- as.data.frame(htr_dge[i, colnames(htr_dge) %in%  tissue_cns_i]) %>%
    gather(key = "Sample", value = "Value") 
  
  long_i <- cbind(long_i,
                  ab_tissues %>%
                    filter(ab_tissues$SAMPID %in% long_i$Sample) %>%
                    dplyr::select(main_structure)) 
  
  summary_i <- data.frame(long_i %>%
                          group_by(main_structure) %>%
                          get_summary_stats(Value, type = "full"))


  summary_i <- summary_i[, -c(2)]
  
  summary_i$Receptor <- rownames (htr_dge[i,] )
  
  summary_i <- summary_i[, c(14, 1:13)]

  htr_stat[[i]] <- summary_i
}

htr_summary_stats <- do.call(rbind, htr_stat)

write.csv(htr_summary_stats, "Output/htr_summary.csv")

#################################################################################

# section 5: an objective standard

#################################################################################

quantile(log2_ab_dge_whole, probs = c(0.25, 0.5, 0.75))

# 1.3 - 2.75 is low

# 2.75 - 4.44 is moderate 

# 4.44- 5.89 is high


##############################################################################

# Striatum 

###############################################################################

Bgl <- c("Str", "GP")

bgl_col <- c("HTR1D" = "#8da0cb",
             "HTR1E" = "#8da0cb",
             "HTR5A" = "#8da0cb",
             "HTR2C" = "#66c2a5",
             "HTR4" = "#fc8d62", 
             "HTR6" = "#fc8d62")


for (i in 1:2) {
  
  tissue_bgl_i <- ab_tissues %>%
    filter (ab_tissues$main_structure %in% Bgl[i]) %>%
    dplyr::select(SAMPID) %>% unlist()
  
  median_i <- median(as.numeric(log2_ab_dge_whole[, colnames(log2_ab_dge_whole) %in%  tissue_bgl_i]))
  
  htr_dge_bgl_i <- as.data.frame(htr_df[, colnames(htr_df) %in%  tissue_bgl_i])
  
  htr_dge_bgl_i <- htr_dge_bgl_i[-c(1, 2, 5, 10), ]
  
  rownames(htr_dge_bgl_i)
  
  htr_dge_bgl_i$Receptor <- rownames(htr_dge_bgl_i)
  
  long_dge_bgl_i <- htr_dge_bgl_i %>%
    gather(key = "Sample", value = "Value", -"Receptor") %>%
    convert_as_factor(Sample) %>%
    convert_as_factor(Receptor)
  
  dge_order <- factor(long_dge_bgl_i$Receptor, levels = c("HTR1E","HTR5A","HTR2C", "HTR1D", "HTR4", "HTR6" ))
  
  expplot_bgl_i <- long_dge_bgl_i %>%
    ggplot() +
    geom_violin(aes(x = dge_order, y = Value, fill = Receptor, col = Receptor)) +
    geom_boxplot(aes(x = dge_order, y = Value, col = Receptor), width = 0.05) +
    scale_y_continuous(breaks =seq(0, 6, .5), limit = c(0, 6)) +
    labs(title = sprintf( "%s", Bgl[i]),
         x = "Receptor",
         y = "log2 cPM") +
    scale_fill_manual(values = bgl_col) +
    scale_color_manual(values = bgl_col) +
    theme_tufte() +
    theme(text = element_text(family = "Arial"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.line = element_line(),
          legend.position = "none",
          axis.title.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 1.3, linetype ="dashed", size = 1, color = "grey") +
    geom_hline(yintercept = median_i,size = 1,linetype ="dashed", color = "grey")
  
  print(expplot_bgl_i)
  
}

################################################################################

# tissue 

################################################################################

receptor_col <- c("HTR1A" = "#8da0cb", 
                  "HTR1B" = "#8da0cb",
                  "HTR1D" = "#8da0cb",
                  "HTR1E" = "#8da0cb",
                  "HTR5A" = "#8da0cb",
                  "HTR2A" = "#66c2a5", 
                  "HTR2C" = "#66c2a5",
                  "HTR4" = "#fc8d62", 
                  "HTR6" = "#fc8d62", 
                  "HTR7" = "#fc8d62")


htr_df <- htr_dge[-c(5), ] # removing htr1f 

rownames(htr_df) # 10 receptors 


for (i in 1:10) {
tissue_i <- ab_tissues %>%
  filter (ab_tissues$main_structure == main_str[i]) %>%
  dplyr::select(SAMPID) %>% unlist()

median_i <- median(as.numeric(log2_ab_dge_whole[, colnames(log2_ab_dge_whole) %in%  tissue_i]))

htr_dge_i <- as.data.frame(htr_df[, colnames(htr_df) %in%  tissue_i])

rownames(htr_dge_i)

htr_dge_i$Receptor <- rownames(htr_dge_i)

long_dge_i <- htr_dge_i %>%
  gather(key = "Sample", value = "Value", -"Receptor") %>%
  convert_as_factor(Sample) %>%
  convert_as_factor(Receptor)
#
dge_order <- factor(long_dge_i$Receptor, levels = c("HTR1E","HTR5A","HTR2A","HTR2C", "HTR1B", "HTR1D",
                                                    "HTR1A", "HTR4", "HTR6", "HTR7"))

expplot_i <- long_dge_i %>%
  ggplot() +
  geom_violin(aes(x = dge_order, y = Value, fill = Receptor, col = Receptor)) +
  geom_boxplot(aes(x = dge_order, y = Value, col = Receptor), width = 0.05) +
  scale_y_continuous(breaks =seq(0, 6, .5), limit = c(0, 6)) +
  labs(title = sprintf( "%s", main_str[i]),
       x = "Receptor",
       y = "log2 cPM") +
  scale_fill_manual(values = receptor_col) +
  scale_color_manual(values = receptor_col) +
  theme_tufte() +
  theme(text = element_text(family = "Arial"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(),
        legend.position = "none",
        axis.title.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 1.3, linetype ="dashed", size = 1,  color = "grey") +
  geom_hline(yintercept = median_i, linetype ="dashed", size = 1, color = "grey")

print(expplot_i)
}

#################################################################################

# section 6: htr1f and htr2b are not detected proof 

#################################################################################

# data will be unfiltered to avoid filtering of these two lowly expressed genes 


ab_count_1 <- read.csv("Data/RNAseqCounts.csv",
                       header = F,
                       check.names = F,
                       row.names = 1)

colnames(ab_count_1) <-  rownames(ab_tissues_1) 


ab_count_2 <- read.csv("Data2/RNAseqCounts.csv",
                       header = F,
                       check.names = F,
                       row.names = 1)


colnames(ab_count_2) <-  rownames(ab_tissues_2) 


ab_count_whole <- cbind(ab_count_1, ab_count_2)


ab_annot_1 <- data.frame(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = rownames(ab_count_1),
                                               columns = c ("GENENAME"),
                                               keytype="GENENAME"))


ab_annot_2 <- data.frame(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys = rownames(ab_count_2),
                                               columns = c ("GENENAME"),
                                               keytype="GENENAME"))


ab_annot_whole <- ab_annot_1 # for DGE

x <- DGEList(counts = ab_count_whole, genes = ab_annot_whole)

x_norm <- calcNormFactors(x) # normalising (but no filtering)

lcpm_norm <- log2(cpm(x_norm$counts) + 1) # log2 + cpm transformation 

exp_median_list <- list()

htr_median_list <- list()

htr_iqr_list <- list()


for (i in 1:11) {
  # step 1 selecting the samples from each tissue 
  
  tissue_i <- ab_tissues %>%
    filter (ab_tissues$sub_structure == tissue_list[i]) %>%
    dplyr::select(SAMPID) %>% unlist()
  
  print(tissue_list[i])
  
  # step 2 divide 
  
  i_dge <- data.frame(lcpm_norm[, colnames(lcpm_norm) %in%  tissue_i])
  
  print(i_dge)
  
  htr_i <- rownames(i_dge[htrgenes,]) # indexing the genes
  
  print(htr_i)
  
  htr_exp_i <- as.data.frame(i_dge %>%
                               subset(rownames(i_dge) %in%  htr_i)) # finding from df
  
  print(htr_exp_i)
  
  sample_num_i <- ncol(htr_exp_i) # for division
  
  # step 3: making a median of of the HTR tissues being present in the samples of each tissue
  
  median_htr_tissue_i <- data.frame(apply((as.matrix(htr_exp_i)), 1, median))
  
  names(median_htr_tissue_i)[names(median_htr_tissue_i) == 'apply..as.matrix.htr_exp_i....1..median.'] <-sprintf( "%s", tissue_list[i])
  
  htr_median_list[[i]] <- median_htr_tissue_i # adding to the empty list
  
  
  # step 4: making a maximum of of the HTR tissues being present in the samples of each tissue
  
  maximum_htr_tissue_i <- data.frame(apply((as.matrix(htr_exp_i)), 1, max))
  
  names(maximum_htr_tissue_i)[names(maximum_htr_tissue_i) == 'apply..as.matrix.htr_exp_i....1..max.'] <-sprintf( "%s", tissue_list[i])
  
  htr_max_list[[i]] <- maximum_htr_tissue_i # adding to the empty list
  
  # step 5: make a median of all matrix
  
  mat_median_i <- data.frame(median(as.matrix(log2_ab_dge_whole[, colnames(log2_ab_dge_whole) %in%  tissue_i])))
  
  names(mat_median_i)[names(mat_median_i) == 'median.as.matrix.log2_ab_dge_whole...colnames.log2_ab_dge_whole...in..'] <-sprintf( "%s", tissue_list[i])
  
  exp_median_list[[i]] <- mat_median_i # adding to the empty list
  
  # step 6: IQR
  
  iqr_htr_tissue_i <- data.frame(apply((as.matrix(htr_exp_i)), 1, iqr))
  
  names(iqr_htr_tissue_i)[names(maximum_htr_tissue_i) == 'apply..as.matrix.htr_exp_i....1..iqr.'] <-sprintf( "%s", tissue_list[i])
  
  htr_iqr_list[[i]] <- maximum_htr_tissue_i # adding to the empty list
}

htr_median_whole <- do.call(cbind, htr_median_list) # binding all values 

htr_maximum_whole <- do.call(cbind, htr_max_list) # binding all values 

htr_iqr_whole <- do.call(cbind, htr_iqr_list) # binding all values 

exp_median_whole <- do.call(cbind, exp_median_list)

htr2b_and_htr1f <- data.frame(rbind(exp_median_whole[1,],
                                    1.3,
                                    htr_median_whole[5,],
                                    htr_median_whole[7,],
                                    htr_maximum_whole[5,],
                                    htr_maximum_whole[7,],
                                    htr_iqr_whole[5,],
                                    htr_iqr_whole[7,]))


htr2b_and_htr1f <- round(htr2b_and_htr1f, digits = 3)

rownames(htr2b_and_htr1f) <- c("Tissue median", "Detection threshold", "HTR1F median", "HTR2B median", 
                               "HTR1F max", "HTR2B max", "HTR1F IQR", "HTR2B IQR")

table(htr2b_and_htr1f[2, ] > htr2b_and_htr1f[5, ])

table(htr2b_and_htr1f[2, ] > htr2b_and_htr1f[6, ])

write.csv(htr2b_and_htr1f, "Output/no_detection.csv")


# plot 

transplot <- data.frame(t(htr2b_and_htr1f))

names(transplot) <- c("Tissue median", "Detection threshold", "HTR1F median", "HTR2B median", 
                      "HTR1F max", "HTR2B max", "HTR1F IQR", "HTR2B IQR")

transplot$Tissues <- rownames(transplot)

min(transplot$`HTR2B IQR`)

transplot %>%
  ggplot() +
  geom_point(aes(x = Tissues, y = `Tissue median`), col = "black", alpha = 0.5, size = 3) +
  geom_line(aes(x = Tissues, y = `Tissue median`), col = "black", group = 1, size = 1) +
  geom_point(aes(x = Tissues, y = `Detection threshold`), col = "red", alpha = 0.5, size = 3) +
  geom_line(aes(x = Tissues, y = `Detection threshold`), col = "red", group = 1, size = 1) +
  geom_point(aes(x = Tissues, y = `HTR2B median`, size = `HTR2B IQR`), col = "#8da0cb", alpha = 0.5) +
  geom_line(aes(x = Tissues, y = `HTR2B median`), col = "#8da0cb", group = 1, size = 1.5) +
  geom_point(aes(x = Tissues, y = `HTR1F median`, size = `HTR1F IQR`), col = "#66c2a5", alpha = 0.5) +
  geom_line(aes(x = Tissues, y = `HTR1F median`), col = "#66c2a5", group = 1, size = 1.5) +
  scale_colour_manual(values = c("black", 
                                 "red", 
                                 "#8da0cb",
                                 "#66c2a5")) +
  scale_size_continuous(range = c(2, 8)) +
  labs(title = "Relative abundance of HTR1F and HT2B in brain tissues",
       x = "Tissue",
       y = "log2(cpm+1)",
       size = "IQR") +
  theme_tufte() +
  theme(text = element_text(family = "Arial"),
        axis.line = element_line(),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        # legend.position ="none",
        plot.title = element_text(hjust = 0.5))


#################################################################################

# PCA

#################################################################################

pca <- prcomp(t(as.matrix( htr_dge)))

tissues <- ab_tissues[, c(1,3)]

autoplot(pca, 
         x = 1,
         y = 2,
         data = tissues, 
         colour = "main_structure", 
         size = 4, 
         main = "PCA of htr") +
  scale_fill_manual(values = main_col) + 
  scale_color_manual(values =  main_col) +
  theme_classic() +
  theme(text = element_text(family = "Arial"))

main_str
