################################################################################

# packages 

################################################################################

# working with data 

library(dplyr) # for wrangling data frames 
library(tidyverse) # tidy data 


# statistics

library(matrixStats) # calculating matrix statistics 
library(rstatix) # for converting to factor 
library(stringr) # strings in a pattern 

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

# heatmap 

library(ComplexHeatmap)
library(GetoptLong)
library(circlize)
library(dendextend)
library(corrplot)
library(factoextra)
library(paletteer)

# pca

library(gplots) # plots 
library(ggfortify) # pca


# anova
library(rstatix)
library(reshape)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(plyr)
library(datarium)
library(paletteer)

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

# hierarchical clustering 

# using normlaised data 

################################################################################

# import data 

################################################################################
# counts 

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

# ab_ tissue

ab_samples_1 <- read.csv("Data/SampleAnnot.csv",
                         header = TRUE,
                         check.names = F,
                         row.names = 1)

ab_tissues_1 <- data.frame(ab_samples_1[ , c(6,7)])

ab_tissues_1$SAMPID <- rownames(ab_tissues_1)

ab_tissues_1$sub_structure <- as.factor(ab_tissues_1$sub_structure)



ab_samples_2 <- read.csv("Data2/SampleAnnot.csv",
                         header = TRUE,
                         check.names = F,
                         row.names = 1)

ab_tissues_2 <- data.frame(ab_samples_2[ , c(6,7)])

ab_tissues_2$SAMPID <- rownames(ab_tissues_2)

ab_tissues_2$sub_structure <- as.factor(ab_tissues_2$sub_structure)

ab_tissues <- rbind(ab_tissues_1, ab_tissues_2)

tissues <- ab_tissues[, c(1:3)]

class(tissues$main_structure) # character 

# htr data 

rownames(log2_ab_dge_whole)

htrgenes <- c("HTR1A",
              "HTR1B",
              "HTR1D",
              "HTR1E",
              "HTR1F",
              "HTR2A",
              "HTR2C",
              "HTR4",
              "HTR5A",
              "HTR6",
              "HTR7")

htr_ind <- rownames(log2_ab_dge_whole[htrgenes,]) # indexing the genes

htr_dge <- as.data.frame(log2_ab_dge_whole %>%
                           subset(rownames(log2_ab_dge_whole) %in%  htr_ind)) # finding from df

htr_dge

################################################################################

# plotting all rows

################################################################################

htr_dge

# data for hc 

mat_3 <- as.matrix(htr_dge)

# annotation for hc 

sample_annot <- ab_tissues[, c(1,3)] # sample id and main str

sample_annot <- sample_annot[,c(2,1)] # reorder 

# add colours by brain str 

sample_annot$colours <- 
ifelse(sample_annot$main_structure == "CgG" , "#9970ab",
  ifelse(sample_annot$main_structure == "Insula" , "#c51b7d",
    ifelse(sample_annot$main_structure == "FL", "#f46d43",
      ifelse(sample_annot$main_structure == "PL" , "#66c2a5" ,
        ifelse(sample_annot$main_structure == "TL" , "#5e4fa2",
          ifelse(sample_annot$main_structure == "PHG" , "#d1e5f0",
            ifelse(sample_annot$main_structure == "OL" , "#e6f598",
              ifelse(sample_annot$main_structure == "Str" , "#6baed6",
                  ifelse(sample_annot$main_structure ==  "GP" , "#fa9fb5",
                     "#9e0142")))))))))


# groups: cortex, cerebellum and basal ganglia

sample_annot$groups <- 
  ifelse(sample_annot$main_structure ==  "GP" , "BGl",
      ifelse(sample_annot$main_structure == "Str" , "BGl",
             ifelse(sample_annot$main_structure == "CbCx" ,"CbCx",
                    # ifelse(sample_annot$main_structure == "PHG" ,"PHG",
                    "Ctx")))

# add colours by groups: cortex, cerebellum and basal ganglia 

sample_annot$groups_col <- 
  ifelse(sample_annot$main_structure ==  "GP" , "#4daf4a",
         ifelse(sample_annot$main_structure == "Str" , "#4daf4a",
                ifelse(sample_annot$main_structure == "CbCx" ,"#9e0142",
                       # ifelse(sample_annot$main_structure == "PHG" , "#ffff99",
                       "#e41a1c")))
  
# donor info 

sample_annot$donor <- ifelse(str_detect(sample_annot$SAMPID, "S01"), "H0351.2001", "H0351.2002")

head(sample_annot)


# annotation for the heatmap 

hb <- HeatmapAnnotation('Brain structures' = sample_annot$main_structure, # data 1 
                        'Brain groups' = sample_annot$groups, # data 2
                        "Donor" = sample_annot$donor, # data 3
                        col = list('Brain structures' = c("CgG" = "#9970ab", # colour 1
                                                 "Ins" = "#c51b7d",
                                                 "FL" = "#f46d43",
                                                 "PL" = "#66c2a5" , #?
                                                 "TL" = "#5e4fa2",
                                                 "PHG" = "#d1e5f0", 
                                                 "OL" = "#e6f598",
                                                 "Str" = "#6baed6",
                                                 "GP" = "#fa9fb5",
                                                 "CbCx" = "#9e0142"),
                                   'Brain groups' = c("Ctx" = "salmon", # colour 2
                                                      "BGl" = "turquoise",
                                                      # "PHG" = "#ffff99",
                                                      "CbCx" = "#9e0142"),
                                   'Donor' = c("H0351.2001" = "#e41a1c",
                                               "H0351.2002" = "#4daf4a")))
                                  

# heatmap colour 

col_heatmap <- colorRamp2(c(0, 7.62), c("white", "#2c7bb6"))

# row dendogram 

row_dend_3 <- as.dendrogram(hclust(dist(mat_3)))


# heatmap 

# pdf("htr_brain_heatmap.pdf", width = 13, height = 9)


# heatmap on unscaled data 

htr <- Heatmap((mat_3),
        col = col_heatmap,
        clustering_distance_columns =  "pearson",
        clustering_method_rows = "single",
        clustering_method_columns  = "single",
        cluster_rows = color_branches(row_dend_3, k = 4),
        row_split = 2, 
        column_title = "Expression of the HTR family in Cortical and Subcortical Regions: Unscaled",
        name = "value",
        row_dend_reorder = T,
        column_names_gp =  gpar(fontsize = 9),
        show_column_names = F,
        top_annotation = hb)

draw(htr)

# dev.off()


# heatmap on scales data 


# mat 3 z score

mat3_zscore <- t(scale(t(mat_3)))

head(mat3_zscore)

min(mat3_zscore) #- 4.72

median(mat3_zscore) # 0.15

max(mat3_zscore) # 4.24

hist(mat3_zscore) # centers around 1

# changing color

col_scale <- colorRamp2(c(-5, 5), c("white", "#2c7bb6"))


htr <- Heatmap((mat_3),
               col = col_scale,
               clustering_distance_columns =  "pearson",
               clustering_method_rows = "single",
               clustering_method_columns  = "single",
               cluster_rows = color_branches(row_dend_3, k = 4),
               row_split = 2, 
               column_title = "Expression of the HTR family in Cortical and Subcortical Regions: Scaled",
               name = "value",
               row_dend_reorder = T,
               column_names_gp =  gpar(fontsize = 9),
               show_column_names = F,
               top_annotation = hb)

draw(htr)

#################################################################################

# test PCA

#################################################################################
plot_main <-  c("Ctx" = "#e41a1c",
                "BGl" = "#4daf4a",
                "PHG" = "#ffff99",
                "CbCx" = "#9e0142")

## Perform pca

pca <- prcomp(t(mat_3))


autoplot(pca, 
         x = 1,
         y = 3,
         # frame = TRUE,
         # frame.type = 'norm',
         shape = F,
         data = sample_annot, 
         colour = "groups", 
         size = 4, 
         main = "PCA of htr") +
  scale_fill_manual(values = plot_main) + 
  scale_color_manual(values =  plot_main) +
  theme_classic() +
  theme(text = element_text(family = "Arial"))

#################################################################################

# Heatmap data by median and % age 

#################################################################################

# matrix of htr

exp_median_list <- list()

exp_prob_list <- list()

tissue_list <- levels(as.factor(ab_tissues$sub_structure))

for (i in 1:length (tissue_list)) # remember to change to 1: length (tissue_list)
{
  # step 1 selecting the samples from each tissue 
  
  tissue_i <- ab_tissues %>%
    filter (ab_tissues$sub_structure == tissue_list[i]) %>%
    dplyr::select(SAMPID) %>% unlist()
  
  print(tissue_list[i])
  
  # step 2 divide 
  
  i_dge <- data.frame(log2_ab_dge_whole[, colnames(log2_ab_dge_whole) %in%  tissue_i])
  
  htr_i <- rownames(i_dge[htrgenes,]) # indexing the genes
  
  htr_exp_i <- as.data.frame(i_dge %>%
                               subset(rownames(i_dge) %in%  htr_i)) # finding from df
  
  sample_num_i <- ncol(htr_exp_i) # for division
  
  # step 3: making a percetange of of the HTR tissues being present in the samples of each tissue
  
  bin_htr_tissue_i <- data.frame(ifelse(htr_exp_i[,1:sample_num_i] > 1.3, 1, 0)) # binarising matrix
  
  htr_probability_i <- (data.frame(rowSums(bin_htr_tissue_i)/sample_num_i)*100) # number of 1s
  
  print(dim(htr_probability_i))
  
  names(htr_probability_i)[names(htr_probability_i) == 'rowSums.bin_htr_tissue_i..sample_num_i'] <-
    sprintf( "%s", tissue_list[i])
  
  exp_prob_list[[i]] <- htr_probability_i # adding to the empty list
  
  print(length(exp_prob_list))
  
  # step 4 making a median of of the HTR tissues being present in the samples of each tissue
  
  median_htr_tissue_i <- data.frame(apply((as.matrix(htr_exp_i)), 1, median))
  
  names(median_htr_tissue_i)[names(median_htr_tissue_i) == 'apply..as.matrix.htr_exp_i....1..median.'] <-sprintf( "%s", tissue_list[i])
  
  exp_median_list[[i]] <- median_htr_tissue_i # adding to the empty list
  
}


htr_median_whole <- do.call(cbind, exp_median_list) # binding all values 

htr_median_whole$Receptor <- rownames(htr_median_whole) # renaming by receptor 


htr_prob_whole <- do.call(cbind, exp_prob_list) # binding all prob values 

htr_prob_whole$Receptor <- rownames(htr_prob_whole) # renaming by receptor 


mat_1 <- as.matrix(htr_median_whole[, 1:29])

mat_2 <- scale(mat_1)

mat_3 <- t(scale(t(mat_1)))

mat_prob <- as.matrix(htr_prob_whole[, 1:29])

cor((mat_1)) # no sd = 0 

cor(mat_2)

hclust(as.dist(1-cor((mat_1)))) # no error 

hclust(as.dist(1-cor((mat_2)))) # no error 

################################################################################

# K means setting 

################################################################################

# elbow 
#  choose a number of clusters so that adding another cluster doesn't improve much better the total WSS.
#  WSS measures the compactness of the clustering and we want it to be as small as possible.

fviz_nbclust(mat_1, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method") # 4

# Silhouette method

# The optimal number of clusters k is the one that maximize the average silhouette over 
# a range of possible values for k (Kaufman and Rousseeuw 1990).                

fviz_nbclust(mat_1, kmeans, method = "silhouette") +
  geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Silhouette method") # 2 


# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
# The estimate of the optimal clusters will be a value that 
# maximize the gap statistic (i.e, that yields the largest gap statistic).         

set.seed(123)

fviz_nbclust(mat_1, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method") #  9

head(mat_2)

# kmeans clustering 

k2 <- kmeans(hclust(dist(mat_1)), 2, nstart = 25)
k3 <- kmeans(hclust(dist(mat_1)), 3, nstart = 25)
k4 <- kmeans(hclust(dist(mat_1)), 4, nstart = 25)
k5 <- kmeans(hclust(dist(mat_1)), 5, nstart = 25)
k6 <- kmeans(hclust(dist(mat_1)), 6, nstart = 25)
k7 <- kmeans(hclust(dist(mat_1)), 7, nstart = 25)

# visualizing 

p2 <- fviz_cluster(k2, data = hclust(dist(mat_1))) + ggtitle("k = 2") + theme_tufte()
p3 <- fviz_cluster(k3, data = hclust(dist(mat_1))) + ggtitle("k = 3") + theme_tufte()
p4 <- fviz_cluster(k4, data = hclust(dist(mat_1))) + ggtitle("k = 4") + theme_tufte()
p5 <- fviz_cluster(k5, data = hclust(dist(mat_1))) + ggtitle("k = 5") + theme_tufte()
p6 <- fviz_cluster(k6, data = hclust(dist(mat_1))) + ggtitle("k = 6") + theme_tufte()
p7 <- fviz_cluster(k7, data = hclust(dist(mat_1))) + ggtitle("k = 7") + theme_tufte()

# comparing different Ks

gridExtra::grid.arrange(p2, p3, p4, p5, nrow = 2)


################################################################################

# preparing heatmap annotations

################################################################################

# heatmap colour

col_fun <- colorRamp2(c(0, 6), c("white", "#0e7f92"))

col_per <- colorRamp2(c(0, 100), c("white", "#0e7f92"))

# label colour 

tissue_colours <- read.csv("allencol.txt",
                           header = F) # contains all ab colours 

colours <- list('tissues' = c("AnG-i" = "#66c2a5",
                              "AnG-s" = "#66c2a5",
                              "Caudate" = "#3288bd",
                              "CbCx" = "#9e0142",
                              "CgG" = "#9970ab",
                              "FuG-its" = "#5e4fa2",
                              "GP" = "#fdae61",
                              "GRe" = "#f46d43",
                              "Insula" = "#c51b7d",
                              "ITG" = "#5e4fa2",
                              "MFG" = "#f46d43",
                              "MTG" = "#5e4fa2",
                              "OrbGyri" = "#f46d43",
                              "orIFG" = "#f46d43",
                              "PCLa-i" = "#f46d43",
                              "PCLa-s" = "#f46d43",
                              "Pcu" = "#66c2a5",
                              "pest_V2" = "#e6f598",
                              "PHG" = "#d1e5f0",
                              "PoG-cs" = "#66c2a5",
                              "PoG-l" = "#66c2a5",
                              "PrG" = "#f46d43",
                              "Putamen" = "#3288bd",
                              "SFG-l" = "#f46d43",
                              "SFG-m" = "#f46d43",
                              "SMG-i" = "#66c2a5",
                              "SPL" = "#66c2a5",
                              "STG" = "#5e4fa2",
                              "str_V1" = "#e6f598"))

ha <- HeatmapAnnotation("tissues" = tissue_colours$V2, col = colours)


################################################################################

# spearman complete heatmaps of median and % exp value in tissues 

################################################################################

# median 

row_dend_1 <- as.dendrogram(hclust(dist(mat_1)))

Heatmap((mat_1),
        col = col_fun,
        clustering_distance_columns = "spearman",
        clustering_method_rows = "complete",
        clustering_method_columns  = "complete",
        cluster_rows = color_branches(row_dend_1, k = 4),
        column_title = "Median Expression Value", #Grouping HTR Suptypes by Distribution in 232 Allen Brain Sample unscaled
        name = "value",
        row_dend_reorder = T,
        column_names_gp =  gpar(fontsize = 9),
        top_annotation = ha) 

# prob 

row_dend_2 <- as.dendrogram(hclust(dist(mat_prob)))

Heatmap((mat_prob),
        col = col_per,
        clustering_distance_columns = "spearman",
        clustering_method_rows = "complete",
        clustering_method_columns  = "complete",
        cluster_rows = color_branches(row_dend_1, k = 4),
        column_title = "Probability of Detection", #Grouping HTR Suptypes by Distribution in 232 Allen Brain Sample
        name = "value",
        row_dend_reorder = T,
        column_names_gp =  gpar(fontsize = 9),
        top_annotation = ha) 

################################################################################

# correlation between median and %age 

################################################################################

# corrplot of tissues based on htr

mat_medians <- mat_1

colnames(mat_medians) <- paste("median", colnames(mat_medians), sep = "_")

mat_percentage <- mat_prob

colnames(mat_percentage) <- paste("percentage", colnames(mat_percentage), sep = "_")

cor_two_mat <-  cor(mat_medians, mat_percentage, method = "pearson")

cor_two_mat

corrplot(cor_two_mat, 
         method = 'color',  # square
         order = 'hclust', # sorted by clustering
         # p.mat = my_Pval,
         # sig.level = c(0.001, 0.01), # levels of asterisk 
         # pch.cex = 0.9, # size of the p vales 
         # insig = 'label_sig', # label
         tl.col = 'black', # labels are black
         tl.cex = 0.8, # size of the labels
         cl.pos = 'b', # annotation bar at bottom
         cl.cex = 0.8, # arranging the text size to fit
         title = "Correlation between detection of expression amd median of expression",
         mar=c(0,0,1,0),
         col = brewer.pal(n = 10, name = 'RdYlBu'))  # colour scheme


# significance for each cell 

head(mat_medians)

mat_comb <- cbind(mat_medians, mat_percentage)

P <- cor.mtest(mat_comb, conf.level = 0.95) # sig 

head(P)

Pval <- P$p

my_Pval <- Pval[c(1:29), c(30:58)]

dim(my_Pval)

################################################################################

# statistics 

################################################################################

# median

d_1 <- as.data.frame(mat_1)

d_1$Genes <- rownames(d_1)

d_1 <- d_1 %>%
  gather(key = "sample", value = "Value", -"Genes")

d_1 <- d_1[ , c(3)] 

# prob 

d_2 <- as.data.frame(mat_prob)

d_2$Genes <- rownames(d_2)

d_2 <- d_2 %>%
  gather(key = "sample", value = "Value", -"Genes")

d_2 <- d_2[ , c(3)] 

length(d_1)

length(d_2)

# pearson 

cor_d_p <- cor(d_1, d_2, method = "pearson") 

cor_d_p # 0.70


# t test 

x <- htr_median_whole[1, 1:29]

y <- htr_median_whole[4, 1:29]

tt <- t.test(x,y)

tt$p.value

