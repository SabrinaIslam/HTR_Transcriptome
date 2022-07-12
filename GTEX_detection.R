
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
library(viridis) # colour 

# heatmap 

library(ComplexHeatmap) # heatmap 
library(GetoptLong) # for images 
library(circlize) # for coloring 
library(dendextend) # dendogram 
library(corrplot)
library(factoextra)



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

# making heatmaps using the percentage of reproducible detection calculated before

################################################################################

# reading the data 

prob_htr <- read.csv("Data/HTR_Prob.csv",
                 header = TRUE, 
                 check.names = FALSE,
                 row.names = 1)

# the data is the percentage of the 12 receptors being detected in 
# any of the 54 tissues 

######################################################################

# heatmap specifications

####################################################################

# we need a dendogram from the row, a colour-code for the tissue annotation 
# and a color-code for the heatmap itself 

# step 1: annotation bar by tissues and setting colour   

# matrix of gtex colours 

tissue_colours <- read.csv("Data/gtexcol.csv",
                           header = F) # gtex tissue colours matrix 

colnames(prob_htr) %in% tissue_colours$V1 # all tissue names are exact same 

# extracting just tissue names

tissues <- tissue_colours$V1   

ann <- data.frame(tissues) # create tissue color data frame 

colnames(ann)

ann # this will be the annotation data 

# mappping a colours to each tissue 

class(ann$tissues) # character 

# creating annoation with ann as df and colours as list 

colAnn <- HeatmapAnnotation("Tissue" = ann$tissues,
                            which = 'col',
                            col = list('Tissue' = c(
                            "Adipose - Subcutaneous" = "#ef883a",                   
                            "Adipose - Visceral (Omentum)" = "#ffa804",             
                            "Adrenal Gland" = "#34de34",                             
                            "Artery - Aorta" = "#ff5551",                           
                            "Artery - Coronary" = "#ffa599",                         
                            "Artery - Tibial" = "#fd0000",                          
                            "Bladder" = "#814540",                                   
                            "Brain - Amygdala" = "#eee03d",                         
                            "Brain - Anterior cingulate cortex (BA24)" = "#eee03d",  
                            "Brain - Caudate (basal ganglia)" = "#eee03d",          
                            "Brain - Cerebellar Hemisphere" = "#eee03d",            
                            "Brain - Cerebellum" = "#eee03d",                      
                            "Brain - Cortex" = "#eee03d",                            
                            "Brain - Frontal Cortex (BA9)" = "#eee03d",             
                            "Brain - Hippocampus" = "#eee03d",                       
                            "Brain - Hypothalamus" = "#eee03d",                     
                            "Brain - Nucleus accumbens (basal ganglia)" = "#eee03d", 
                            "Brain - Putamen (basal ganglia)" = "#eee03d",          
                            "Brain - Spinal cord (cervical c-1)" = "#eee03d",        
                            "Brain - Substantia nigra" = "#eee03d",                 
                            "Breast - Mammary Tissue" = "#30ccd1",                   
                            "Cells - Cultured fibroblasts" = "#abedff",             
                            "Cells - EBV-transformed lymphocytes" = "#e2baf4",       
                            "Cervix - Ectocervix" = "#ffe3e6",                      
                            "Cervix - Endocervix" = "#e5d4ed",                      
                            "Colon - Sigmoid" =  "#f8dbba",                         
                            "Colon - Transverse" = "#e5cdaa",                        
                            "Esophagus - Gastroesophageal Junction" = "#8d7450",    
                            "Esophagus - Mucosa" = "#522302",                        
                            "Esophagus - Muscularis"  = "#bc9885",                  
                            "Fallopian Tube" =  "#ffd3ff",                           
                            "Heart - Atrial Appendage" = "#9c00fc",                 
                            "Heart - Left Ventricle" = "#6a0091",                  
                            "Kidney - Cortex" = "#89ffee",                         
                            "Kidney - Medulla" = "#99fde0",                         
                            "Liver" = "#bac29d",                                    
                            "Lung" = "#c9fc80",                                     
                            "Minor Salivary Gland" = "#99BB88",                     
                            "Muscle - Skeletal" = "#d5d7f6",                         
                            "Nerve - Tibial" = "#ffd600",                          
                            "Ovary" = "#fce3fc",                                     
                            "Pancreas" = "#9b5522",                                 
                            "Pituitary" = "#abff93",                                 
                            "Prostate"= "#e1dcdc",                                 
                            "Skin - Not Sun Exposed (Suprapubic)" = "#0000fa",       
                            "Skin - Sun Exposed (Lower leg)" = "#7977fe",           
                            "Small Intestine - Terminal Ileum" ="#525026" ,          
                            "Spleen" = "#7a8750",                                   
                            "Stomach" = "#fdde8c",                                   
                            "Testis" = "#aaaaaa",                                  
                            "Thyroid" = "#006701",                                  
                            "Uterus" = "#f96cff",                                   
                            "Vagina" = "#ee6697",                                    
                            "Whole Blood" = "#ff00bb")), # this will be the annotation list ,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))


# step 2: organising the dendogram 

row_dend <- as.dendrogram(hclust(dist(as.matrix(prob_htr[, 1:54]))))

# the dendogram ia mased on the hierarchical clustering of the data 

# step 3: color for heatmap 

col_fun <- colorRamp2(c(0, 100), c("white", "#0e7f92")) 

##############################################################################

# K means setting 

##############################################################################

# we want to specify how many clusters we want to see according to the dendogram 

M1 <- (as.matrix(prob_htr[, 1:54])) # matrix 

M2 <- t(as.matrix(prob_htr[, 1:54])) # matrix 

# kmeans clustering 

k3 <- kmeans(hclust(dist(M1)), 3, nstart = 25)
k4 <- kmeans(hclust(dist(M1)), 4, nstart = 25)
k5 <- kmeans(hclust(dist(M1)), 5, nstart = 25)
k6 <- kmeans(hclust(dist(M1)), 6, nstart = 25)
k7 <- kmeans(hclust(dist(M1)), 7, nstart = 25)

# elbow 
#  choose a number of clusters so that adding another cluster doesn't improve much better the total WSS.
#  WSS measures the compactness of the clustering and we want it to be as small as possible.
        fviz_nbclust(M1, kmeans, method = "wss") +
        geom_vline(xintercept = 4, linetype = 2)+
        labs(subtitle = "Elbow method")


# optimal k is 4
        
# Silhouette method

# The optimal number of clusters k is the one that maximize the average silhouette over 
# a range of possible values for k (Kaufman and Rousseeuw 1990).                
  
        fviz_nbclust(M1, kmeans, method = "silhouette")+
        labs(subtitle = "Silhouette method")

# optimal ks are 2 and 4
        
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
# The estimate of the optimal clusters will be a value that 
# maximize the gap statistic (i.e, that yields the largest gap statistic).         
        
        set.seed(123)
        
        fviz_nbclust(M1, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
        labs(subtitle = "Gap statistic method")

# optimal k is 4
        
# visulaising 

p3 <- fviz_cluster(k3, data = hclust(dist(M1))) + ggtitle("k = 3") + theme_clean()
p4 <- fviz_cluster(k4, data = hclust(dist(M1))) + ggtitle("k = 4") + theme_clean()
p5 <- fviz_cluster(k5, data = hclust(dist(M1))) + ggtitle("k = 5") + theme_clean()
p6 <- fviz_cluster(k6, data = hclust(dist(M1))) + ggtitle("k = 6") + theme_clean()
p7 <- fviz_cluster(k7, data = hclust(dist(M1))) + ggtitle("k = 4") + theme_clean()

p5 <- fviz_cluster(k5, data = hclust(dist(M1))) + 
  scale_color_brewer('Cluster', palette='Set2') + 
  scale_fill_brewer('Cluster', palette='Set2') +
  scale_shape_manual('Cluster', values=c(16,16,16,16,16)) +
  ggtitle("HTR family clusters") + theme_clean()

p5

# comparing 

gridExtra::grid.arrange(p3, p4, p5, p6, p7, nrow = 2)

# 4 k seems reasonable 

####################################################

# quick heatmaps

####################################################

heatmap(as.matrix(prob_htr[, 1:54]))

##################################################

# proper heatmaps

# we will test out all linkage methods: single, average, and 

# which we will apply on both row and column 

# as well as both distance calculation (pearson and spearman)


# Pearson, single 

Heatmap(as.matrix(prob_htr[, 1:54]), 
        col = col_fun,
        clustering_distance_columns = "pearson",
        clustering_method_rows = "single",
        clustering_method_columns  = "single",
        cluster_rows = color_branches(row_dend, k = 4),
        name = "value",
        column_title = "Pearson Single",
        column_names_gp =  gpar(fontsize = 9),
        top_annotation = colAnn)

# HTR2B, (HTR5A, HTR1E), (HTR2A, HTR1B, HTR7), (HTR2C, HTR1D, HTR6,) HTR1A, (HTR1F, HTR4)

# Pearson, average 

Heatmap(as.matrix(prob_htr[, 1:54]), 
        col = col_fun,
        clustering_distance_columns = "pearson",
        clustering_method_rows = "average",
        clustering_method_columns  = "average",
        name = "value",
        cluster_rows = color_branches(row_dend, k = 4),
        column_title = "Pearson Average",
        column_names_gp =  gpar(fontsize = 9),
        top_annotation = colAnn)

# HTR2B, (HTR5A, HTR1E), (HTR2A, HTR1B, HTR7), (HTR2C, HTR1D, HTR6,) HTR1A, (HTR1F, HTR4)

# Pearson, complete

Heatmap(as.matrix(prob_htr[, 1:54]), 
        col = col_fun,
        clustering_distance_columns = "pearson",
        clustering_method_rows = "complete",
        clustering_method_columns  = "complete",
        cluster_rows = color_branches(row_dend, k = 4),
        name = "value",
        column_title = "Pearson Complete",
        column_names_gp =  gpar(fontsize = 9),
        top_annotation = colAnn)

# HTR2B, (HTR5A, HTR1E), (HTR2A, HTR1B, HTR7), (HTR2C, HTR1D, HTR6,) HTR1A, (HTR1F, HTR4)


##########################################

# Spearman, single 

Heatmap(as.matrix(prob_htr[, 1:54]), 
        col = col_fun,
        clustering_distance_columns =  "spearman",
        clustering_method_rows = "single",
        clustering_method_columns  = "single",
        name = "value",
        column_title = "Spearman Single",
        cluster_rows = color_branches(row_dend, k = 4),
        column_names_gp =  gpar(fontsize = 9),
        top_annotation = colAnn)

# HTR2B, (HTR5A, HTR1E), (HTR2A, HTR1B, HTR7), (HTR2C, HTR1D, HTR6,) HTR1A, (HTR1F, HTR4)

# Spearman, average 

Heatmap(as.matrix(prob_htr[, 1:54]), 
        col = col_fun,
        clustering_distance_columns =  "spearman",
        clustering_method_rows = "average",
        clustering_method_columns  = "average",
        cluster_rows = color_branches(row_dend, k = 4),
        name = "value",
        column_title = "Spearman Average",
        column_names_gp =  gpar(fontsize = 9),
        top_annotation = colAnn)

# HTR2B, (HTR5A, HTR1E), (HTR2A, HTR1B, HTR7), (HTR2C, HTR1D, HTR6,) HTR1A, (HTR1F, HTR4)

# Spearman complete

pdf("prob_heatmap.pdf", width = 13, height = 9)

ht <- Heatmap(as.matrix(prob_htr[, 1:54]), 
        col = col_fun,
        clustering_distance_columns =  "spearman",
        clustering_method_rows = "complete",
        clustering_method_columns  = "complete",
        cluster_rows = color_branches(row_dend, k = 4),
        name = "value",
        column_title = "HTR Distribution in Human Tissues",
        column_names_gp =  gpar(fontsize = 9),
        top_annotation = colAnn)
draw(ht,
     heatmap_legend_side = "right", 
     annotation_legend_side = "left")

dev.off()

################################################################################

# PPT 

################################################################################

colAnn2 <- HeatmapAnnotation("Tissue" = ann$tissues,
                             "BrainVsNonBrain" = ann$tissues,
                             which = 'col',
                             col = list('Tissue' = c("Adipose - Subcutaneous" = "#ef883a",                   
                                                     "Adipose - Visceral (Omentum)" = "#ef883a",             
                                                     "Adrenal Gland" = "#34de34",                             
                                                     "Artery - Aorta" = "#fd0000",                           
                                                     "Artery - Coronary" = "#fd0000",                         
                                                     "Artery - Tibial" = "#fd0000",                          
                                                     "Bladder" = "#814540",                                   
                                                     "Brain - Amygdala" = "#eee03d",                         
                                                     "Brain - Anterior cingulate cortex (BA24)" = "#eee03d",  
                                                     "Brain - Caudate (basal ganglia)" = "#eee03d",          
                                                     "Brain - Cerebellar Hemisphere" = "#eee03d",            
                                                     "Brain - Cerebellum" = "#eee03d",                      
                                                     "Brain - Cortex" = "#eee03d",                            
                                                     "Brain - Frontal Cortex (BA9)" = "#eee03d",             
                                                     "Brain - Hippocampus" = "#eee03d",                       
                                                     "Brain - Hypothalamus" = "#eee03d",                     
                                                     "Brain - Nucleus accumbens (basal ganglia)" = "#eee03d", 
                                                     "Brain - Putamen (basal ganglia)" = "#eee03d",          
                                                     "Brain - Spinal cord (cervical c-1)" = "#eee03d",        
                                                     "Brain - Substantia nigra" = "#eee03d",                 
                                                     "Breast - Mammary Tissue" = "#30ccd1",                   
                                                     "Cells - Cultured fibroblasts" = "#abedff",             
                                                     "Cells - EBV-transformed lymphocytes" = "#e2baf4",       
                                                     "Cervix - Ectocervix" = "#e5d4ed",                      
                                                     "Cervix - Endocervix" = "#e5d4ed",                      
                                                     "Colon - Sigmoid" =  "#f8dbba",                         
                                                     "Colon - Transverse" = "#f8dbba",                        
                                                     "Esophagus - Gastroesophageal Junction" = "#8d7450",    
                                                     "Esophagus - Mucosa" = "#8d7450",                        
                                                     "Esophagus - Muscularis"  = "#8d7450",                  
                                                     "Fallopian Tube" =  "#ffd3ff",                           
                                                     "Heart - Atrial Appendage" = "#9c00fc",                 
                                                     "Heart - Left Ventricle" = "#9c00fc",                  
                                                     "Kidney - Cortex" = "#99fde0",                         
                                                     "Kidney - Medulla" = "#99fde0",                         
                                                     "Liver" = "#c9fc80",                                    
                                                     "Lung" = "#ceddc2",                                     
                                                     "Minor Salivary Gland" = "#99BB88",                     
                                                     "Muscle - Skeletal" = "#d5d7f6",                         
                                                     "Nerve - Tibial" = "#ffd600",                          
                                                     "Ovary" = "#ffd1ff",                                     
                                                     "Pancreas" = "#9b5522",                                 
                                                     "Pituitary" = "#abff93",                                 
                                                     "Prostate"= "#e1dcdc",                                 
                                                     "Skin - Not Sun Exposed (Suprapubic)" = "#7977fe",       
                                                     "Skin - Sun Exposed (Lower leg)" = "#7977fe",           
                                                     "Small Intestine - Terminal Ileum" ="#f8dbba" ,          
                                                     "Spleen" = "#7a8750",                                   
                                                     "Stomach" = "#fdde8c",                                   
                                                     "Testis" = "#aaaaaa",                                  
                                                     "Thyroid" = "#006701",                                  
                                                     "Uterus" = "#f96cff",                                   
                                                     "Vagina" = "#ee6697",                                    
                                                     "Whole Blood" = "#ff00bb"), # this will be the annotation list 
                                        'BrainVsNonBrain' = c("Adipose - Subcutaneous" = "#fd0000",                   
                                                              "Adipose - Visceral (Omentum)" = "#fd0000",             
                                                              "Adrenal Gland" = "#fd0000",                             
                                                              "Artery - Aorta" = "#fd0000",                           
                                                              "Artery - Coronary" = "#fd0000",                         
                                                              "Artery - Tibial" = "#fd0000",                          
                                                              "Bladder" = "#fd0000",                                   
                                                              "Brain - Amygdala" = "#eee03d",                         
                                                              "Brain - Anterior cingulate cortex (BA24)" = "#eee03d",  
                                                              "Brain - Caudate (basal ganglia)" = "#eee03d",          
                                                              "Brain - Cerebellar Hemisphere" = "#eee03d",            
                                                              "Brain - Cerebellum" = "#eee03d",                      
                                                              "Brain - Cortex" = "#eee03d",                            
                                                              "Brain - Frontal Cortex (BA9)" = "#eee03d",             
                                                              "Brain - Hippocampus" = "#eee03d",                       
                                                              "Brain - Hypothalamus" = "#eee03d",                     
                                                              "Brain - Nucleus accumbens (basal ganglia)" = "#eee03d", 
                                                              "Brain - Putamen (basal ganglia)" = "#eee03d",          
                                                              "Brain - Spinal cord (cervical c-1)" = "#eee03d",        
                                                              "Brain - Substantia nigra" = "#eee03d",                 
                                                              "Breast - Mammary Tissue" = "#fd0000",                   
                                                              "Cells - Cultured fibroblasts" = "#fd0000",             
                                                              "Cells - EBV-transformed lymphocytes" = "#fd0000",       
                                                              "Cervix - Ectocervix" = "#fd0000",                      
                                                              "Cervix - Endocervix" = "#fd0000",                      
                                                              "Colon - Sigmoid" =  "#fd0000",                         
                                                              "Colon - Transverse" = "#fd0000",                        
                                                              "Esophagus - Gastroesophageal Junction" = "#fd0000",    
                                                              "Esophagus - Mucosa" = "#fd0000",                        
                                                              "Esophagus - Muscularis"  = "#fd0000",                  
                                                              "Fallopian Tube" =  "#fd0000",                           
                                                              "Heart - Atrial Appendage" = "#fd0000",                 
                                                              "Heart - Left Ventricle" = "#fd0000",                  
                                                              "Kidney - Cortex" = "#fd0000",                         
                                                              "Kidney - Medulla" = "#fd0000",                         
                                                              "Liver" = "#fd0000",                                    
                                                              "Lung" = "#fd0000",                                     
                                                              "Minor Salivary Gland" = "#fd0000",                     
                                                              "Muscle - Skeletal" = "#fd0000",                         
                                                              "Nerve - Tibial" = "#fd0000",                          
                                                              "Ovary" = "#fd0000",                                     
                                                              "Pancreas" = "#fd0000",                                 
                                                              "Pituitary" = "#fd0000",                                 
                                                              "Prostate"= "#fd0000",                                 
                                                              "Skin - Not Sun Exposed (Suprapubic)" = "#fd0000",       
                                                              "Skin - Sun Exposed (Lower leg)" = "#fd0000",           
                                                              "Small Intestine - Terminal Ileum" ="#fd0000" ,          
                                                              "Spleen" = "#fd0000",                                   
                                                              "Stomach" = "#fd0000",                                   
                                                              "Testis" = "#fd0000",                                  
                                                              "Thyroid" = "#fd0000",                                  
                                                              "Uterus" = "#fd0000",                                   
                                                              "Vagina" = "#fd0000",                                    
                                                              "Whole Blood" = "#fd0000")),
                             annotation_width = unit(c(1, 4), 'cm'),
                             gap = unit(1, 'mm'))



pdf(qq("poster.pdf"), width = 13, height = 9)

ht <- Heatmap(as.matrix(prob_htr[, 1:54]), 
              col = col_fun,
              clustering_distance_columns =  "spearman",
              clustering_method_rows = "complete",
              clustering_method_columns  = "complete",
              cluster_rows = color_branches(row_dend, k = 4),
              name = "value",
              column_title = "HTR Distribution in Human Tissues",
              show_column_names = F,
              column_names_gp =  gpar(fontsize = 9),
              top_annotation = colAnn2)


draw(ht,
     show_heatmap_legend = FALSE, 
     show_annotation_legend = FALSE)

dev.off()

# HTR2B, (HTR5A, HTR1E), (HTR2A, HTR1B, HTR7), (HTR2C, HTR1D, HTR6,) HTR1A, (HTR1F, HTR4)

#################################################################################

# subsetting 

################################################################################

# brain

ann_brain <- ann[c(8:20), ]

# colurs list 

colours_brain <- list('Tissues' = c("Brain - Amygdala" = "#eee03d",                         
                                  "Brain - Anterior cingulate cortex (BA24)" = "#eee03d",  
                                  "Brain - Caudate (basal ganglia)" = "#eee03d",          
                                  "Brain - Cerebellar Hemisphere" = "#eee03d",            
                                  "Brain - Cerebellum" = "#eee03d",                      
                                  "Brain - Cortex" = "#eee03d",                            
                                  "Brain - Frontal Cortex (BA9)" = "#eee03d",             
                                  "Brain - Hippocampus" = "#eee03d",                       
                                  "Brain - Hypothalamus" = "#eee03d",                     
                                  "Brain - Nucleus accumbens (basal ganglia)" = "#eee03d", 
                                  "Brain - Putamen (basal ganglia)" = "#eee03d",          
                                  "Brain - Spinal cord (cervical c-1)" = "#eee03d",        
                                  "Brain - Substantia nigra" = "#eee03d"))

# creating annoation with ann df and colours list 

colAnn_brain <- HeatmapAnnotation(Tissues = ann_brain,
                                which = 'col',
                                col = colours_brain,
                                annotation_width = unit(c(1, 4), 'cm'),
                                gap = unit(1, 'mm'))

# dendogram 

row_dend_brain <- as.dendrogram(hclust(dist(as.matrix(prob_htr[, c(8:20)]))))

# heatmap 

pdf(qq("brain_heatmap.pdf"), width = 9, height = 9)

ht_brain <- Heatmap(as.matrix(prob_htr[, c(8:20)]), 
               col = col_fun,
               clustering_distance_columns =  "spearman",
               clustering_method_rows = "complete",
               clustering_method_columns  = "complete",
               cluster_rows = color_branches(row_dend_brain, k = 4),
               name = "value",
               column_title = "HTR distribution in the brain",
               show_column_names = T,
               column_names_gp =  gpar(fontsize = 9),
               top_annotation = colAnn_brain)

draw(ht_brain)

dev.off()


# subsetting the reproductive

# ann

ann_repro <- ann[c(24, 25, 31, 41, 44, 50, 52, 53), ]

# colurs list 

colours_repro <- list('Tissues' = c("Cervix - Ectocervix" = "#ffe3e6",                      
                                    "Cervix - Endocervix" = "#e5d4ed",                      
                                    "Fallopian Tube" =  "#ffd3ff",
                                    "Ovary" = "#ffd1ff",
                                    "Prostate"= "#e1dcdc",
                                    "Testis" = "#aaaaaa",
                                    "Uterus" = "#f96cff",                                   
                                    "Vagina" = "#ee6697"))

# creating annoation with ann df and colours list 

colAnn_repro <- HeatmapAnnotation(Tissues = ann_repro,
                                which = 'col',
                                col = colours_repro,
                                annotation_width = unit(c(1, 4), 'cm'),
                                gap = unit(1, 'mm'))

# dendogram 

row_dend_repro <- as.dendrogram(hclust(dist(as.matrix(prob_htr[, c(24, 25, 31, 41, 44, 50, 52, 53)]))))

# heatmap 

pdf(qq("repro_heatmap.pdf"), width = 9, height = 9)

ht_repro <- Heatmap(as.matrix(prob_htr[, c(24, 25, 31, 41, 44, 50, 52, 53)]), 
                  col = col_fun,
                  clustering_distance_columns =  "spearman",
                  clustering_method_rows = "complete",
                  clustering_method_columns  = "complete",
                  cluster_rows = color_branches(row_dend_repro, k = 4),
                  name = "value",
                  column_title = "HTR distribution in reproductive tissues",
                  show_column_names = T,
                  column_names_gp =  gpar(fontsize = 9),
                  top_annotation = colAnn_repro)

draw(ht_repro)

dev.off()


# subsetting ens

# ann

ann_ens <- ann[c(26, 27, 28, 29, 30, 36, 38, 42, 47, 48, 49), ]

# colurs list 

colours_ens <- list('Tissues' = c("Colon - Sigmoid" =  "#f8dbba",                         
                                  "Colon - Transverse" = "#e5cdaa",                        
                                  "Esophagus - Gastroesophageal Junction" = "#8d7450",    
                                  "Esophagus - Mucosa" = "#522302",                        
                                  "Esophagus - Muscularis"  = "#bc9885",
                                  "Liver" = "#c9fc80",
                                  "Minor Salivary Gland" = "#d5d7f6",
                                  "Pancreas" = "#9b5522",
                                  "Small Intestine - Terminal Ileum" ="#525026" ,          
                                  "Spleen" = "#7a8750",                                   
                                  "Stomach" = "#fdde8c"))

# creating annoation with ann df and colours list 

colAnn_ens <- HeatmapAnnotation(Tissues = ann_ens,
                                  which = 'col',
                                  col = colours_ens,
                                  annotation_width = unit(c(1, 4), 'cm'),
                                  gap = unit(1, 'mm'))

# dendogram 

row_dend_ens <- as.dendrogram(hclust(dist(as.matrix(prob_htr[, c(26, 27, 28, 29, 30, 36, 38, 42, 47, 48, 49)]))))

# heatmap 

pdf(qq("ens_heatmap.pdf"), width = 9, height = 9)

ht_ens <- Heatmap(as.matrix(prob_htr[, c(26, 27, 28, 29, 30, 36, 38, 42, 47, 48, 49)]), 
                    col = col_fun,
                    clustering_distance_columns =  "spearman",
                    clustering_method_rows = "complete",
                    clustering_method_columns  = "complete",
                    cluster_rows = color_branches(row_dend_ens, k = 4),
                    name = "value",
                    column_title = "HTR distribution in the digestive tract",
                    show_column_names = T,
                    column_names_gp =  gpar(fontsize = 9),
                    top_annotation = colAnn_ens)

draw(ht_ens)

dev.off()


# subsetting cvd 

# ann

ann_artery <- ann[c(4,5,6, 32, 33), ]

# colurs list 

colours_artery <- list('Tissues' = c("Artery - Aorta" = "#ff5551",                           
                                  "Artery - Coronary" = "#ffa599",                         
                                  "Artery - Tibial" = "#fd0000", 
                                  "Heart - Atrial Appendage" = "#9c00fc",                 
                                  "Heart - Left Ventricle" = "#6a0091"))

# creating annoation with ann df and colours list 

colAnn_artery <- HeatmapAnnotation(Tissues = ann_artery,
                                which = 'col',
                                col = colours_artery,
                                annotation_width = unit(c(1, 4), 'cm'),
                                gap = unit(1, 'mm'))

# dendogram 

row_dend_artery <- as.dendrogram(hclust(dist(as.matrix(prob_htr[, c(4,5,6, 32, 33)]))))

# heatmap 

pdf(qq("artery_heatmap.pdf"), width = 9, height = 9)

ht_artery <- Heatmap(as.matrix(prob_htr[, c(4,5,6, 32, 33)]), 
                  col = col_fun,
                  clustering_distance_columns =  "spearman",
                  clustering_method_rows = "complete",
                  clustering_method_columns  = "complete",
                  cluster_rows = color_branches(row_dend_artery, k = 4),
                  name = "value",
                  column_title = "HTR distribution in the heart and the arteries",
                  show_column_names = T,
                  column_names_gp =  gpar(fontsize = 9),
                  top_annotation = colAnn_artery)

draw(ht_artery)

dev.off()


# fatty subsetting 

# ann

ann_fat <- ann[c(1,2,21), ]

# colurs list 

colours_fat <- list('Tissues' = c("Adipose - Subcutaneous" = "#ef883a",                   
                                  "Adipose - Visceral (Omentum)" = "#ffa804",
                                  "Breast - Mammary Tissue" = "#30ccd1"))

# creating annoation with ann df and colours list 

colAnn_fat <- HeatmapAnnotation(Tissues = ann_fat,
                                   which = 'col',
                                   col = colours_fat,
                                   annotation_width = unit(c(1, 4), 'cm'),
                                   gap = unit(1, 'mm'))

# dendogram 

row_dend_fat <- as.dendrogram(hclust(dist(as.matrix(prob_htr[, c(1,2,21)]))))

# heatmap 

pdf(qq("fat.pdf"), width = 9, height = 9)

ht_fat <- Heatmap(as.matrix(prob_htr[, c(1,2,21)]), 
                  col = col_fun,
                  clustering_distance_columns =  "spearman",
                  clustering_method_rows = "complete",
                  clustering_method_columns  = "complete",
                  cluster_rows = color_branches(row_dend_fat, k = 4),
                  name = "value",
                  column_title = "HTR Suptype Distribution in Fatty Tissue",
                  show_column_names = ,
                  column_names_gp =  gpar(fontsize = 9),
                  top_annotation = colAnn_fat)

draw(ht_fat)

dev.off()


###############################################################################

# correlation

# correlation of tissues

cor_M1 <- cor(M1, M1, method = "spearman")  # corr
        
P1 <- cor.mtest(M1, conf.level = 0.95) # significance 
        
# pdf(qq("corr_tissue.pdf"), width = 9, height = 9)
        
corrplot(cor_M1, 
  method = 'color',  # square 
  order = 'hclust', # sorted by clustering
  tl.col = 'black', # labels are black
  tl.cex = 0.5, # size of the labels
  cl.pos = 'b', # annotation bar at bottom
  cl.cex = 0.5, # arranging the text size to fit,
  title = "Tissues clustered by HTR distribution",
  mar=c(0,0,1,0),
  col = brewer.pal(n = 10, name = 'RdYlBu')) # colour scheme 
        
# dev.off()


# correlation of receptors 
        
        
cor_M2 <- cor(M2, M2, method = "pearson") # corr
        
P2 <- cor.mtest(M2, conf.level = 0.95) # sig 
        
# pdf(qq("corr_htr.pdf"), width = 9, height = 9)
        
corrplot(cor_M2, 
  method = 'color', # square 
  order = 'hclust', # sorted by clustering
  p.mat = P2$p, # p values
  sig.level = c(0.001, 0.01, 0.05), # levels of asterisk 
  pch.cex = 0.9, # size of the p vales 
  insig = 'label_sig', # label
  tl.col = 'black', # labels are black 
  tl.cex = 0.7, # size of the labels 
  cl.pos = 'b', # annotation bar at bottom 
  title = "HTR clustered by tissue distribution",
  mar=c(0,0,1,0),
  col = brewer.pal(n = 10, name = 'RdYlBu')) # color scheme 
        
        
# dev.off()
        
