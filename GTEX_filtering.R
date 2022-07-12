
################################################################################

# packages 

################################################################################

# data manipulation 

library(dplyr) # for wrangling data frames 
library(data.table) # fread
library(tidyverse) # tidy data 

# visualisation 

library(ggplot2) # plotting 
library(gplots) # plotting data 
library(RColorBrewer) # build color-pallates for plots 
library(styler) # codes looks nicer
library(ggthemes) # themes 
library(viridis) # colour 

# statistics 

library(spgs) # tuning point test 
library(rstatix)
library(reshape)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(plyr)
library(datarium)
library(paletteer)


## RNA-seq specific packages 

library(limma) # for expression data 
library(edgeR) # for RNA- seq data 
library(stringr) # to create a matrix using split
library(affy) # plotDensity polynomial fitted plots  


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

# data: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz 

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

# setting the threshold for detection  

################################################################################

# importing data 

###############################################################################

# row number 1 is header which is true, 

# check.names set to F to make sure names are not changed by R

gtex_count <- read.csv("Data/bigGTEX.csv",
                    header = TRUE,
                    check.names = F,
                    row.names = 1)


gtex_tissue <- read.csv("Data/bigGTEXTissues.csv",
                           header = TRUE,
                           check.names = F,
                           row.names = 1)

# quick look

dim(gtex_count) # 56200 by 1280

str(gtex_count)

# making the ensemble ID rownames

rownames(gtex_count) <- gtex_count$Name

# making "Name" and "Description" factor 

gtex_count$Name <- as.factor(gtex_count$Name)

gtex_count$Description <- as.factor(gtex_count$Description)

# just the annotation : the Names and Descriptions 

gtex_anno <-  as.data.frame(gtex_count[,c(1,2)])

rownames(gtex_anno) <- rownames(gtex_count) # proper rownames 

# just the expression values

# make the data matrix numeric 

# gtex_exp <- data.frame(sapply(gtex_count, function(x) as.numeric(as.character(x))))

gtex_exp <- gtex_count[,-c(1,2)]

rownames(gtex_exp) <- rownames(gtex_count) # proper rownames 

class(gtex_exp) # data frame 

dim(gtex_exp) # 56200 by 779 

gtex_exp[1:3, 1:4]  


# very important lol, do we have the same samples 

colnames(gtex_count) %in%  gtex_tissue$SAMPID # T

colnames(gtex_exp) %in%  gtex_tissue$SAMPID # T

# color palette

plot_col <- brewer.pal(12, "Paired") # color palette "paired" 

plot_sample_col <- c("#ffbfb0",
                    "#00bf5d",
                    "#da2baa",
                    "#65a600",
                    "#c563eb",
                    "#0e7f92",
                    "#bf0016",
                    "#628eff",
                    "#ab5800",
                    "#553eb0",
                    "#fcb96f",
                    "#6b3f8a",
                    "#00a794",
                    "#ff4d58",
                    "#74d6d6",
                    "#a30852",
                    "#724a12",
                    "#ff69a7",
                    "#ae9754",
                    "#c49bcb",
                    "#2ca25f",
                    "#8856a7",
                    "#c994c7",
                    "#dd1c77",
                    "#f03b20")

length(plot_sample_col)

# long data

sum_gtex_exp <- data.frame(value=apply(gtex_exp,2,sum))

sum_gtex_exp$key <- rownames(sum_gtex_exp)

# barplot 

ggplot() +
geom_bar(data = sum_gtex_exp[25:44, ], 
         aes(x = key, y =value, fill = key),
         stat="identity") +
labs(title = "library size in millions",
     x = "Samples",
     y = "Library size") +
scale_fill_manual(values = plot_sample_col) +  
theme_tufte() +
theme(axis.line = element_line(),  
      axis.text.x = element_text(angle = 90, vjust=0.6, hjust=0),
      panel.background = element_blank(),
      legend.position = "null")

################################################################################
  
# transforming data   
  
################################################################################  

head(gtex_exp)

# making dge object 

gtex_dge <- DGEList(counts = gtex_exp, genes = gtex_anno)

names(gtex_dge) # sample, counts, and genes 

# the samples 

head(gtex_dge$samples) # names of the samples, group, lib.size, normalizing factor

length(gtex_dge$samples$lib.size) # 779

class(gtex_dge$samples$lib.size) # numeric 

# all normalizing factor is 1, lib.size `app. 1000000

# the genes 

head(gtex_dge$genes) # genes ensembl id and common name 

# the counts 

head(gtex_dge$counts) # the expression matrix

class(gtex_dge$counts) # matrix

# log2 of dge

log2_gtex_dge <- log2((gtex_dge$counts) + 1)


# checking library sizes of sample 77-100 with box plot

par(mfrow = c(1, 1))

barplot(gtex_dge$samples$lib.size[77:100],
      col=plot_sample_col,
      las = 2)
title("Barplot of library sizes")

  
###############################################################################
  
# sample statistics 

###############################################################################
  
# library statistics 
  
lib_mean <- mean(gtex_dge$samples$lib.size) * 1e-6

lib_median <- median(gtex_dge$samples$lib.size) * 1e-6

c(lib_mean, lib_median) # both mean and median are 1?


# normalsing factor of the library == 1 

norm_data <- as.data.frame(gtex_dge$samples[, -c(2,3)])

norm_data$sample <- rownames(gtex_dge$samples)

colnames(norm_data)


names(norm_data) <- c("norm.factors", "sample")

norm_data$norm.factor <- as.numeric(norm_data$norm.factors)

norm_data <- norm_data[c(2,1)]

hist(as.matrix(gtex_dge$samples$norm.factors), 
     col=plot_sample_col[3], 
     border="white",
     main="Normalising factors of the library", 
     xlab="noralising factors", 
     ylab="number", 
     las=1, 
     cex.axis=0.7) # 1

# all norm.factors are indeed 1

cutoff_library <- log2(10/lib_median + 2/lib_mean)

cutoff_library # 3.58

# count properties 

count_mean <- mean(log2_gtex_dge) # 1.31

count_median <- median(log2_gtex_dge) 0.065

c(count_mean, count_median)

# median of tissue means 

mean_tissue <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", 
                     stringsAsFactors = F,
                     header = T, 
                     sep = "\t", 
                     quote = "",
                     skip = 2)  

head(mean_tissue)

mean_tissue <- mean_tissue [, -c(1,2)] # removing "names" and "description"

head(mean_tissue)

class(mean_tissue$`Adipose - Subcutaneous`) # numeric 

row_sub <- apply(mean_tissue, 1, function(row) all(row !=0 )) # index 

mean_tissue <- mean_tissue[row_sub, ] # removes the 0

dim(mean_tissue)

my_median <- data.frame(apply(mean_tissue,2,median)) # median calculation 

my_median$Tissue <- rownames(my_median) # making median a column 

colnames(my_median)[which(names(my_median) == "apply.mean_tissue..2..median.")] <- "Median" 

my_median <- my_median[c(2,1)]

head(my_median)

ggplot(my_median, aes(Tissue, Median))+
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 8)) 

median(my_median$Median) # median of tissue means is 6.526132
  

# conclusion: mean and median of counts are very low, 

# possibly due to a lot of zeros 

################################################################################
  
# density 

################################################################################
  
# single sample 
  
par(mfrow = c(1,1))

plot(density(log2_gtex_dge[, 755]), 
col= plot_sample_col[1], 
xlab="logcpm",
ylab = "density",
lwd=2, 
ylim=c(0,2.5), 
las=2, 
main="Density plot of sample 755")
grid()
abline(v= cutoff_library, lty=3)
abline(v = count_median, col = "red")

# majority of genes are lowly expressed / unexpressed 

# majority of logCPM values are slight positive 

# whole data 

# we need the plot density function from affy 

whole_density <- plotDensity(log2_gtex_dge, 
        lty=1, 
        col= plot_sample_col[6], 
        xlab="logcpm",
        ylab = "density",
        lwd=2,
        ylim=c(0,15),
        main="Density plot of whole dataset")
grid()
abline(v= cutoff_library, lty=3)
abline(v = count_median, col = "red")

# same pattern although we do see some counts that can peak very high 

# its a polynomial fit,looks smoother than the reality 

# hard to form conlusions on setting threshold

# conclusion: we need to look tissue-by-tissue to get an idea 
 
################################################################################

# density by tissues 

################################################################################

levels(as.factor(gtex_tissue$SMTSD))

tissue_list <- levels(as.factor(gtex_tissue$SMTSD))

# poly fitted density plot 

for (i in 1:length(tissue_list)) 
{
 tissue_i <- gtex_tissue %>%
   filter (SMTSD == tissue_list[i]) %>%
   dplyr::select(SAMPID) %>% unlist  
 print(tissue_list[i])
 
 poly_i <- plotDensity(as.data.frame(log2_gtex_dge[, colnames(log2_gtex_dge) %in%  tissue_i]),
                       lty=1,
                       col= plot_sample_col[1],
                       xlab="logcpm",
                       ylab = "density",
                       lwd=2,
                       ylim=c(0,5),
                       main=sprintf( "Logcpm in %s Sample", tissue_list[i]))
 grid()
 abline(v= cutoff_library, lty=3)
 abline(v = count_median, col = "red")
 
 
 poly_i
}

# its not bimodal and the detection threshold is still hard to determine 

################################################################################

# investigating the count data 

################################################################################

# percentage of counts between intervals

head(summary(log2_gtex_dge)) # lots of min 0

length(log2_gtex_dge [log2_gtex_dge == 0]) / length(log2_gtex_dge) * 100 # 46.73 are 0 

length(log2_gtex_dge [log2_gtex_dge > 0]) / length(log2_gtex_dge) * 100 # 53.26 are values that are not 0 

count_median <- median(log2_gtex_dge)

counts_above_median <- length(log2_gtex_dge [log2_gtex_dge >= count_median]) / length(log2_gtex_dge) * 100

counts_above_median # 50% counts are above median 

quantile(log2_gtex_dge, probs = c(0.05, 0.95)) # 95% falls under 5.32 

valid_counts_above_median <- 
 length(log2_gtex_dge [log2_gtex_dge > 0]) / length(log2_gtex_dge) * 100 - (length(log2_gtex_dge [log2_gtex_dge >= count_median]) / length(log2_gtex_dge) * 100) 

valid_counts_above_median # 3.27 % that are not 0 above the median 

# how many rows are all zeroes?

# conclusion: we need to filter out the uninformative counts

# before we can detect the threshold 

################################################################################

 ### filtering less than 0 ###################

###############################################################################

gtex_thresh <- gtex_exp > 0

head(gtex_thresh) # logical 

gtex_keep <- (rowSums(gtex_thresh) > 25) # at least 25 TRUES i

table(gtex_keep)

summary(gtex_keep) # genes that have 0 in all samples 

gtex_exp_filter <- gtex_exp [gtex_keep, ]

dim(gtex_exp_filter)

# using the same index to create new anno, because current anno has rows for all genes  

## the row.names of the count_keep (not the column 1)

gtex_keep_genes <- row.names(gtex_exp_filter)

gtex_anno_keep <- gtex_anno %>% rownames() %in% gtex_keep_genes 

table(gtex_anno_keep) # 49497  T

gtex_anno_filter <- as.data.frame(gtex_anno[gtex_anno_keep, ])

head(gtex_anno_filter)

# new dge here 

gtex_fil_dge <- DGEList(counts = gtex_exp_filter, genes = gtex_anno_filter)

names(gtex_fil_dge) # sample, counts, and genes 

# log2 of dge

log2_gtex_fil_dge <- log2((gtex_fil_dge$counts) + 1)

dim(log2_gtex_fil_dge) #49497, 779 

100 - ((nrow(log2_gtex_fil_dge)/nrow(log2_gtex_dge))*100) #11.9% filtered 
  
# conclusion: not enough filtering 

##############################################################################
  
### filter by expression with edgeR 
  
###############################################################################
  
# filter genes, while keeping as many genes as possible with worthwhile counts

samples <- gtex_tissue

head(samples)

# a design matrix just on tissue types 

samples$SMTSD <- as.factor(samples$SMTSD)

samples$SAMPID <- as.factor(samples$SAMPID)

design <- model.matrix(~0 + samples$SMTSD)

colnames(design) <- levels(samples$SMTSD)

design

# filtering 

keep <- filterByExpr(gtex_fil_dge, design)

keep_gtex_fil_dge <- gtex_fil_dge[keep, , keep.lib.sizes = FALSE]

dim(keep_gtex_fil_dge$counts) # 23277

log2_keep_gtex_fil_dge <- log2((keep_gtex_fil_dge$counts) + 1)

dim(log2_keep_gtex_fil_dge)

100 - ((nrow(log2_keep_gtex_fil_dge)/nrow(log2_gtex_dge))*100) # 58.858 % removed when it should be 30% 

# plotting density plots to see the effect of filtering 

tissue_list_1 <- levels(as.factor(gtex_tissue$SMTSD))[1:27] # from adipose to breast 

tissue_list_2 <- levels(as.factor(gtex_tissue$SMTSD))[28:54] # from cervix to vagina 

for (i in 1:length(tissue_list_1)) 
{
tissue_i_1 <- gtex_tissue %>%
  filter (SMTSD == tissue_list_1[i]) %>%
  dplyr::select(SAMPID) %>% unlist  
print(tissue_list_1[i])

poly_i <- plotDensity(as.data.frame(log2_gtex_dge[, colnames(log2_gtex_dge) %in%  tissue_i_1]),
                      lty=1,
                      col= plot_sample_col[1],
                      xlab="log2(tpm+1)",
                      ylab = "density",
                      lwd=2,
                      ylim=c(0,5),
                      main=sprintf( "Logcpm in %s Sample", tissue_list_1[i]))
grid()
abline(v = 1.15, col = "red")


poly_i

poly_i <- plotDensity(as.data.frame(log2_gtex_fil_dge[, colnames(log2_gtex_fil_dge) %in%  tissue_i_1]),
                      lty=1,
                      col= plot_sample_col[1],
                      xlab="log2(tpm+1)",
                      ylab = "density",
                      lwd=2,
                      ylim=c(0,5),
                      main=sprintf( "Logcpm in Filtered %s Sample", tissue_list_1[i]))
grid()
abline(v = count_median, col = "red")


poly_i

poly_i <- plotDensity(as.data.frame(log2_keep_gtex_fil_dge[, colnames(log2_keep_gtex_fil_dge) %in%  tissue_i_1]),
                      lty=1,
                      col= plot_sample_col[1],
                      xlab="log2(tpm+1)",
                      ylab = "density",
                      lwd=2,
                      ylim=c(0,0.6),
                      main=sprintf( "Logcpm in filterByExpression Filtered %s Sample", tissue_list_1[i]))
grid()
abline(v = 1.15, col = "red")

poly_i
}



for (i in 1:length(tissue_list_2)) 
{
tissue_i_2 <- gtex_tissue %>%
  filter (SMTSD == tissue_list_2[i]) %>%
  dplyr::select(SAMPID) %>% unlist  
print(tissue_list_2[i])

poly_i <- plotDensity(as.data.frame(log2_gtex_dge[, colnames(log2_gtex_dge) %in%  tissue_i_2]),
                      lty=1,
                      col= plot_sample_col[1],
                      xlab="logcpm",
                      ylab = "density",
                      lwd=2,
                      ylim=c(0,5),
                      main=sprintf( "Logcpm in %s Sample", tissue_list_2[i]))
grid()


poly_i

poly_i <- plotDensity(as.data.frame(log2_gtex_fil_dge[, colnames(log2_gtex_fil_dge) %in%  tissue_i_2]),
                      lty=1,
                      col= plot_sample_col[1],
                      xlab="logcpm",
                      ylab = "density",
                      lwd=2,
                      ylim=c(0,5),
                      main=sprintf( "Logcpm in Filtered %s Sample", tissue_list_2[i]))
grid()
abline(v = 1.15, col = "red")


poly_i

poly_i <- plotDensity(as.data.frame(log2_keep_gtex_fil_dge[, colnames(log2_keep_gtex_fil_dge) %in%  tissue_i_2]),
                      lty=1,
                      col= plot_sample_col[1],
                      xlab="logcpm",
                      ylab = "density",
                      lwd=2,
                      ylim=c(0,0.6),
                      main=sprintf( "Logcpm in filterByExpression Filtered %s Sample", tissue_list_2[i]))
grid()
abline(v = 1.15, col = "red")


poly_i
}

# DT is 1.15

# variable tissues: 

  # thyroid samples, stomach, small intestine, liver, kidney, left ventricle, atrial appendage, eso mus, eso gi, 
  # tranverse colon, sigmoid colon, ectocervix, breast, substantia nigra
  # putamen, nucleus accubens, hypothalamus, hippocampus, frontal cortex, cortex, cer hem, caudate basal ganglia, 
  # acc, amugdala, adrenal 

# conclusion: a lot of the genes are getting filtered by edgeR

#################################################################################

# why so much filtered? 
  
################################################################################  

ordinary_fil <- 19000/ 28000 # 68%
  
my_data_fil <- nrow(log2_keep_gtex_fil_dge)/nrow(log2_gtex_dge) # 41%

ratio <- 56200/ 28000 # 2

if_my_was_ordinary <- nrow(log2_keep_gtex_fil_dge)/28000 # 83 % 

# because I have too many empty genes 
  
###################################################################################
  
# proof with new filtering
  
# filtering keep from 1.15

gtex_thresh_1.15 <- log2_keep_gtex_fil_dge > 1.15

head(gtex_thresh_1.15) # logical 

gtex_keep_1.15 <- (rowSums(gtex_thresh_1.15) > 250) # at least 1/4 samples TRUES i

table(gtex_keep_1.15)

summary(gtex_keep_1.15) # genes that have 0 in all samples 

log2_keep_gtex_1.15_dge <- log2_keep_gtex_fil_dge [gtex_keep_1.15, ]

dim(log2_keep_gtex_1.15_dge) # 18050

# using the same index to create new anno, because current anno has rows for all genes  

## the row.names of the count_keep (not the column 1)

gtex_keep_genes_1.15 <- row.names(log2_keep_gtex_1.15_dge)

gtex_anno_keep_1.15 <- gtex_anno %>% rownames() %in% gtex_keep_genes_1.15 

table(gtex_anno_keep_1.15) # 37208  T

gtex_anno_filter_1.15 <- as.data.frame(gtex_anno[gtex_anno_keep_1.15, ])

head(gtex_anno_filter_1.15)

100 - ((nrow(log2_keep_gtex_fil_dge)/nrow(log2_gtex_dge))*100) # 58.89 % filtered 

100 - ((nrow(log2_keep_gtex_1.15_dge)/nrow(log2_gtex_dge))*100) # 67.68 % filtered 

nrow(log2_keep_gtex_1.15_dge) - nrow(log2_keep_gtex_fil_dge)

for (i in 1:length(tissue_list)) 
{
tissue_i <- gtex_tissue %>%
  filter (SMTSD == tissue_list[i]) %>%
  dplyr::select(SAMPID) %>% unlist  
print(tissue_list[i])

poly_i <- plotDensity(as.data.frame(log2_keep_gtex_1.15_dge[, colnames(log2_keep_gtex_1.15_dge) %in%  tissue_i]),
                      lty=1,
                      col= plot_sample_col[1],
                      xlab="log2(tpm+1)",
                      ylab = "density",
                      lwd=2,
                      ylim=c(0,0.6),
                      main=sprintf( "Logtpm in %s Sample after >1.15", tissue_list[i]))
grid()

poly_i
}

  
# we get normal curves

################################################################################

# looking at individual libraries rather than polynomial fit

################################################################################
  
par(mfrow = c(1,1))

for (i in 1:length(tissue_list)) 

{
tissue_i <- gtex_tissue %>%
  filter (SMTSD == tissue_list[i]) %>%
  dplyr::select(SAMPID) %>% unlist  

print(tissue_list[i])

samplenames_i <- as.vector(as.character(tissue_i))

nsamples_i <- (length(tissue_i))

print(samplenames_i)

plot(density(as.data.frame(log2_keep_gtex_1.15_dge[, colnames(log2_keep_gtex_1.15_dge) %in%  tissue_i])[,1]),
     col=plot_sample_col[1],
     lwd=2,
     xlab="logcpm",
     ylab = "density",
     ylim=c(0,0.8),
     las=2,
     main=sprintf( "Logcpm in %s Sample after >1.15", tissue_list[i]))

for (j in 2:nsamples_i)
{
  den_j <- density(as.data.frame(log2_keep_gtex_1.15_dge[, colnames(log2_keep_gtex_1.15_dge) %in%  tissue_i])[,j])
  lines(den_j$x,
        den_j$y,
        col=plot_sample_col[j],            
        lwd=2)
}

legend("topright", samplenames_i, text.col= plot_sample_col, bty="n")

}

  
#################################################################################
  
  # vooms 
  
################################################################################    

# the probability of getting r events in a large sample. 
# in rna seq,  the probability of pulling out a particular transcript is very small
# It is appropriate for data where mean == variance.

# decreasing trend between the means and variances 
# variance values plateau at high expression values

par(mfrow = c(2,2))

v1 <- voom(log2_gtex_dge, design, plot=TRUE)
v1



v2 <- voom(log2_gtex_fil_dge, design, plot=TRUE)
v2

v3 <- voom(log2_keep_gtex_fil_dge, design, plot=TRUE)
v3

v4 <- voom(log2_keep_gtex_1.15_dge, design, plot=TRUE)
v4
  
# the mean-variance plot flattens out

# conclusion: 1.15 is indeed the DT 
  
################################################################################
  
# boxplots of counts 

################################################################################  
  
par(mfrow = c(1,1))

for (i in 1:54) 
{
tissue_i <- gtex_tissue %>%
  filter (SMTSD == tissue_list[i]) %>%
  dplyr::select(SAMPID) %>% unlist  

print(tissue_list[i])

samp_i <- as.data.frame(log2_gtex_dge[, colnames(log2_gtex_dge) %in%  tissue_i])

long_samp_i <- samp_i %>%
  gather(key = "Sample", value = "Value") %>%
  convert_as_factor(Sample)

box_i <- long_samp_i %>%
  ggplot() +
  geom_boxplot(aes(x = Sample, 
                   y = Value, 
                   fill = Sample),
               outlier.shape = NA) +
  scale_y_continuous(breaks =seq(0, 10, .5), limit = c(0, 10)) +
  coord_cartesian(xlim = NULL, ylim = c(0, 2.5)) +
  labs(title = sprintf( "%s Sample Counts", tissue_list[i]),
       x = "Samples",
       y = "Log2 CPM") +
  scale_fill_manual(values = plot_sample_col) +
  # geom_hline(yintercept = median(long_samp_i$Value), linetype ="dashed", color = "green", size = 1) +
  theme_tufte() +
  theme(axis.line = element_line(),
        axis.text.x = element_text(angle = 90, vjust=0.6, hjust=0, size = 10),
        axis.title.y.left = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        legend.position = "null")

print(box_i)

# median_i <- data.frame(apply(samp_i, 2, median))
# 
# var_i <- data.frame(apply(samp_i, 2, var))
# 
# mat_i <- cbind(median_i, var_i)
# 
# print(mat_i)
}

# samples are not normliased   
  
################################################################################
  
  # min, max, sd statsitics 

################################################################################  
  
for (i in 1:54) 
{
tissue_i <- gtex_tissue %>%
  filter (SMTSD == tissue_list[i]) %>%
  dplyr::select(SAMPID) %>% unlist  

print(tissue_list[i])

samp_i <- as.data.frame(log2_gtex_dge[, colnames(log2_gtex_dge) %in%  tissue_i])

median_i <- data.frame(apply(samp_i, 2, median))

var_i <- data.frame(apply(samp_i, 2, var))

max_i <- data.frame(apply(samp_i, 2, max))

min_i <- data.frame(apply(samp_i, 2, min))

mat_i <- cbind(min_i, median_i, max_i, var_i)

colnames(mat_i) <- c("minimum", "median", "max", "variance")

round(mat_i, digits = 2)

print(mat_i)

}

# a bonus check on libraries: they are all variable 


  
  