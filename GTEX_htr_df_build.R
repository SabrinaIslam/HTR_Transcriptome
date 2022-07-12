################################################################################

# packages 

################################################################################

library(dplyr) # for wrangling data frames 

library(data.table) # fread

library(tidyverse) # tidy data 

library(matrixStats) # matrix stats 


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

# shrorten the dataframe just to have the HTR genes 

################################################################################

# making df 

# meta is the gtex meta data 

meta <- read.csv("Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                 header = T, 
                 sep = "\t",     
                 quote = "")

# indexing the rna seq samples 

rnaseq <- meta %>% 
  filter(SMAFRZE == "RNASEQ") %>% 
  dplyr::select(SAMPID)

dim(rnaseq)

# read table 

i_read <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
                select = c("Name","Description", rnaseq),
                stringsAsFactors = F,
                header = T, 
                sep = "\t", 
                quote = "",
                check.names = F,
                skip = 2)  

# htr genes 

htrgenes <- c("ENSG00000158748.3",
              "ENSG00000179546.4",
              "ENSG00000135914.5",
              "ENSG00000179097.5",
              "ENSG00000178394.4",
              "ENSG00000164270.17",
              "ENSG00000135312.6",
              "ENSG00000168830.7",
              "ENSG00000157219.3",
              "ENSG00000148680.15",
              "ENSG00000102468.10",
              "ENSG00000147246.9")

# filtered just the htr genes from the read table  

htr_read <- data.frame(i_read %>%
                         filter (i_read$Name %in% htrgenes))


htr_read$Description # checks 

rownames(htr_read) <- htr_read$Name # set rownames 

rownames(htr_read) #checks 

# writing csv 

write.csv(htr_read, "Data/GTEX_HTR.csv") # data

# reading the data as a test 

htr_exp <- read.csv("Data/GTEX_HTR.csv", 
                    row.names = 1, 
                    header = T)

rownames(htr_exp) # checks 

colnames(htr_exp) # checks 

# making annotation 

sample_names <- colnames(htr_exp[, -c(1,2)]) # are the gene annotations 

sample_names <- gsub("\\.", "-", sample_names) # replacing the "." with "-"

sample_tissues <- meta %>% 
  filter(SAMPID %in% sample_names) %>%
  dplyr::select(c("SAMPID", "SMTSD")) # selecting these samples from meta data 

head(sample_tissues$SAMPID) #checks 

write.csv(sample_tissues, "Data/SampleTissues.csv") # writes as a file 


# statistics 

htr_exp <- as.matrix(htr_exp[, c(3:17384)])

head(htr_exp)

quantile(htr_exp, probs = c(0.1, 0.25, 0.5, 0.75, 0.95, 0.99, 1)) #50% under 0.1

mean(htr_exp) # 1.17

# min without 0 

for (i in 1:12) {
  
  htr_mins <-  min(htr_exp[i,][htr_exp[i,] != 0])   
  print(htr_mins)
}

htr_mins <- c(0.005414,
              0.004105,
              0.009073,
              0.005277,
              0.006509,
              0.0009282,
              0.003172,
              0.003594,
              0.003133,
              0.003555,
              0.003492,
              0.001683)

htr_max <- rowMaxs(htr_exp)

htr_median <- rowMedians(htr_exp)

htr_means <- rowMeans(htr_exp)

htr_sds <- rowSds(htr_exp)

# makinmatrix

htr_stats <- as.data.frame(cbind(htr_mins, htr_max, htr_means, htr_median, htr_sds))

htr_stats <- rename(htr_stats, 
                    min = htr_mins, 
                    max = htr_max, 
                    mean = htr_means, 
                    median = htr_median, 
                    SD = htr_sds )


htr_stats

