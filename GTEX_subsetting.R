################################################################################

# packages 

################################################################################

# statistics 

library(matrixStats)

# data manipulation 

library(dplyr) # for wrangling data frames 
library(data.table) # fread
library(tidyverse) # tidy data 

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

# splitting the GTEX data into system-specific data sets 


######################################################################

# reading metadata

########################################################################## 

## header = F to make sure the top row does not become column names

## t separated

meta <- read.csv("Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                 header = F, 
                 sep = "\t",     
                 quote = "")

colnames(meta)

## has a weird column header labeled V1 to V63 

## the second row should be column header

## making the right column header 

  ## this row will become the header 

  colnames(meta) <- as.character(unlist(meta[1,]))

  ## removing the header as a row 

  meta <- meta[-1, ]

##############################################################################
  
# subset criterion
  
############################################################################### 
 
   # tissues as vectors 

  cns <- c('Brain - Putamen (basal ganglia)',
         'Brain - Nucleus accumbens (basal ganglia)',
         'Brain - Caudate (basal ganglia)',
         'Brain - Cerebellum',
         'Brain - Cerebellar Hemisphere',
         'Brain - Anterior cingulate cortex (BA24)',
         'Brain - Frontal Cortex (BA9)',
         'Brain - Cortex',
         'Brain - Hypothalamus',
         'Brain - Hippocampus',
         'Brain - Amygdala',
         'Brain - Substantia nigra',
         'Brain - Spinal cord (cervical c-1)',
         'Nerve - Tibial') #14 

  ens <- c('Small Intestine - Terminal Ileum',
         'Colon - Sigmoid',
         'Colon - Transverse',
         'Esophagus - Muscularis',
         'Esophagus - Mucosa',
         'Esophagus - Gastroesophageal Junction',
         'Stomach',
         'Minor Salivary Gland',
         'Liver',
         "Spleen", 
         "Pancreas") # 11

  cv <- c('Artery - Tibial',
        'Artery - Coronary',
        'Artery - Aorta',
        'Heart - Atrial Appendage',
        'Heart - Left Ventricle') # 5

  endocrine <- c('Pituitary',
               'Thyroid',
               'Adrenal Gland') # 3
  
  fem_reproductive <- c("Uterus",
    "Cervix - Endocervix",
    "Fallopian Tube",
    "Cervix - Ectocervix",
    "Ovary",
    "Vagina") # 6
  
  mus_reproductive <- c("Prostate",
                        "Testis") #2
  
  excretory <- c("Bladder",
    "Kidney - Medulla",
    "Kidney - Cortex") # 3
  
  adipose <- c("Adipose - Subcutaneous",
    "Breast - Mammary Tissue",
    "Adipose - Visceral (Omentum)") # 3
  
  skin <- c("Skin - Not Sun Exposed (Suprapubic)",
            "Skin - Sun Exposed (Lower leg)") #2

  mus <- c("Muscle - Skeletal") #1 
  
  resp <- c("Lung") #1 
  
  lymphocyte <- c("Cells - EBV-transformed lymphocytes") #1 
  
  whole_blood <- c("Whole Blood") #1 
  
  fibroblast <- c("Cells - Cultured fibroblasts") # 1
  
  
# sample ids of those tissues + RNA seq data 

    ## rnaseq will be the SMTSD %in% tissues I am interested and only RNA seq    

  ## cns 
  
  rnaseq_cns <- meta %>% 
  filter(SMTSD %in% cns, SMAFRZE == "RNASEQ") %>% 
  dplyr::select(SAMPID)
  
  dim(rnaseq_cns) # 3261 cns samples 

  ## ens 
  
  rnaseq_ens <- meta %>% 
    filter(SMTSD %in% ens, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_ens) # 3727 ens samples 
  
  ## cv 
  
  rnaseq_cv <- meta %>% 
    filter(SMTSD %in% cv, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)

  dim(rnaseq_cv) # 2196 samples 
  
  # endocrine
  
  rnaseq_endo <- meta %>% 
    filter(SMTSD %in% endocrine, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_endo) # 1194 
  
  # female repro 
  
  rnaseq_fem_repro <- meta %>% 
    filter(SMTSD %in% fem_reproductive, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_fem_repro) # 506
  
  # male repro
   
  rnaseq_mus_repro <- meta %>% 
    filter(SMTSD %in% mus_reproductive, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_mus_repro) # 606
  
  # excretory 
  
  rnaseq_excre <- meta %>% 
    filter(SMTSD %in% excretory, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_excre) # 110
  
  # adipose
  
  rnaseq_adi <- meta %>% 
    filter(SMTSD %in% adipose, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_adi) # 1663
  
  # skin 
  
  rnaseq_skin <- meta %>% 
    filter(SMTSD %in% skin, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_skin) # 1305
  
  # muscle
  
  rnaseq_muscle <- meta %>% 
    filter(SMTSD %in% mus, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_muscle) # 803
  
  
  # resp <- c("Lung") 
  
  rnaseq_lung <- meta %>% 
    filter(SMTSD %in% resp, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_lung) # 578
  
  
  # lymphocyte  
  
  rnaseq_lcyte <- meta %>% 
    filter(SMTSD %in% lymphocyte, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_lcyte) # 174
  
  
  # whole_blood  
  
  rnaseq_blood <- meta %>% 
    filter(SMTSD %in% whole_blood, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_blood) # 775
  
  # fibroblast 
  
  rnaseq_fblast <- meta %>% 
    filter(SMTSD %in% fibroblast, SMAFRZE == "RNASEQ") %>% 
    dplyr::select(SAMPID)
  
  dim(rnaseq_fblast) # 504
  
################################################################################

  ## reading the GTEX data     
  
############################################################################### 
  
  
  ## header false 
  
  ## fread allows fishing out target columns 
  
  ## skip the first two rows 
  
  ## fread function 


  ## warning: no method or default for coercing "character" to "SAMPID"
  
  # cns
  
  read_cns <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
              select = c("Name","Description", rnaseq_cns), 
              stringsAsFactors = F,
              header = T, 
              sep = "\t", 
              quote = "",
              skip = 2)  
  
  dim(read_cns) #3261
 
  # ens 
  
  read_ens <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                    select = c("Name","Description", rnaseq_ens), 
                    stringsAsFactors = F,
                    header = T, 
                    sep = "\t", 
                    quote = "",
                    skip = 2)  

  dim(read_ens) ## 3727
  
  
  # cv 
  
  read_cv <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                    select = c("Name","Description", rnaseq_cv), 
                    stringsAsFactors = F,
                    header = T, 
                    sep = "\t", 
                    quote = "",
                    skip = 2)  
  
  dim(read_cv) # 2196
  
  # endo
  
  read_endo <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                   select = c("Name","Description", rnaseq_endo), 
                   stringsAsFactors = F,
                   header = T, 
                   sep = "\t", 
                   quote = "",
                   skip = 2)  
  
  dim(read_endo) # 1194
  
  # repro fem 
  
  read_fem_repro <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                     select = c("Name","Description", rnaseq_fem_repro), 
                     stringsAsFactors = F,
                     header = T, 
                     sep = "\t", 
                     quote = "",
                     skip = 2)  
  
  dim(read_fem_repro) # 506
  
  # repro mus 
  
  read_mus_repro <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                      select = c("Name","Description", rnaseq_mus_repro), 
                      stringsAsFactors = F,
                      header = T, 
                      sep = "\t", 
                      quote = "",
                      skip = 2)  
  
  dim(read_mus_repro) # 606
  
  # excre
  
  read_excre <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                     select = c("Name","Description", rnaseq_excre), 
                     stringsAsFactors = F,
                     header = T, 
                     sep = "\t", 
                     quote = "",
                     skip = 2)  
  
  dim(read_excre) # 110
  
  # adi
  
  read_adi <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                     select = c("Name","Description", rnaseq_adi), 
                     stringsAsFactors = F,
                     header = T, 
                     sep = "\t", 
                     quote = "",
                     skip = 2)  
  
  dim(read_adi) # 1663

  # lung
  
  read_lung <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                      select = c("Name","Description", rnaseq_lung), 
                      stringsAsFactors = F,
                      header = T, 
                      sep = "\t", 
                      quote = "",
                      skip = 2)  
  
  dim(read_lung) # 578
  
  # muscle 
  
  read_muscle <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                      select = c("Name","Description", rnaseq_muscle), 
                      stringsAsFactors = F,
                      header = T, 
                      sep = "\t", 
                      quote = "",
                      skip = 2)  
  
  dim(read_muscle) # 803
  
  # skin
  
  read_skin <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                      select = c("Name","Description", rnaseq_skin), 
                      stringsAsFactors = F,
                      header = T, 
                      sep = "\t", 
                      quote = "",
                      skip = 2)  
  
  dim(read_skin) # 1307
  
  
  # lymphocyte 
  
  read_lcyte <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                      select = c("Name","Description", rnaseq_lcyte), 
                      stringsAsFactors = F,
                      header = T, 
                      sep = "\t", 
                      quote = "",
                      skip = 2)  
  
  dim(read_lcyte) # 174
  
  str(read_lcyte)
  
  # whole blood 
  
  read_blood <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                      select = c("Name","Description", rnaseq_blood), 
                      stringsAsFactors = F,
                      header = T, 
                      sep = "\t", 
                      quote = "",
                      skip = 2)  
  
  dim(read_blood) # 755
  
  # fibroblast 
  
  read_fblast <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                      select = c("Name","Description", rnaseq_fblast), 
                      stringsAsFactors = F,
                      header = T, 
                      sep = "\t", 
                      quote = "",
                      skip = 2)  
  
  dim(read_fblast) # 504 
  
  # writing csv 
  
  write.csv(read_cns, "Data/GTEXCNS.csv", sep = "\t") #cns 
  
  write.csv(read_ens, "Data/GTEXENS.csv", sep = "\t") #ens
  
  write.csv(read_cv, "Data/GTEXCV.csv", sep = "\t") # cv
  
  write.csv(read_endo, "Data/GTEXENDO.csv", sep = "\t") # endo 
  
  write.csv(read_fem_repro, "Data/GTEXFEM.csv", sep = "\t") # fem
  
  write.csv(read_mus_repro, "Data/GTEXMUS.csv", sep = "\t") # mus 
  
  write.csv(read_excre, "Data/GTEXEXCRE.csv", sep = "\t") # excre
  
  write.csv(read_adi, "Data/GTEXADI.csv", sep = "\t") # adi 
  
  write.csv(read_muscle, "Data/GTEXMUSCLE.csv", sep = "\t") # muscle
  
  write.csv(read_skin, "Data/GTEXSKIN.csv", sep = "\t") # skin 
  
  write.csv(read_lung, "Data/GTEXLUNG.csv", sep = "\t") # lung 
  
  write.csv(read_blood, "Data/GTEXBLOOD.csv", sep = "\t") # blood
  
  write.csv(read_lcyte, "Data/GTEXLCYTE.csv", sep = "\t") # lcyte
  
  write.csv(read_fblast, "Data/GTEXFBLAST.csv", sep = "\t") # fibroblast  
  