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

# making a small data set with 25 randomly selected sample from each tissue  


######################################################################

############## reading metadata ###################################

########################################################################## 

## header = F to make sure the top row doesnt become columm names

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

colnames(meta) # nice 

levels(as.factor(meta$SMTSD))

###############################################################################

##### reading data #################

###############################################################################


## The header is set to TRUE so that the first row is read as the header row 

## ((header is set to TRUE if and only if the first row contains one fewer field than the number of columns))

## The row.names is set to 1 so that the first column in read as the column with the the row-IDs. 

## check.names=FALSE to make sure "-" don't get replaced with "."

# ######## CNS #############

  cns_data <- read.csv("Data/GTEXCNS.csv",
                       header = TRUE, 
                       check.names=FALSE,
                       row.names = 1)
  
  head(cns_data) # checks 
  
  colnames(cns_data) # correct 
  
  str(cns_data) # numeric variables except Name and Description are characters 
  
 
 ############### ENS #######################
  
  
  ens_data <- read.csv("Data/GTEXENS.csv",
                       header = TRUE, 
                       check.names=FALSE,
                       row.names = 1)
  
  head(ens_data) # checks 
  
  colnames(ens_data) # correct 
  
  str(ens_data) # numeric variables except Name and Description are characters 
  

 ################# cv ###################
  
  cv_data <- read.csv("Data/GTEXCV.csv",
                      header = TRUE,
                      check.names=FALSE,
                      row.names = 1)
  
  head(cv_data) # checks 
  
  colnames(cv_data) # correct 
  
  str(cv_data) # numeric variables except Name and Description are characters 
  
  ############### ENdo #######################
  
  endo_data <- read.csv("Data/GTEXENDO.csv",
                       header = TRUE, 
                       check.names=FALSE,
                       row.names = 1)
  
  
  head(endo_data) # checks 
  
  colnames(endo_data) # correct 
  
  str(endo_data) # numeric variables except Name and Description are characters 
  
  
  ############### FEM #######################
  
  femm_repro_data <- read.csv("Data/GTEXFEM.csv",
                        header = TRUE, 
                        check.names=FALSE,
                        row.names = 1)
  
  head(femm_repro_data) # checks 
  
  colnames(femm_repro_data) # correct 
  
  str(femm_repro_data) # numeric variables except Name and Description are characters 
  
  
  ############### mus #######################
  
  mus_repro_data <- read.csv("Data/GTEXMUS.csv",
                              header = TRUE, 
                              check.names=FALSE,
                              row.names = 1)
  
  head(mus_repro_data) # checks 
  
  colnames(mus_repro_data) # correct 
  
  str(mus_repro_data) # numeric variables except Name and Description are characters 
  
  ############### EXCRE #######################
  
  excre_data <- read.csv("Data/GTEXEXCRE.csv",
                        header = TRUE, 
                        check.names=FALSE,
                        row.names = 1)
  
  
  head(excre_data) # checks 
  
  colnames(excre_data) # correct 
  
  str(excre_data) # numeric variables except Name and Description are characters
  

  ############### ADI #######################
  
  adi_data <- read.csv("Data/GTEXADI.csv",
                        header = TRUE, 
                        check.names=FALSE,
                        row.names = 1)
  
  
  head(adi_data) # checks 
  
  colnames(adi_data) # correct 
  
  str(adi_data) # numeric variables except Name and Description are characters 
  
  ############### skin #######################
  
  skin_data <- read.csv("Data/GTEXSKIN.csv",
                       header = TRUE, 
                       check.names=FALSE,
                       row.names = 1)
  
  
  head(skin_data) # checks 
  
  colnames(skin_data) # correct 
  
  str(skin_data) # numeric variables except Name and Description are characters 

  ############### lung #######################
  
  lung_data <- read.csv("Data/GTEXLUNG.csv",
                        header = TRUE, 
                        check.names=FALSE,
                        row.names = 1)
  
  
  head(lung_data) # checks 
  
  colnames(lung_data) # correct 
  
  str(lung_data) # numeric variables except Name and Description are characters 
  
  ############### muscle #######################
  
  muscle_data <- read.csv("Data/GTEXMUSCLE.csv",
                        header = TRUE, 
                        check.names=FALSE,
                        row.names = 1)
  
  
  head(muscle_data) # checks 
  
  colnames(muscle_data) # correct 
  
  str(muscle_data) # numeric variables except Name and Description are characters 
 
  ############### blood #######################
  
  blood_data <- read.csv("Data/GTEXBLOOD.csv",
                        header = TRUE, 
                        check.names=FALSE,
                        row.names = 1)
  
  
  head(blood_data) # checks 
  
  colnames(blood_data) # correct 
  
  str(blood_data) # numeric variables except Name and Description are characters 
  
   
  ############### lymphocyte #######################
  
  lcyte_data <- read.csv("Data/GTEXLCYTE.csv",
                        header = TRUE, 
                        check.names=FALSE,
                        row.names = 1)
  
  
  head(lcyte_data) # checks 
  
  colnames(lcyte_data) # correct 
  
  str(lcyte_data) # numeric variables except Name and Description are characters 
  
  ############### fibroblast #######################
  
  fblast_data <- read.csv("Data/GTEXFBLAST.csv",
                        header = TRUE, 
                        check.names=FALSE,
                        row.names = 1)
  
  
  head(fblast_data) # checks 
  
  colnames(fblast_data) # correct 
  
  str(fblast_data) # numeric variables except Name and Description are characters 
  
##################################################################  
  
########## subsetting from tissues #############################
  
###############################################################
  
  #######  'Brain - Cortex' ######################
  
  # all cortex samples 
  
    cortex <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Cortex') %>%
    dplyr::select("SAMPID") 
    
  dim(cortex) # 255 
  
  # random 25 sample vector  
  
  sampleindex_cortex <- sample(cortex$SAMPID, 25)
  
  sample_cortex <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_cortex))
  
  head(sample_cortex) # checks
  
  
  #######'Brain - Putamen (basal ganglia)'######################
  
  # all putamen samples 
  
  putamen <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Putamen (basal ganglia)') %>%
    dplyr::select("SAMPID") 
  
  dim(putamen) # 205 
  
  # random 25 sample vector  
  
  sampleindex_putamen <- sample(putamen$SAMPID, 25)
  
  sample_putamen <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_putamen))
  
  head(sample_putamen) # checks
  
  
  #######'Brain - Nucleus accumbens (basal ganglia)'######################

  # all putamen samples 
  
  n_accubens <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Nucleus accumbens (basal ganglia)') %>%
    dplyr::select("SAMPID") 
  
  dim(n_accubens) # 246 
  
  # a vector randomly sampled 20   
  
  sampleindex_naccubens <- sample(n_accubens$SAMPID, 25)
  
  sample_naccubens <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_naccubens))
  
  head(sample_naccubens) # checks
  
  
  #######'Brain - Caudate (basal ganglia)'######################
  
  caudate <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Caudate (basal ganglia)') %>%
    dplyr::select("SAMPID") 
  
  dim(caudate)  #246
  
  # random 25 sample vector  
  
  sampleindex_caudate <- sample(caudate$SAMPID, 25)
  
  sample_caudate <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_caudate))
  
  head(sample_caudate) # checks

  
  #######'Brain - Cerebellum'######################
  
  cerebellum <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Cerebellum') %>%
    dplyr::select("SAMPID") 
  
  dim(cerebellum)  #241
  
  # random 25 sample vector  
  
  sampleindex_cerebellum <- sample(cerebellum$SAMPID, 25)
  
  sample_cerebellum <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_cerebellum))
  
  head(sample_cerebellum) # checks
  
 
  #######'Brain - Cerebellar Hemisphere'######################
  
  cerebellum_hemi <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Cerebellar Hemisphere') %>%
    dplyr::select("SAMPID") 
  
  dim(cerebellum_hemi)  #215
  
  # random 25 sample vector  
  
  sampleindex_cerebellum_hemi<- sample(cerebellum_hemi$SAMPID, 25)
  
  sample_cerebellum_hemi <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_cerebellum_hemi))
  
  head(sample_cerebellum_hemi) # checks
  
  
  #######'Brain - Anterior cingulate cortex (BA24)'######################
  
  # all accortex samples 
  
  ac_cortex <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Anterior cingulate cortex (BA24)') %>%
    dplyr::select("SAMPID") 
  
  dim(ac_cortex) # 176 
  
  # random 25 sample vector  
  
  sampleindex_ac_cortex <- sample(ac_cortex$SAMPID, 25)
  
  sample_ac_cortex <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_ac_cortex))
  
  head(sample_ac_cortex) # checks
  
  
  #######'Brain - Frontal Cortex (BA9)'######################
  
  fro_cortex <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Frontal Cortex (BA9)') %>%
    dplyr::select("SAMPID") 
  
  dim(fro_cortex) # 209 
  
  # random 25 sample vector  
  
  sampleindex_fro_cortex <- sample(fro_cortex$SAMPID, 25)
  
  sample_fro_cortex <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_fro_cortex))
  
  head(sample_fro_cortex) # checks
  
  
  #######'Brain - Hypothalamus'######################
  
  hypothalamus <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Hypothalamus') %>%
    dplyr::select("SAMPID") 
  
  dim(hypothalamus) # 202
  
  # random 25 sample vector  
  
  sampleindex_hypothalamus <- sample(hypothalamus$SAMPID, 25)
  
  sample_hypothalamus <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_hypothalamus))
  
  head(sample_hypothalamus) # checks
  
  
  #######'Brain - Hippocampus'######################
  
  hippocampus <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Hippocampus') %>%
    dplyr::select("SAMPID") 
  
  dim(hippocampus) # 197
  
  # random 25 sample vector  
  
  sampleindex_hippocampus <- sample(hippocampus$SAMPID, 25)
  
  sample_hippocampus <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_hippocampus))
  
  head(sample_hippocampus) # checks
  
  
  #######'Brain - Amygdala'######################
  
  amygdala <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Amygdala') %>%
    dplyr::select("SAMPID") 
  
  dim(amygdala) # 152
  
  # random 25 sample vector  
  
  sampleindex_amygdala <- sample(amygdala$SAMPID, 25)
  
  sample_amygdala <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_amygdala))
  
  head(sample_amygdala) # checks
  
  
  #######'Brain - Substantia nigra'######################
  
  substantia_nigra <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Substantia nigra') %>%
    dplyr::select("SAMPID") 
  
  dim(substantia_nigra) # 139
  
  # random 25 sample vector  
  
  sampleindex_substantia_nigra <- sample(substantia_nigra$SAMPID, 25)
  
  sample_substantia_nigra <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_substantia_nigra))
  
  head(sample_substantia_nigra) # checks
  
  
  #######'Brain - Spinal cord (cervical c-1)'######################
  
  spinal_cord <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Brain - Spinal cord (cervical c-1)') %>%
    dplyr::select("SAMPID") 
  
  dim(spinal_cord) # 159
  
  # random 25 sample vector  
  
  sampleindex_spinal_cord <- sample(spinal_cord$SAMPID, 25)
  
  sample_spinal_cord <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_spinal_cord))
  
  head(sample_spinal_cord) # checks
  
  #######'Nerve - Tibial'######################
  
  nerve_tibia <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Nerve - Tibial') %>%
    dplyr::select("SAMPID") 
  
  dim(nerve_tibia) # 619
  
  # random 25 sample vector  
  
  sampleindex_nerve_tibia <- sample(nerve_tibia$SAMPID, 25)
  
  sample_nerve_tibia <- cns_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_nerve_tibia))
  
  head(sample_nerve_tibia) # checks
  
  
  #######'Small Intestine - Terminal Ileum' ######################
  
  ileum <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Small Intestine - Terminal Ileum') %>%
    dplyr::select("SAMPID") 
  
  dim(ileum) # 187
  
  # random 25 sample vector  
  
  sampleindex_ileum <- sample(ileum$SAMPID, 25)
  
  sample_ileum <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_ileum))
  
  head(sample_ileum) # checks
  
  
  #######'Colon - Sigmoid'######################
  
  colon_sigmoid <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Colon - Sigmoid') %>%
    dplyr::select("SAMPID") 
  
  dim(colon_sigmoid) # 373
  
  # random 25 sample vector  
  
  sampleindex_colon_sigmoid <- sample(colon_sigmoid$SAMPID, 25)
  
  sample_colon_sigmoid <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_colon_sigmoid))
  
  head(sample_colon_sigmoid) # checks
  
  
  #######'Colon - Transverse'######################
  
  colon_transverse <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Colon - Transverse') %>%
    dplyr::select("SAMPID") 
  
  dim(colon_transverse) # 406
  
  # random 25 sample vector  
  
  sampleindex_colon_transverse <- sample(colon_transverse$SAMPID, 25)
  
  sample_colon_transverse <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_colon_transverse))
  
  head(sample_colon_transverse) # checks
  
  
  #######'Esophagus - Muscularis' ######################
  
  esophagus <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Esophagus - Muscularis') %>%
    dplyr::select("SAMPID") 
  
  dim(esophagus) # 515
  
  # random 25 sample vector  
  
  sampleindex_esophagus <- sample(esophagus$SAMPID, 25)
  
  sample_esophagus <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_esophagus))
  
  head(sample_esophagus) # checks
  
  #######'Esophagus - GIs' ######################
  
  gastroesophagus <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Esophagus - Gastroesophageal Junction') %>%
    dplyr::select("SAMPID") 
  
  dim(gastroesophagus) # 375
  
  # random 25 sample vector  
  
  sampleindex_gastroesophagus <- sample(gastroesophagus$SAMPID, 25)
  
  sample_gastroesophagus <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_gastroesophagus))
  
  head(sample_gastroesophagus) # checks
  
  
  #######'Esophagus - Mucosa' ######################
  
  mucosa <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Esophagus - Mucosa') %>%
    dplyr::select("SAMPID") 
  
  dim(mucosa) # 555
  
  # random 25 sample vector  
  
  sampleindex_mucosa <- sample(mucosa$SAMPID, 25)
  
  sample_mucosa <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_mucosa))
  
  head(sample_mucosa) # checks
  
  
  #######'Stomach'######################
  
  stomach <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Stomach') %>%
    dplyr::select("SAMPID") 
  
  dim(stomach) # 359
  
  # random 25 sample vector  
  
  sampleindex_stomach <- sample(stomach$SAMPID, 25)
  
  sample_stomach <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_stomach))
  
  head(sample_stomach) # checks
  
  #######'Spleen'######################
  
  spleen <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Spleen') %>%
    dplyr::select("SAMPID") 
  
  dim(spleen) # 241
  
  # random 25 sample vector  
  
  sampleindex_spleen <- sample(spleen$SAMPID, 25)
  
  sample_spleen <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_spleen))
  
  head(sample_spleen) # checks
  
  #######'Pancreas'######################
  
  pancreas <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Pancreas') %>%
    dplyr::select("SAMPID") 
  
  dim(pancreas) # 328
  
  # random 25 sample vector  
  
  sampleindex_pancreas <- sample(pancreas$SAMPID, 25)
  
  sample_pancreas <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_pancreas))
  
  head(sample_pancreas) # checks
  
  ####### 'Liver'######################
  
  liver <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Liver') %>%
    dplyr::select("SAMPID") 
  
  dim(liver) # 226
  
  # random 25 sample vector  
  
  sampleindex_liver <- sample(liver$SAMPID, 25)
  
  sample_liver <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_liver))
  
  head(sample_liver) # checks
  
  
  ####### 'Salivary Gland'######################
  
  sgland <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Minor Salivary Gland') %>%
    dplyr::select("SAMPID") 
  
  dim(sgland) # 162
  
  # random 25 sample vector  
  
  sampleindex_sgland <- sample(sgland$SAMPID, 25)
  
  sample_sgland <- ens_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_sgland))
  
  head(sample_sgland) # checks
  
  
  #######'Artery - Tibial'######################
  
  artery_tibia <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Artery - Tibial') %>%
    dplyr::select("SAMPID") 
  
  dim(artery_tibia) # 663
  
  # random 25 sample vector  
  
  sampleindex_artery_tibia <- sample(artery_tibia$SAMPID, 25)
  
  sample_artery_tibia <- cv_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_artery_tibia))
  
  head(sample_artery_tibia) # checks
  
  
  #######'Artery - Coronary'######################
  
  artery_coronary <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Artery - Coronary') %>%
    dplyr::select("SAMPID") 
  
  dim(artery_coronary) # 240
  
  # random 25 sample vector  
  
  sampleindex_artery_coronary <- sample(artery_coronary$SAMPID, 25)
  
  sample_artery_coronary <- cv_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_artery_coronary))
  
  head(sample_artery_coronary) # checks
  
  
  ####### 'Artery - Aorta'######################
  
  artery_aorta <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Artery - Aorta') %>%
    dplyr::select("SAMPID") 
  
  dim(artery_aorta) # 432
  
  # random 25 sample vector  
  
  sampleindex_artery_aorta <- sample(artery_aorta$SAMPID, 25)
  
  sample_artery_aorta <- cv_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_artery_aorta))
  
  head(sample_artery_aorta) # checks
  
  
  #######'Heart - Atrial Appendage' ######################
  
  heart_atria <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Heart - Atrial Appendage') %>%
    dplyr::select("SAMPID") 
  
  dim(heart_atria) # 429
  
  # random 25 sample vector  
  
  sampleindex_heart_atria <- sample(heart_atria$SAMPID, 25)
  
  sample_heart_atria <- cv_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_heart_atria))
  
  head(sample_heart_atria) # checks
  
  
  #######'Heart - Left Ventricle'######################
  
  heart_lv <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Heart - Left Ventricle') %>%
    dplyr::select("SAMPID") 
  
  dim(heart_lv) # 432
  
  # random 25 sample vector  
  
  sampleindex_heart_lv <- sample(heart_lv$SAMPID, 25)
  
  sample_heart_lv <- cv_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_heart_lv))
  
  head(sample_heart_lv) # checks
  
  
  ####### 'Pituitary'######################
  
  pituitary <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Pituitary') %>%
    dplyr::select("SAMPID") 
  
  dim(pituitary) # 283
  
  # random 25 sample vector  
  
  sampleindex_pituitary <- sample(pituitary$SAMPID, 25)
  
  sample_pituitary <- endo_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_pituitary))
  
  head(sample_pituitary) # checks
  
  
  ####### 'Thyroid'######################

  thyroid <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Thyroid') %>%
    dplyr::select("SAMPID") 
  
  dim(thyroid) # 653
  
  # random 25 sample vector  
  
  sampleindex_thyroid <- sample(thyroid$SAMPID, 25)
  
  sample_thyroid <- endo_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_thyroid))
  
  head(sample_thyroid) # checks
  
  
  ####### 'Adrenal Gland'######################
  
  adrenal <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Adrenal Gland') %>%
    dplyr::select("SAMPID") 
  
  dim(adrenal) # 258
  
  # random 25 sample vector  
  
  sampleindex_adrenal <- sample(adrenal$SAMPID, 25)
  
  sample_adrenal <- endo_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_adrenal))
  
  head(sample_adrenal) # checks
  
  ####### 'Uterus'######################
  
  uterus <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Uterus') %>%
    dplyr::select("SAMPID") 
  
  dim(uterus) # 142
  
  # random 25 sample vector  
  
  sampleindex_uterus <- sample(uterus$SAMPID, 25)
  
  sample_uterus <- femm_repro_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_uterus))
  
  head(sample_uterus) # checks
  
  
  ####### 'Ovary'######################

  ovary <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Ovary') %>%
    dplyr::select("SAMPID") 
  
  dim(ovary) # 180
  
  # random 25 sample vector  
  
  sampleindex_ovary <- sample(ovary$SAMPID, 25)
  
  sample_ovary <- femm_repro_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_ovary))
  
  head(sample_ovary) # checks
  
  ####### 'Vagina'######################
  
  vagina <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Vagina') %>%
    dplyr::select("SAMPID") 
  
  dim(vagina) # 156
  
  # random 25 sample vector  
  
  sampleindex_vagina <- sample(vagina$SAMPID, 25)
  
  sample_vagina <- femm_repro_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_vagina))
  
  head(sample_vagina) # checks
  
  ####### 'Cervix - Endocervix'######################
  
  endocervix <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Cervix - Endocervix') %>%
    dplyr::select("SAMPID") 
  
  dim(endocervix) # 10
  
  sampleindex_endocervix <- sample(endocervix$SAMPID, 10)
  
  sample_endocervix <- femm_repro_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_endocervix))
  
  head(sample_endocervix) # checks
  
  
  ####### 'Cervix - Ectocervix'######################
  
  ectocervix <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Cervix - Ectocervix') %>%
    dplyr::select("SAMPID") 
  
  dim(ectocervix) # 9
  
  sampleindex_ectocervix <- sample(ectocervix$SAMPID, 9)
  
  sample_ectocervix <- femm_repro_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_ectocervix))
  
  head(sample_ectocervix) # checks

  ####### 'Fallopian Tube'######################
  
  ftube <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Fallopian Tube') %>%
    dplyr::select("SAMPID") 
  
  dim(ftube) # 9
  
  sampleindex_ftube <- sample(ftube$SAMPID, 9)
  
  sample_ftube <- femm_repro_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_ftube))
  
  head(sample_ftube) # checks
  
  
  ####### 'Testes'######################
  
  testes <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Testis') %>%
    dplyr::select("SAMPID") 
  
  dim(testes) # 361
  
  # random 25 sample vector  
  
  sampleindex_testes <- sample(testes$SAMPID, 25)
  
  sample_testes <- mus_repro_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_testes))
  
  head(sample_testes) # checks
  
  ####### 'protrates'######################
  
  prostate <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Prostate') %>%
    dplyr::select("SAMPID") 
  
  dim(prostate) # 245
  
  # random 25 sample vector  
  
  sampleindex_prostate <- sample(prostate$SAMPID, 25)
  
  sample_prostate <- mus_repro_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_prostate))
  
  head(sample_testes) # checks
  
  ####### 'Bladder'######################
  
  bladder <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Bladder') %>%
    dplyr::select("SAMPID") 
  
  dim(bladder) # 21
  
  # random 25 sample vector  
  
  sampleindex_bladder <- sample(bladder$SAMPID, 21)
  
  sample_bladder <- excre_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_bladder))
  
  head(sample_bladder) # checks
  
  ####### 'Kidney - Cortex' #############
  
  kidney_cortex <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Kidney - Cortex') %>%
    dplyr::select("SAMPID") 
  
  dim(kidney_cortex) # 25
  
  # random 25 sample vector  
  
  sampleindex_kidney_cortex <- sample(kidney_cortex$SAMPID, 25)
  
  sample_kidney_cortex <- excre_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_kidney_cortex))
  
  head(sample_kidney_cortex) # checks
  
  ####### 'Kidney - Medulla' #############
  
  kidney_medulla <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Kidney - Medulla') %>%
    dplyr::select("SAMPID") 
  
  dim(kidney_medulla) # 4
  
  sampleindex_kidney_medulla <- sample(kidney_medulla$SAMPID, 4)
  
  sample_kidney_medulla <- excre_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_kidney_medulla))
  
  head(sample_kidney_medulla) # checks
  
  
  ########'Adipose - Subcutaneous'################

  adipose_sub <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Adipose - Subcutaneous') %>%
    dplyr::select("SAMPID") 
  
  dim(adipose_sub) # 663
  
  # random 25 sample vector  
  
  sampleindex_adipose_sub <- sample(adipose_sub$SAMPID, 25)
  
  sample_adipose_sub <- adi_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_adipose_sub))
  
  head(sample_adipose_sub) # checks
  
  ########'Breast - Mammary Tissue'################ 
  
 breast <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Breast - Mammary Tissue') %>%
    dplyr::select("SAMPID") 
  
  dim(breast) # 459
  
  # random 25 sample vector  
  
  sampleindex_breast <- sample(breast$SAMPID, 25)
  
  sample_breast <- adi_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_breast))
  
  head(sample_breast) # checks 
  
  ########'Adipose - Visceral (Omentum)'################

  adipose_vis <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Adipose - Visceral (Omentum)') %>%
    dplyr::select("SAMPID") 
  
  dim(adipose_vis) # 541
  
  # random 25 sample vector  
  
  sampleindex_adipose_vis <- sample(adipose_vis$SAMPID, 25)
  
  sample_adipose_vis <- adi_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_adipose_vis))
  
  head(sample_adipose_vis) # checks  
  
  ########'muscle'################
  
  muscle <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Muscle - Skeletal') %>%
    dplyr::select("SAMPID") 
  
  dim(muscle) # 803
  
  # random 25 sample vector  
  
  sampleindex_muscle <- sample(muscle$SAMPID, 25)
  
  sample_muscle <- muscle_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_muscle))
  
  head(sample_muscle) # checks  
  
  ########'No Sun Skin'################
  
  skin_no_sun <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Skin - Not Sun Exposed (Suprapubic)') %>%
    dplyr::select("SAMPID") 
  
  dim(skin_no_sun) # 604
  
  # random 25 sample vector  
  
  sampleindex_skin_no_sun <- sample(skin_no_sun$SAMPID, 25)
  
  sample_skin_no_sun <- skin_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_skin_no_sun))
  
  head(sample_skin_no_sun) # checks  
  
  ########'Sun skin'################
  
  skin_sun <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Skin - Sun Exposed (Lower leg)') %>%
    dplyr::select("SAMPID") 
  
  dim(skin_sun) # 1701
  
  # random 25 sample vector  
  
  sampleindex_skin_sun <- sample(skin_sun$SAMPID, 25)
  
  sample_skin_sun <- skin_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_skin_sun))
  
  head(sample_skin_sun) # checks  
  
  ########'Lung'################
  
  lung <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Lung') %>%
    dplyr::select("SAMPID") 
  
  dim(lung) # 578
  
  # random 25 sample vector  
  
  sampleindex_lung <- sample(lung$SAMPID, 25)
  
  sample_lung <- lung_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_lung))
  
  head(sample_lung) # checks  
  
  ########'Lymphocyte'################
  
  lcyte <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Cells - EBV-transformed lymphocytes') %>%
    dplyr::select("SAMPID") 
  
  dim(lcyte) # 174
  
  # random 25 sample vector  
  
  sampleindex_lcyte <- sample(lcyte$SAMPID, 25)
  
  sample_lcyte <- lcyte_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_lcyte))
  
  head(sample_lcyte) # checks  
  
  ########'Whole blood'################
  
  blood <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Whole Blood') %>%
    dplyr::select("SAMPID") 
  
  dim(blood) # 755
  
  # random 25 sample vector  
  
  sampleindex_blood <- sample(blood$SAMPID, 25)
  
  sample_blood <- blood_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_blood))
  
  head(sample_blood) # checks  
  
  ########'FBLAST'################
  
  fblast <- meta %>%
    filter(SMAFRZE == "RNASEQ") %>%
    filter(SMTSD %in% 'Cells - Cultured fibroblasts') %>%
    dplyr::select("SAMPID") 
  
  dim(fblast) # 504
  
  # random 25 sample vector  
  
  sampleindex_fblast <- sample(fblast$SAMPID, 25)
  
  sample_fblast <- fblast_data %>% 
    dplyr::select(c("Name","Description",  sampleindex_fblast))
  
  head(sample_fblast) # checks  
  
############################################################################### 

# joining all  

############################################################################### 

  
  # new data frame 25 columns from each tissue  with a few exceptions 


  bigGTEX <- 
    sample_cortex %>% 
    full_join(sample_putamen, by = c("Name", "Description")) %>% 
    full_join(sample_naccubens, by = c("Name", "Description"))  %>%
    full_join(sample_caudate, by = c("Name", "Description")) %>%
    full_join(sample_cerebellum, by = c("Name", "Description")) %>%
    full_join(sample_cerebellum_hemi, by = c("Name", "Description")) %>%
    full_join(sample_ac_cortex, by = c("Name", "Description")) %>%
    full_join(sample_fro_cortex, by = c("Name", "Description")) %>%
    full_join(sample_hypothalamus, by = c("Name", "Description")) %>%
    full_join(sample_hippocampus, by = c("Name", "Description")) %>%
    full_join(sample_amygdala, by = c("Name", "Description")) %>%
    full_join(sample_substantia_nigra, by = c("Name", "Description")) %>%
    full_join(sample_spinal_cord, by = c("Name", "Description")) %>%
    full_join(sample_nerve_tibia, by = c("Name", "Description")) %>%
    full_join(sample_ileum, by = c("Name", "Description")) %>%
    full_join(sample_colon_sigmoid, by = c("Name", "Description")) %>%
    full_join(sample_colon_transverse, by = c("Name", "Description")) %>%
    full_join(sample_esophagus, by = c("Name", "Description")) %>%
    full_join(sample_mucosa, by = c("Name", "Description")) %>%
    full_join(sample_gastroesophagus, by = c("Name", "Description")) %>%
    full_join(sample_stomach, by = c("Name", "Description")) %>%
    full_join(sample_pancreas, by = c("Name", "Description")) %>%
    full_join(sample_spleen, by = c("Name", "Description")) %>%
    full_join(sample_sgland, by = c("Name", "Description")) %>%
    full_join(sample_liver, by = c("Name", "Description"))  %>%
    full_join(sample_artery_tibia, by = c("Name", "Description")) %>%
    full_join(sample_artery_coronary, by = c("Name", "Description")) %>%
    full_join(sample_artery_aorta, by = c("Name", "Description")) %>%
    full_join(sample_heart_atria, by = c("Name", "Description")) %>%
    full_join(sample_heart_lv, by = c("Name", "Description")) %>%
    full_join(sample_lung, by = c("Name", "Description")) %>%
    full_join(sample_pituitary, by = c("Name", "Description")) %>%
    full_join(sample_thyroid, by = c("Name", "Description")) %>%
    full_join(sample_adrenal, by = c("Name", "Description")) %>%
    full_join(sample_uterus, by = c("Name", "Description")) %>%
    full_join(sample_ovary, by = c("Name", "Description")) %>%
    full_join(sample_ftube, by = c("Name", "Description"))  %>%
    full_join(sample_vagina, by = c("Name", "Description")) %>%
    full_join(sample_endocervix, by = c("Name", "Description")) %>%
    full_join(sample_ectocervix, by = c("Name", "Description")) %>%
    full_join(sample_testes, by = c("Name", "Description")) %>%
    full_join(sample_prostate, by = c("Name", "Description")) %>%
    full_join(sample_bladder, by = c("Name", "Description")) %>%
    full_join(sample_kidney_cortex, by = c("Name", "Description")) %>%
    full_join(sample_kidney_medulla, by = c("Name", "Description")) %>%
    full_join(sample_breast, by = c("Name", "Description")) %>%
    full_join(sample_adipose_sub, by = c("Name", "Description")) %>%
    full_join(sample_adipose_vis, by = c("Name", "Description"))  %>%
    full_join(sample_skin_no_sun, by = c("Name", "Description"))  %>%
    full_join(sample_skin_sun, by = c("Name", "Description"))  %>%
    full_join(sample_blood, by = c("Name", "Description"))  %>%
    full_join(sample_lcyte, by = c("Name", "Description"))  %>%
    full_join(sample_fblast, by = c("Name", "Description"))  %>%
    full_join(sample_muscle, by = c("Name", "Description"))  

 dim(bigGTEX) # 1280 from 17382 
 
  write.csv(bigGTEX, "Data/bigGTEX.csv", sep = "\t") 
   
#############################################################################
  
  # tissue information 

##############################################################################    
  
  bigGTEX_demo <- read.csv("Data/bigGTEX.csv",
                             check.names = F,
                             header = T,
                             row.names = 1)
  
  dim(bigGTEX_demo) # 1280
  
  bigGTEX_sample_names <- colnames(bigGTEX_demo)

  bigGTEX_sample_names <- bigGTEX_sample_names[-c (1,2)] 
  
  length(bigGTEX_sample_names) # 1278
  
  class(bigGTEX_sample_names) # char
  
  bigGTEX_sample_tissues <- meta %>% 
    filter(SAMPID %in% bigGTEX_sample_names) %>%
    dplyr::select(c("SAMPID", "SMTSD"))
  
  dim(bigGTEX_sample_tissues) # 1278
  
  # 
  # misssingsamples <- bigGTEX_sample_names [! (bigGTEX_sample_names %in% bigGTEX_sample_tissues$SAMPID)]
  # 
  # misssingsamples <- as.data.frame(misssingsamples)
  # 
  # try <- c("GTEX-1EKGG-1226-SM-7IGNO",
  #          "GTEX-12WSI-0226-SM-5GCNA",
  #          "GTEX-147GR-1326-SM-7IGLB",
  #          "GTEX-R53T-0326-SM-48FEC" ,
  #          "GTEX-SJXC-1226-SM-4DM78" ,
  #           "GTEX-1IKJJ-1526-SM-C1YSC",
  #           "GTEX-14XAO-0226-SM-68728",
  #           "GTEX-XXEK-1126-SM-4BRUX")
  # 
  # meta %>% 
  #   filter(SAMPID %in% try) %>%
  #   dplyr::select(c("SAMPID", "SMTSD"))
  # 
  
  write.csv(bigGTEX_sample_tissues, "Data/bigGTEXTissues.csv", sep = "\t") 
  

  ###################################
  
  mean_tissue <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", 
                       stringsAsFactors = F,
                       header = T, 
                       sep = "\t", 
                       quote = "",
                       skip = 2)  
  
  mean_tissue <- mean_tissue [, -c(1,2)]
  
  class(mean_tissue$`Adipose - Subcutaneous`)
  
  row_sub <- apply(mean_tissue, 1, function(row) all(row !=0 ))
    
  mean_tissue <- mean_tissue[row_sub, ]
  
  my_median <- data.frame(apply(mean_tissue,2,median))
  
  my_median$Tissue <- rownames(my_median)
  
  colnames(my_median)[which(names(my_median) == "apply.mean_tissue..2..median.")] <- "Median" 
  
  ggplot(my_median, aes(Tissue, Median))+
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 8)) 
  
  median(my_median$Median)
  
  