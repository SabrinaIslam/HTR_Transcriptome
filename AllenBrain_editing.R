################################################################################

# packages 

################################################################################

# RNA-seq specific packages  

library(dplyr)
library(tidyverse)

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

# editing the DEG output 

# changing column order

# adding proper column names 

# rounding the digits 

# adding "*" as significance level indicator

# removing duplicates of pairwise t test 

################################################################################

# edited the .csv s in shel 

# combine.sh 

# used grep in bash to grep receptors each into 11 receptor.csv

# concat.sh 

# cat an output.csv with all receptors

# edited the outout.csv to drop column 1 and 8 with cut -d -f --complement 

# saved the new files as results.csv 

################################################################################

raw_output <- read.csv("limma_output/2.bash_output/results.csv",
                       header = F) # results 

htr_abundance_brain <- raw_output[, c(1,7,2,3,4,5,6)] # changing the order 

names(htr_abundance_brain) <- c("Receptor",
                                "Comparison",
                                "log2 fold-change",
                                "Average Expression",
                                "t-statistic",
                                "P-value",
                                "Adjusted P-value") # colnames 
# rounding the digits 

htr_abundance_brain$`log2 fold-change` <- round(htr_abundance_brain$`log2 fold-change`, digits = 3)

htr_abundance_brain$`Average Expression` <- round(htr_abundance_brain$`Average Expression`, digits = 3)

# adding * 

htr_abundance_brain$Significance <- ifelse(htr_abundance_brain$`Adjusted P-value` < 0.001, "***",
                                                ifelse(htr_abundance_brain$`Adjusted P-value` < 0.01, "**",
                                                       ifelse(htr_abundance_brain$`Adjusted P-value` < 0.05, "*",
                                                              ".")))

# library(kableExtra)
# lx <- htr_abundance_brain_edit %>% kable("latex") %>% head()

remove <- c("CgGVsCbCx",
            "FLVsCbCx",
            "FLVsCgG",
            "GPVsCbCx",
            "GPVsCgG",
            "GPVsFL",
            "InsVsCbCx",
            "InsVsCgG",
            "InsVsFL",
            "InsVsGP",
            "OLVsCbCx",
            "OLVsCgG",
            "OLVsFL",
            "OLVsGP",
            "OLVsIns",
            "PHGVsCbCx",
            "PHGVsCgG",
            "PHGVsFL",
            "PHGVsGP",
            "PHGVsIns",
            "PHGVsOL",
            "PLVsCbCx",
            "PLVsCgG",
            "PLVsFL",
            "PLVsGP",
            "PLVsIns",
            "PLVsOL",
            "PLVsPHG",
            "StrVsCbCx",
            "StrVsCgG",
            "StrVsFL",
            "StrVsGP",
            "StrVsIns",
            "StrVsOL",
            "StrVsPHG",
            "StrVsPL",
            "TLVsCbCx",
            "TLVsCgG",
            "TLVsFL",
            "TLVsGP",
            "TLVsIns",
            "TLVsOL",
            "TLVsPHG",
            "TLVsPL",
            "TLVsStr") # are the duplicate pairwise t test 

length(remove)

results_final <- htr_abundance_brain[-c(which(results_2$Comparison %in% remove)), ] # removing duplicates 

rownames(results_final) <- c(1:450) # renaming rows 
  
head(results_final) # all good 

write.csv(results_final, "result_final.csv")
