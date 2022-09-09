# HTR_Transcriptome

## Data 

- https://human.brain-map.org/static/download
- https://gtexportal.org/home/tissueSummaryPage

## Analysis Based on

- https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html for heatmaps
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/ for data pre-processing and DEG
- https://f1000research.com/articles/9-1444 for design matrix 
- http://genomicsclass.github.io/book/ for statistics 

## File name: GTEX_subsetting

**Code function:**

splitting the GTEX data into system-specific data sets

**Location in Results:**

Classifying the HTRs by Shared Pattern of Expression (Thesis Chapter Two: Materials and Methods)

## File name: GTEX_propersampling

**Code function:**

making a small data set with 25 randomly selected sample from each tissue 

**Location in Results:**

Classifying the HTRs by Shared Pattern of Expression (Thesis Chapter Two: Materials and Methods)

## File name: GTEX_filtering

**Code function:**

setting the threshold for detection  

**Location in Results:**

Classifying the HTRs by Shared Pattern of Expression (Thesis Chapter Two: Materials and Methods)

## File name: GTEX_detection_data

**Code function:**

creating the probabilityXtissue matrix for HTR detection percentage  

**Location in Results:**

HTRs Can be Grouped by Tissue Distribution (Thesis Chapter Four: The Specialised Expression Patterns of HTRs: A Family of Diverse Roles in Tissues)

## File name: GTEX_detection

**Code function:**

making heatmaps using the percentage of reproducible detection calculated before  

**Location in Results:**

HTRs Can be Grouped by Tissue Distribution (Thesis Chapter Four: The Specialised Expression Patterns of HTRs: A Family of Diverse Roles in Tissues)

## File name: GTEX_htr_df_build

**Code function:**

shorten the dataframe just to have the HTR genes  

**Location in Results:**

HTRs Can be Grouped by Tissue Distribution (Thesis Chapter Four: The Specialised Expression Patterns of HTRs: A Family of Diverse Roles in Tissues)

## File name: AllenBrain_filtering

**Code function:**

evaluation the threshold of expression detection 

**Location in Results:**

Classifying the HTRs in Brain by Relative Abundance of Expression (Thesis Chapter Two: Materials and Methods)

## File name: AllenBrain_compare_donor

**Code function:**

compare the two donors to see if they can be treated as replicates 

**Location in Results:**

Classifying the HTRs in Brain by Relative Abundance of Expression (Thesis Chapter Two: Materials and Methods)

## File name: AllenBrain_qc

**Code function:**

go through quality-control checkpoints

**Location in Results:**

Classifying the HTRs in Brain by Relative Abundance of Expression (Thesis Chapter Two: Materials and Methods)

## File name: AllenBrain_abundance

**Code function:**

look at the relative abundance of the HTR in brain tissues

**Location in Results:**

Five Patterns of HTR Expression Exists in the Brain (Thesis Chapter Four: The Specialised Expression Patterns of HTRs: A Family of Diverse Roles in Tissues)

## File name: AllenBrain_lm

**Code function:**

fit the data to the linear model, look for variability across tissues, account for donor variability 

**Location in Results:**

Five Patterns of HTR Expression Exists in the Brain (Thesis Chapter Four: The Specialised Expression Patterns of HTRs: A Family of Diverse Roles in Tissues)

