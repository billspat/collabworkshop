## TITLE:         Workshop main file
## AUTHORS:       Pat Bills

## DATA:          NEON organismal data pull (https://drive.google.com/drive/folders/10vjvS6E3RoJfmP-w8nCVDKMPqP20yvfC)
## PROJECT:       Macrosystems biodiversity project Workshop  
## DATE:          initiated: Nov 3, 2020

library(readr)
library(dplyr)

source('google_drive.R')

# very very basic function 
count_by_taxon <- function(file_path){
    data <- read.csv(file_path)
    count(data, taxonID)
}

download_path = 'downloads'
download_neon_taxonomy_lists()


# goals for this script
# 1. download CSV data in google_drive.R script
# 2. very basic summary functions for each file downloaded, each person creates a summary function
#    for example: 

#    fish_trapping_summary <- function(){ 
#       fish_trapping <- count_by_taxon(file.path(download_path,"fish_trapping.csv"))
#       return(fish_trapping)
#    }
