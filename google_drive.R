## TITLE:         Google Drive download utility
## AUTHOR:        Pat Bills billspat@msu.edu
## DATA:          None
## PROJECT:       Macrosystems biodiversity project Workshop : 
## DATE:          initiated: Nov 3, 2020

library(googledrive)
# downloads CVS from a google drive folder by folder ID in to a data folder (downloads by defaul)
# for working with small size csv files 
download_csv_files <-function(data_folder_id, download_path) {
    data_folder = drive_get(id=as_id(data_folder_id))
    all_files = drive_ls(data_folder)
    csv_files <- all_files[grep(".csv$", all_files$name),]
    
    if (!dir.exists( download_path )){ dir.create( download_path )}
    for(row in 1:nrow(csv_files)) {
        
        drive_download( csv_files[row,],  file.path( download_path, csv_files[row,]$name), overwrite = TRUE  ) 
    }
}

# example function with default for NEON csv data
# obtain this folder ID by browsing google drive folders, and copying the last part of the url
# for this example,  "https://drive.google.com/drive/folders/10vjvS6E3RoJfmP-w8nCVDKMPqP20yvfC"

# MacrosystemsBiodiversity/data/NEON_Biodiversity/L0/organismal_data_latest_pull
download_neon_taxonomy_lists<-function(data_folder_id = "10vjvS6E3RoJfmP-w8nCVDKMPqP20yvfC", download_path = "downloads") {
    download_csv_files(data_folder_id, download_path)
}

