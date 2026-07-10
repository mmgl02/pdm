# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)

# Links
txt_file <- "data_raw/envidatS3paths.txt"

# Where to save
dest_dir <- path_raw_bioclim
dir.create(dest_dir, showWarnings = FALSE)

# Read links and suppress space
links <- trimws(readLines(txt_file))

# Download
for (link in links) {

  file_name <- basename(link)
  
  dest_file <- file.path(dest_dir, file_name)
  
  download.file(url = link, destfile = dest_file, mode = "wb")
  
  cat("Downloaded :", file_name, "\n")
}