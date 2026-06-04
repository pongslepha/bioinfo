# GEO Raw and Supplementary Data Downloader
# Saves raw data (.tar) or individual supplementary files of GSE39872, GSE134033, and GSE37114 to 00.data/GEO

# 1. Define Accession list
gse_ids <- c("GSE39872", "GSE134033", "GSE37114")

# 2. Output directory
# Determine script directory dynamically to make it robust
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(sub("--file=", "", file_arg))
  } else {
    # Fallback to default if run interactively in R console
    return("xx.script/01.download_geo.R")
  }
}

script_dir <- dirname(normalizePath(get_script_path(), mustWork = FALSE))
dest_dir <- file.path(script_dir, "../00.data/GEO")

if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
  message("Created output directory: ", dest_dir)
}

# 3. Download function with automatic fallback
download_geo_data <- function(gse_id, dest_dir) {
  # Extract digits to construct the GEO FTP/HTTP URL structure
  digits <- gsub("GSE", "", gse_id)
  num_digits <- nchar(digits)
  
  if (num_digits <= 3) {
    subfolder <- "GSEnnn"
  } else {
    subfolder <- paste0("GSE", substr(digits, 1, num_digits - 3), "nnn")
  }
  
  # Base URLs
  suppl_dir_url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", subfolder, gse_id)
  raw_tar_url <- sprintf("%s%s_RAW.tar", suppl_dir_url, gse_id)
  raw_tar_dest <- file.path(dest_dir, paste0(gse_id, "_RAW.tar"))
  
  message("\n==========================================")
  message(sprintf("Processing Accession: %s", gse_id))
  message("==========================================")
  
  # Step A: Attempt to download primary _RAW.tar file
  if (file.exists(raw_tar_dest) && file.info(raw_tar_dest)$size > 1000) {
    message(sprintf("Primary RAW tar for %s already exists in %s. Skipping download.", gse_id, dest_dir))
    download_success <- TRUE
  } else {
    message(sprintf("Attempting to download primary RAW tar from: %s", raw_tar_url))
    
    download_success <- FALSE
    tryCatch({
      # Check headers or try downloading directly (libcurl handles HTTP errors nicely)
      download.file(raw_tar_url, destfile = raw_tar_dest, method = "libcurl", mode = "wb", quiet = TRUE)
      
      # Check if download succeeded and file is non-empty
      if (file.exists(raw_tar_dest) && file.info(raw_tar_dest)$size > 1000) {
        message(sprintf("Successfully downloaded primary RAW tar: %s", basename(raw_tar_dest)))
        download_success <- TRUE
      } else {
        # Remove invalid/small files (e.g. 404 error HTML pages saved as tar)
        if (file.exists(raw_tar_dest)) file.remove(raw_tar_dest)
      }
    }, error = function(e) {
      # Remove file on error
      if (file.exists(raw_tar_dest)) file.remove(raw_tar_dest)
    })
  }
  
  # Step B: Fallback if _RAW.tar is not found (404)
  if (!download_success) {
    message(sprintf("Notice: Primary _RAW.tar not found for %s. Scanning directory for supplementary files...", gse_id))
    
    html <- tryCatch({
      readLines(suppl_dir_url, warn = FALSE)
    }, error = function(e) {
      message(sprintf("Error: Could not access directory %s", suppl_dir_url))
      NULL
    })
    
    if (!is.null(html)) {
      # Parse all href links from the HTML index page
      lines_with_href <- grep('href="', html, value = TRUE)
      hrefs <- gsub('.*href="([^"]+)".*', "\\1", lines_with_href)
      
      # Filter links to keep only files (exclude parent directories, absolute URLs, etc.)
      files_to_download <- hrefs[!grepl("/$", hrefs) & !grepl("^http", hrefs) & !grepl("^/", hrefs) & hrefs != ""]
      # De-duplicate
      files_to_download <- unique(files_to_download)
      
      if (length(files_to_download) > 0) {
        message(sprintf("Found %d supplementary file(s) for %s:", length(files_to_download), gse_id))
        for (f in files_to_download) {
          message(sprintf("  - %s", f))
        }
        
        for (file_name in files_to_download) {
          file_url <- paste0(suppl_dir_url, file_name)
          file_dest <- file.path(dest_dir, file_name)
          
          if (file.exists(file_dest) && file.info(file_dest)$size > 1000) {
            message(sprintf("File %s already exists. Skipping download.", file_name))
            next
          }
          
          message(sprintf("Downloading %s...", file_name))
          
          tryCatch({
            download.file(file_url, destfile = file_dest, method = "libcurl", mode = "wb")
            message(sprintf("Successfully downloaded: %s", file_name))
          }, error = function(e) {
            warning(sprintf("Failed to download file %s: %s", file_name, e$message))
          })
        }
      } else {
        message(sprintf("No supplementary files found in the directory for %s.", gse_id))
      }
    }
  }
}

# 4. Run downloads
for (gse_id in gse_ids) {
  download_geo_data(gse_id, dest_dir)
}
