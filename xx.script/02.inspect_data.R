# GEO Data Extractor and Inspector
# Extracts downloaded raw data, decompresses any .gz files, and checks the format of all datasets

# 1. Paths Setup
# Determine script directory dynamically to make it robust
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(sub("--file=", "", file_arg))
  } else {
    # Fallback to default if run interactively in R console
    return("xx.script/02.inspect_data.R")
  }
}

script_dir <- dirname(normalizePath(get_script_path(), mustWork = FALSE))
dest_dir <- file.path(script_dir, "../00.data/GEO")

message("Data directory: ", normalizePath(dest_dir, mustWork = FALSE))

# Decompress a .gz file using pure base R connections (fast & package-free)
decompress_gz <- function(infile) {
  outfile <- sub("\\.gz$", "", infile)
  if (infile == outfile) return(FALSE) # Not a .gz file
  
  tryCatch({
    in_conn <- gzfile(infile, "rb")
    out_conn <- file(outfile, "wb")
    
    # Read and write chunks of 1MB
    while (length(buf <- readBin(in_conn, "raw", n = 1048576)) > 0) {
      writeBin(buf, out_conn)
    }
    close(in_conn)
    close(out_conn)
    
    # Remove the original .gz file to save space
    file.remove(infile)
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Failed to decompress %s: %s", infile, e$message))
    try(close(in_conn), silent=TRUE)
    try(close(out_conn), silent=TRUE)
    return(FALSE)
  })
}

# 2. Function to extract tar files
extract_tar_file <- function(gse_id, dest_dir) {
  tar_path <- file.path(dest_dir, paste0(gse_id, "_RAW.tar"))
  extract_dir <- file.path(dest_dir, gse_id)
  
  message("\n==========================================")
  message(sprintf("Extracting & Inspecting: %s", gse_id))
  message("==========================================")
  
  if (!file.exists(tar_path)) {
    warning(sprintf("Error: %s does not exist.", tar_path))
    return(NULL)
  }
  
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir, recursive = TRUE)
  }
  
  tryCatch({
    # untar extracts the tar archive
    untar(tar_path, exdir = extract_dir)
    message(sprintf("Successfully extracted tar to: %s", extract_dir))
    
    # Find all .gz files inside
    gz_files <- list.files(extract_dir, pattern = "\\.gz$", full.names = TRUE)
    
    if (length(gz_files) > 0) {
      message(sprintf("Found %d compressed (.gz) file(s). Decompressing...", length(gz_files)))
      for (gz_file in gz_files) {
        decompress_gz(gz_file)
      }
      message("All files successfully decompressed!")
    } else {
      message("No compressed (.gz) files found or already decompressed.")
    }
    
    # List and inspect uncompressed files
    uncompressed_files <- list.files(extract_dir)
    message(sprintf("Current files inside %s (Uncompressed): %d file(s)", gse_id, length(uncompressed_files)))
    message("First 5 files & previews:")
    
    for (i in 1:min(5, length(uncompressed_files))) {
      f <- uncompressed_files[i]
      f_path <- file.path(extract_dir, f)
      f_size <- file.info(f_path)$size
      
      # Determine preview type based on file extension
      if (grepl("\\.(txt|bed|tsv|csv|sam)$", f, ignore.case = TRUE)) {
        message(sprintf("  - %s (%.1f KB) [Text Data]", f, f_size / 1024))
        tryCatch({
          con <- file(f_path, "rt")
          preview_lines <- readLines(con, n = 2)
          close(con)
          for (line in preview_lines) {
            message(sprintf("    > %s", substr(line, 1, 100)))
          }
        }, error = function(e) {
          message(sprintf("    > [Preview error: %s]", e$message))
        })
      } else if (grepl("\\.cel$", f, ignore.case = TRUE)) {
        message(sprintf("  - %s (%.1f KB) [Microarray CEL Binary Raw Data]", f, f_size / 1024))
      } else {
        message(sprintf("  - %s (%.1f KB) [Other Format]", f, f_size / 1024))
      }
    }
    if (length(uncompressed_files) > 5) {
      message(sprintf("  ... and %d more files.", length(uncompressed_files) - 5))
    }
  }, error = function(e) {
    message(sprintf("Error extracting %s_RAW.tar: %s", gse_id, e$message))
  })
}

# 3. Function to preview GSE134033 data
inspect_gse134033 <- function(dest_dir) {
  csv_gz_path <- file.path(dest_dir, "GSE134033_10d_RNA-seq_genes_fpkm_annotation.csv.gz")
  xlsx_path <- file.path(dest_dir, "GSE134033_CLIP-seq_gene_expression.xlsx")
  
  message("\n==========================================")
  message("Inspecting GSE134033 Supplementary Files...")
  message("==========================================")
  
  # Preview CSV.GZ
  if (file.exists(csv_gz_path)) {
    csv_size <- file.info(csv_gz_path)$size
    message(sprintf("GSE134033 CSV.GZ Size: %.2f MB", csv_size / (1024*1024)))
    
    tryCatch({
      message("Previewing first 5 lines of GSE134033 FPKM csv.gz:")
      con <- gzfile(csv_gz_path, "rt")
      lines <- readLines(con, n = 5)
      close(con)
      
      for (i in 1:length(lines)) {
        line <- lines[i]
        message(sprintf("  Line %d: %s", i, substr(line, 1, 120)))
      }
    }, error = function(e) {
      message(sprintf("Error reading CSV.GZ: %s", e$message))
    })
  } else {
    message("GSE134033 FPKM CSV.GZ not found.")
  }
  
  # Check XLSX exists
  if (file.exists(xlsx_path)) {
    xlsx_size <- file.info(xlsx_path)$size
    message(sprintf("GSE134033 XLSX Size: %.2f MB", xlsx_size / (1024*1024)))
    message("Notice: Excel spreadsheet detected (GSE134033_CLIP-seq_gene_expression.xlsx).")
  } else {
    message("GSE134033 XLSX not found.")
  }
}

# 4. Execute extraction and inspection
extract_tar_file("GSE39872", dest_dir)
extract_tar_file("GSE37114", dest_dir)
inspect_gse134033(dest_dir)
