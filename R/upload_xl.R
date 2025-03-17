#' upload_xl
#' 
#' @description Imports and saves raw data generated from the N2a assay stored in `.xlsx` files. Loops through all sheets and extracts file and sheet names (96-well plate layout).
#' 
#' @param range Character. A string specifying the cell range of data to import, in the format "C46:N53".
#' 
#' @return A list of imported data (tibbles) from the file, stored in the global environment.
#' 
#' @examples
#' upload_xl(range = "C46:N53")
#' 
#' @export



upload_xl <- function(range = "C46:N53") {
  cat("Please select the Excel file.\n")
  filename <- file.choose()
  
  if (tolower(tools::file_ext(filename)) != "xls" && tolower(tools::file_ext(filename)) != "xlsx") {
    stop("Invalid filetype, please select .xlsx, or .xls-based file. Import aborted.")
  }
  
   if (!grepl("^([A-Z]+[0-9]+:[A-Z]+[0-9]+)$", range)) {
    stop("Invalid range. Use format 'C46:N53'.")
  }
  
  data_list <- list()
  
  for (sheet in excel_sheets(filename)) {
    cat("Reading sheet:", sheet, "\n")
    
    df <- tryCatch({
      suppressMessages(read_excel(filename, sheet = sheet, range = range, col_names = FALSE))
    }, error = function(e) {
      stop("Error reading sheet:", sheet, ". Import aborted.")
    })
    
    if (nrow(df) == 0 || ncol(df) == 0) {
      stop("Selected range contains no data in sheet: ", sheet, ". Import aborted.")
    }
    
    df[] <- lapply(df, function(x) {
      x <- gsub(",", ".", x)  
      numeric_x <- suppressWarnings(as.numeric(x))  
      
      if (all(is.na(numeric_x))) {  
        stop("All values in a column are non-numeric, import aborted.")
      }
      
      return(numeric_x)
    })
    
    if (nrow(df) != 8 || ncol(df) != 12) {
      stop("Datasets outside expected dimensions (8x12). Import aborted.")
    }
    
    data_list[[sheet]] <- df
  }
  
  if (!exists("upload_xl_result", envir = parent.frame())) {
    return(data_list)  
  }
  
  return(data_list)
}



  

