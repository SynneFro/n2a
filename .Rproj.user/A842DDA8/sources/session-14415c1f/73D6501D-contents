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
  
  if (!grepl("^([A-Z]+[0-9]+:[A-Z]+[0-9]+)$", range)) {
    stop("Invalid range. Use format 'C46:N53'.")
  }
  
  data_list <- list()
  
  for (sheet in excel_sheets(filename)) {
    cat("Reading sheet:", sheet, "\n")
    df <- tryCatch({
      suppressMessages(read_excel(filename, sheet = sheet, range = range, col_names = FALSE))
    }, error = function(e) {
      message("Error reading sheet:", sheet)
      return(NULL)
    })
    
    if (nrow(df) == 0 || ncol(df) == 0) {
      stop("Selected range contains no data in sheet: ", sheet)
    }
    
    data_list[[sheet]] <- df
  }
  
  cat("Data imported and stored.\n")
  return(data_list)
}



  

