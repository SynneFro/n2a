#' n2a_report
#'
#' @description Renders a report of N2a dose-response analysis: one dose-response
#'   plot and parameter summary per plate, followed by the full results table.
#'   Wraps an R Markdown template; requires the 'rmarkdown' and 'knitr' packages.
#'
#' @param datalist List of imported data (from `upload_xl()` or `upload_txt()`).
#' @param stock Numeric. Stock concentration of sample or standard (passed to `n2dr()`).
#' @param dose Numeric. Sample aliquot volume (passed to `n2dr()`).
#' @param well.vol Numeric. Final well volume in muL (default 230).
#' @param tissue Character. "dry" or "liquid" (default "liquid").
#' @param output Character. Output file path (default "n2a_report.html").
#' @param title Character. Report title (default "N2a Dose-Response Report").
#'
#' @return (Invisibly) the path to the rendered report.
#'
#' @examples
#' \dontrun{
#' plates <- upload_xl()
#' n2a_report(plates, stock = 20, dose = 1.5)
#' }
#'
#' @export
n2a_report <- function(datalist, stock, dose, well.vol = 230, tissue = "liquid",
                       output = "n2a_report.html",
                       title = "N2a Dose-Response Report") {
  
  if (!requireNamespace("rmarkdown", quietly = TRUE) ||
      !requireNamespace("knitr", quietly = TRUE)) {
    stop("The 'rmarkdown' and 'knitr' packages are required for n2a_report(). ",
         "Install them with install.packages(c('rmarkdown','knitr')).")
  }
  if (missing(stock) || missing(dose)) {
    stop("Both 'stock' and 'dose' must be provided.")
  }
  
  template <- system.file("rmarkdown", "n2a_report.Rmd", package = "n2a")
  if (template == "") {
    stop("Report template not found. Reinstall the n2a package.")
  }
  
  output <- normalizePath(output, mustWork = FALSE)
  
  rmarkdown::render(
    input = template,
    output_file = output,
    params = list(
      datalist = datalist,
      stock = stock,
      dose = dose,
      well.vol = well.vol,
      tissue = tissue,
      title = title
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  
  cat(sprintf(">> Report saved to: %s\n", output))
  invisible(output)
}
