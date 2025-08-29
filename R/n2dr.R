#' n2dr
#' 
#' @description Dose-response analysis by the 4PL model of N2a assay data. Generates descriptive summaries of parameter estimates and plots curves.
#' 
#' @param datalist List of data (imported through `upload_xl()` or `upload_txt()`).
#' @param stock Numeric. Stock concentration of sample or standard.
#' @param dose Numeric. Sample aliquot volume of sample or standard.
#' @param tissue Character. Tissue concentration conversion, either "dry" (ng/muL) or "liquid" (ng/mL) (default is "liquid").
#' @param x.axis Character. Title for the x-axis (default is "").
#' @param y.axis Character. Title for the y-axis (default is "% cell survival").
#' @param subset Numeric. A subset of datasets to analyze (e.g., `subset = 1:3`). Default is `NULL`, which loops through all.
#' 
#' @return A list containing dose-response curves (4PL) and parameter estimates.
#' 
#' @examples
#' n2dr(exampledata, stock = 2.5, dose = 20, subset = 1)
#' 
#' @export


n2dr <- function(datalist, stock, dose, tissue = "liquid", 
                 x.axis = "", y.axis = "", subset = NULL) {
  if (missing(stock) || missing(dose)) {
    stop("Both 'stock' and 'dose' must be provided.")
  }
  suppressWarnings({ 
    
    if (length(datalist) == 0 || all(sapply(datalist, function(x) dim(x)[1] == 0))) {
      stop("Empty dataset: No data available for analysis.")
    }
    
    names_list <- names(datalist)
    
    if (is.null(subset)) {
      subset <- seq_along(datalist) 
    } else {
      subset <- as.integer(subset)  
      if (any(subset < 1 | subset > length(datalist))) {
        stop("Subset indices are out of range for the datalist.")
      }
    }
    
    datalist <- datalist[subset]
    names_list <- names_list[subset]
    
    stock <- if (length(stock) == length(datalist)) stock else rep(stock, length.out = length(datalist))
    dose <- if (length(dose) == length(datalist)) dose else rep(dose, length.out = length(datalist))
    x.axis.title <- rep(x.axis, length.out = length(datalist))
    y.axis.title <- rep(y.axis, length.out = length(datalist))
    
    calc_conc <- function(tissue, dose_value, stock_value) {
      if (tissue == "dry") {
        return(((dose_value * stock_value) / 200) * 10 * (0.5^(0:7)))
      } else if (tissue == "liquid") {
        return((((dose_value * stock_value) / 200) * 10 * (0.5^(0:7)) / 230) * 1000)
      } else {
        stop("tissue must be defined as either 'dry' or 'liquid'")
      }
    }
    
    norm_values <- function(data, indices, control_mean) {
      norm <- (t(data[indices, 4:11]) / control_mean) * 100
      return(list(mean = rowMeans(norm, na.rm = TRUE), std = apply(norm, 1, sd, na.rm = TRUE)))
    }
    
    estimate_start <- function(x, y) {
      min_val <- min(y, na.rm = TRUE)
      max_val <- max(y, na.rm = TRUE)
      e50 <- median(x)
      slope <- -1
      return(c(min = min_val, max = max_val, EC50 = e50, Hillslope = slope))
    }
    
    for (i in seq_along(datalist)) {
      data <- datalist[[i]]
      dose_value <- dose[i]
      stock_value <- stock[i]
      x_title <- x.axis.title[i]
      y_title <- y.axis.title[i]
      main_title <- names_list[i]
      
      conc <- calc_conc(tissue, dose_value, stock_value)
      control_mean_min <- mean(as.numeric(as.matrix(data[2:4, 2:3])), na.rm = TRUE)
      control_mean_plus <- mean(as.numeric(as.matrix(data[5:7, 2:3])), na.rm = TRUE)
      norm_min <- norm_values(data, 2:4, control_mean_min)
      norm_plus <- norm_values(data, 5:7, control_mean_plus)
      ov <- mean((control_mean_plus / control_mean_min) * 100, na.rm = TRUE)
      
      fit_4pl <- function(conc, norm_plus) {
        fourpl <- function(x, min, max, EC50, Hillslope) {
          min + (max - min) / (1 + (x / EC50)^(-Hillslope))
        }
        
        start_vals <- estimate_start(conc, norm_plus$mean)
        
        lower_bounds <- c(min = -10, max = -10, EC50 = 0, Hillslope = -5)  
        upper_bounds <- c(min = 200, max = 200, EC50 = 50, Hillslope = 5)  
        
        obj_func <- function(par) {
          if (any(is.nan(par)) || any(is.infinite(par))) {
            return(1e10)  
          }
          
          y_pred <- fourpl(conc, par[1], par[2], par[3], par[4])
          
          if (any(is.nan(y_pred)) || any(is.infinite(y_pred))) {
            return(1e10) 
          }
          
          sum((norm_plus$mean - y_pred)^2)
        }
        
        fit <- nlminb(start = start_vals, objective = obj_func, lower = lower_bounds, upper = upper_bounds)
        
        if (any(is.nan(fit$par)) || any(is.infinite(fit$par))) {
          stop("Optimization failed due to NaNs or Infs in the parameters.")
        }
        
        hessian <- numDeriv::hessian(func = obj_func, x = fit$par)
        cov_matrix <- solve(hessian)
        SE <- sqrt(diag(cov_matrix))  
        
        return(list(fit = fit, SE = SE))
      }
      
      fit_mean <- fit_4pl(conc, norm_plus)
      
      params_mean <- fit_mean$fit$par
      SE_params <- fit_mean$SE
      
      y_max <- max(c(norm_plus$mean + norm_plus$std, norm_min$mean + norm_min$std), na.rm = TRUE)
      y_min <- min(c(norm_plus$mean - norm_plus$std, norm_min$mean - norm_min$std), na.rm = TRUE)
      y_padding <- 0.15 * (y_max - y_min)
      
      par(bty = "l", mar = c(5, 5, 4, 4) + 0.1, lwd = 2, tck = -0.02, cex.axis = 1.2)
      
      plot(conc, norm_plus$mean, log = "x", type = "n",
           xlab = x_title, ylab = y_title,
           ylim = c(y_min - y_padding, y_max + 1.2 * y_padding),
           main = main_title, cex.lab = 1.5, cex.main = 1.5, las = 1)
      
      axis(1)
      maj_ticks <- axTicks(1)
    
      minor_ticks <- unlist(sapply(maj_ticks, function(x) x * (2:9)))
      minor_ticks <- minor_ticks[minor_ticks > min(conc) & minor_ticks < max(conc)]
      axis(1, at = minor_ticks, labels = FALSE, tcl = -0.3)
      
      if (any(norm_min$std > 0, na.rm = TRUE)) {
        arrows(conc, norm_min$mean - norm_min$std, conc, norm_min$mean + norm_min$std, 
               length = 0.1, angle = 90, code = 3, col = "gray67")
      }
      if (any(norm_plus$std > 0, na.rm = TRUE)) {
        arrows(conc, norm_plus$mean - norm_plus$std, conc, norm_plus$mean + norm_plus$std, length = 0.1, angle = 90, code = 3)
      }
      
            xval <- seq(min(conc), max(conc), length.out = 1000)
      yval <- params_mean[1] + (params_mean[2] - params_mean[1]) / (1 + (xval / params_mean[3])^(-params_mean[4]))
      lines(xval, yval, lwd = 3) 
      
      points(conc, norm_min$mean, pch = 16, col = "gray67", cex = 1.3)
      points(conc, norm_plus$mean, pch = 21, col = "black", bg = "white", cex = 1.6)
      
      legend_y <- max(norm_min$mean, na.rm = TRUE) + 1.5 * y_padding 
      legend("topright", 
             legend = c("-OV", "+OV"), 
             pch = c(16, 21),  
             col = c("gray67", "black"), 
             pt.bg = c(NA, "white"),  
             bty = "n", 
             cex = 1.2, 
             xpd = TRUE, 
             y.intersp = 0.8, 
             text.col = "black")
      
      n2.sum <- data.frame(
        "Parameters" = c("EC50", "max", "min", "Hillslope"),
        "Mean" = round(c(params_mean[3], params_mean[2], params_mean[1], params_mean[4]), 2),
        "SE" = round(SE_params, 2),  
        stringsAsFactors = FALSE
      )
      
      
      cat(sprintf("Summary for: %s\n", main_title))
      cat(sprintf("stock: %s, dose: %s\n", stock_value, dose_value))
      print_summary <- function(n2.sum, indentation) {
        cat(sprintf("%s%-12s %-6s %-6s\n", indentation, "[Parameters]", "[Mean]", " SE"))
        cat(sprintf("%s%s\n", indentation, strrep("-", 9 + 9 + 9)))
        for (row in seq_len(nrow(n2.sum))) {
          cat(sprintf("%s%-12s %-6.2f  %-3.2f\n", indentation,
                      n2.sum$Parameters[row], n2.sum$Mean[row], n2.sum$SE[row]))
        }
      }
      
    
      print_summary(n2.sum, paste(rep(" ", 10), collapse = ""))
      cat(sprintf("%s%s\n", paste(rep(" ", 10), collapse = ""), strrep("-", 9 + 9 + 9)))
      if (ov >= 60 && ov <= 80) {
        cat(sprintf("OV%%: %.2f -> Within optimal range (60-80%%)\n", round(ov, 2)))
      } else {
        cat(sprintf("* OV%%: %.2f -> Outside optimal range (60-80%%)\n", round(ov, 2)))
      }
      
      cv_min <- norm_min$std / norm_min$mean * 100
      cv_plus <- norm_plus$std / norm_plus$mean * 100
      
      if (any(cv_min > 20, na.rm = TRUE) || any(cv_plus > 20, na.rm = TRUE)) {
        if (any(cv_min > 20, na.rm = TRUE)) {
          high_cv_rows <- which(cv_min > 20)
          cat(sprintf("* CV above 20%% | Sample %s (-OV)\n", paste(high_cv_rows, collapse = ", ")))
        }
        if (any(cv_plus > 20, na.rm = TRUE)) {
          high_cv_rows <- which(cv_plus > 20)
          cat(sprintf("* CV above 20%% | Sample %s (+OV)\n", paste(high_cv_rows, collapse = ", ")))
        }
      }
      for (i in seq_along(norm_min$mean)) {
        other_means <- norm_min$mean[-i]
        if (all(norm_min$mean[i] < 0.3 * other_means)) {
          cat("* Matrix effect detected\n")
          
        }
      }
      
    }
  }) 
}


