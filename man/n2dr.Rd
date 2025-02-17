\name{n2dr}
\alias{n2dr}
\title{Dose-response analysis using the 4PL model for N2a assay data}
\description{
  Dose-response analysis by the 4PL model of N2a assay data. Generates descriptive summaries of parameter estimates and plots curves.
}
\usage{
n2dr(datalist, stock, dose, tissue = "liquid", 
                 x.axis = "", y.axis = "", subset = NULL)
}
\arguments{
  \item{datalist}{
    A list of data imported through \code{upload_xl()} or \code{upload_txt()}.
  }
  \item{stock}{
    Numeric. Stock concentration of sample or standard.
  }
  \item{dose}{
    Numeric. Sample aliquot volume of sample or standard.
  }
  \item{tissue}{
    Character. Tissue concentration conversion, either "dry" (ng/muL) or "liquid" (ng/mL), default is "liquid".
  }
  \item{x.axis}{
    Character. Title for the x-axis, default is "".
  }
  \item{y.axis}{
    Character. Title for the y-axis, default is "\% Cell Survival".
  }
  \item{subset}{
    Numeric. A subset of datasets to analyze (e.g., \code{subset = 1:3}). Default is \code{NULL}, which loops through all datasets.
  }
}
\details{
  The function performs a dose-response analysis using the four-parameter logistic (4PL) model on the provided N2a assay data. It generates parameter estimates, plots dose-response curves, and provides descriptive summaries including EC50 values and other model parameters. The function also allows customization of axis labels and subset selection for analysis.
}
\value{
  A list containing dose-response curves (4PL) and parameter estimates.
  \item{fit}{A list with the optimized parameters from the 4PL model fitting.}
  \item{SE}{Standard errors for the parameter estimates.}
}

\author{
  Synne Thorbjørnsen Frøstrup
}
\note{
  The function requires that both \code{stock} and \code{dose} are provided for analysis.
}

