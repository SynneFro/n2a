\name{upload_xl}
\alias{upload_xl}
\title{Import and Save Raw Data from N2a Assay Excel Files}
\description{
  Imports and saves raw data generated from the N2a assay stored in `.xlsx` files. Loops through all sheets and extracts file and sheet names (96-well plate layout).
}
\usage{
  upload_xl(range = "C46:N53")
}
\details{
  The function reads an `.xlsx` file selected by the user and imports data from all sheets within the specified cell range. The default range is "C46:N53". The data is imported into a list, with each sheet being stored as a separate element.
}

\arguments{
  \item{range}{A character string specifying the cell range of data to import, in the format "C46:N53".}
}

\author{
  Synne Thorbjørnsen Frøstrup
}
\note{
  The function assigns the imported data to a variable in the global environment, named after the `.xlsx` file (without the extension).
}

