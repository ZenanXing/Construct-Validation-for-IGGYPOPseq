
# Function to check and install packages
install_if_needed <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
}

# Install BiocManager if not installed
install_if_needed('BiocManager')
if (!requireNamespace('cpp11', quietly = TRUE)) {
  BiocManager::install("cpp11")
} 

# Install CRAN packages if not installed
#packages1 <- c('cpp11', 'tzdb', 'timechange', 'systemfonts', 'textshaping', 'ragg', 'gtable', 'purrr')
#packages2 <- c('tidyverse', 'ggplot2', 'cowplot', 'googledrive', 'googlesheets4', 'readr', 'readxl', 'tidyr', 'broom', 'dbplyr', 'haven', 'modelr')
#install_if_needed(packages1)
#install_if_needed(packages2)

packages <- c("waldo", "tzdb", "systemfonts", "timechange", "textshaping", "purrr", "gtable",
              "tidyverse", "Biostrings", "openxlsx", "seqinr",   
              "stringr", "ggplot2", "cowplot", "rmarkdown")  
install_if_needed(packages)


# Install Biostrings if not installed
if (!requireNamespace('Biostrings', quietly = TRUE)) {
  BiocManager::install('Biostrings')
} 

