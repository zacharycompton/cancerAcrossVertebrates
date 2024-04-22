required_libraries <- c(
  "ape", "nlme", "rms", "phytools",
  "geiger", "caper", "tidyverse", "cowplot",
  "ggrepel", "ggsci", "patchwork", "gee",
  "poolr", "OUwie", "ggpubr", "gridGraphics",
  "devtools", "ggtree", "plotrix", "picante",
  "car", "gridExtra", "ggplot2", "grid",
  "lemon", "sqldf", "broom.mixed","huxtable","tidymodels","rr2"
)

# Function to check and install libraries
install_if_missing <- function(lib) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib, dependencies = TRUE)
    library(lib, character.only = TRUE)
  } else {
    message(paste("Library", lib, "is already installed."))
  }
}

# Check and install required libraries
invisible(sapply(required_libraries, install_if_missing))
