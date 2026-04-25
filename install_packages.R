# install_packages.R
# Run this script once to install all required R packages before launching app.R
#
#   Rscript install_packages.R

pkgs <- c(
  "shiny",
  "shinydashboard",
  "shinyWidgets",
  "colourpicker",
  "bio3d",
  "ggplot2",
  "ggrepel",
  "plotly",
  "DT",
  "dplyr",
  "tidyr",
  "scales",
  "fmsb"
)

failed <- character(0)
for (pkg in pkgs) {
  tryCatch({
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    message("  OK  ", pkg)
  }, error = function(e) {
    message("FAIL  ", pkg, ": ", conditionMessage(e))
    failed <<- c(failed, pkg)
  })
}

if (length(failed) > 0) {
  stop("The following packages could not be installed: ",
       paste(failed, collapse = ", "))
} else {
  message("All packages installed successfully.")
}
