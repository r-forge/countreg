.onAttach <- function(lib, pkg) {
  if(interactive()) {
    packageStartupMessage(
      "Starting from 'countreg' 0.3-0 the package was refactored:\n",
      "* Rootograms and other visualizations are in 'topmodels' on R-Forge.\n",
      "  https://topmodels.R-Forge.R-project.org/articles/topmodels.html\n",
      "* Distribution functions (d/p/q/r) are in 'distributions3' on CRAN.\n",
      "  https://www.zeileis.org/news/user2022/")
  }
}
