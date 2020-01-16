require("shiny")
require("knitr")
require("kableExtra")
require("data.table")
require("dplyr")
require("ggplot2")
require("plotly")
require("ggpubr")
require("foreach")
require("glue")
require("WGCNA")
require("impute")
require("limma")
require("itertools")
require("parallel")
require("doSNOW")
require("latex2exp")
source("../code/plots.R")
source("../code/utils.R", chdir=TRUE)
source("../code/DistributedLmFit.R", chdir=TRUE)
set.seed(1234)

knitr::opts_knit$set(
  root.dir = getwd()
)

knitr::opts_chunk$set(
  context   = "render",
  echo      = FALSE, 
  include   = FALSE, 
  out.width = "100%", 
  fig.width = 10, 
  fig.height= 6)


RECOMPUTE <<- FALSE

touch <- function() {
  RECOMPUTE <<- TRUE
}

cache <- function(condition=RECOMPUTE, foo=foo, fname, verbose=FALSE, ...) {
  # If condition evaluates to TRUE, then 
  # recompute using the call to function
  if (condition | !file.exists(fname)) {
    if(verbose)
      message("Recomputing.")
    res <- foo(...)
    saveRDS(res, fname)
    RECOMPUTE <<- FALSE
  } else {
    if (verbose)
      message("Using cached result.")
    res <- readRDS(fname)
  }
  return (res)
}

