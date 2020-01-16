quantileNormalize <- function( x ) {
  # Normalizes raw intensity signals using 
  # quantile normalization procedure
  # INPUT
  # x - a matrix of raw intensities
  # OUTPUT
  # a matrix of quantile normalized intensities
  nloci <- dim(x)[1]
  nsamp <- dim(x)[2]

  xi <- matrix(nrow=nloci,ncol=nsamp)
  xq <- matrix(nrow=nloci,ncol=nsamp)


  for ( i in 1:nsamp) xi[,i] <- rank( x[,i], ties.method="first")
  cat("1\n")

  for ( i in 1:nsamp) xq[,i] <- sort( x[,i], na.last=TRUE )
  cat("2\n")

  m <- matrixStats::rowMedians( xq )
  cat("3\n")

  for ( i in 1:nsamp) xq[,i] <- m[  xi[,i] ] 
  cat("4\n")
  colnames(xq) <- colnames(x)
  rownames(xq) <- rownames(x)
  return(xq)
}

getCountMatrix <- function(dataObj, field = c("coverage", "numCs", "numTs")) {
  stopifnot("counts" %in% names(dataObj))
  field <- match.arg(field)
  n <- length(dataObj$counts)
  dt <- dataObj$counts[[1]] %>% S3Part(., strictS3 = TRUE) %>% setDT %>% .[, list(chr, start, end, strand, get(field))] %>% 
    setkey(chr, start, end, strand) %>% setnames(5, dataObj$sample$ID[1])
  for (i in 2:n) {
    cat(i, " ")
    dt.y <- dataObj$counts[[i]] %>% S3Part(., strictS3 = TRUE) %>% setDT %>% .[, list(chr, start, end, strand, get(field))] %>% 
      setkey(chr, start, end, strand) %>% setnames(5, dataObj$sample$ID[i])
    dt <- merge(dt, dt.y, all=TRUE)
  }
  cat("\n")
  dt
}

matrixDensityPlot <- function(X, id.vars = 1:4,
  subset = 10000,
  sampleKey
) 
{
  stopifnot("ID" %in% colnames(sampleKey))
  stopifnot("Color" %in% colnames(sampleKey))

  sample(nrow(X), min(nrow(X), subset)) %>% 
  X[., ] %>%
  reshape2::melt(id.vars=id.vars) %>% 
  setDT %>%
  setnames("variable", "ID") %>%
  merge(sampleKey, by="ID") %>%
  ggplot(., aes(value, color=Color, group=ID)) + 
  geom_density() + 
  scale_color_brewer("", palette = "Set2") + 
  scale_x_log10() + 
  theme_bw(base_size=10)
}

quantileNormalizePadlocks <- function(dataObj, totalCoverageThreshold = 30, context = "CpG") {
  message("Extracting C count matrix")
  qnCs <- getCountMatrix(dataObj, "numCs")
  message("Extracting T count matrix")
  qnTs <- getCountMatrix(dataObj, "numTs")

  # Expect the same coordinates and the same column names
  stopifnot(all(qnCs[, paste(chr, start, end, strand)] ==  qnTs[, paste(chr, start, end, strand)]))
  stopifnot(all(colnames(qnCs) == colnames(qnTs)))

  # Now normalize the matrices
  coords <- qnCs[, 1:4, with=FALSE]
  message("Quantile normalizing Cs")
  qnCs <- qnCs[, -(1:4), with=FALSE] %>% as.matrix %>% quantileNormalize
  message("Quantile normalizing Ts")
  qnTs <- qnTs[, -(1:4), with=FALSE] %>% as.matrix %>% quantileNormalize

  # Put back into the data normalized signals
  message("Writing normalized results to methylRawList")
  n <- length(dataObj$counts)
  for (i in 1:n) {
    j <- which((qnCs + qnTs)[,i] >= totalCoverageThreshold)
    dataObj$counts[[i]] <- new("methylRaw", 
      .Data = data.table(
        coords[j,], 
        coverage = (qnCs + qnTs)[j,i], 
        numCs = qnCs[j,i], 
        numTs = qnTs[j,i]),
      sample.id = dataObj$sample$ID[i],
      assembly = "hg19", 
      context = context, 
      resolution = "base",
      names = c("chr", "start", "end", "strand", "coverage", "numCs", "numTs"),
      row.names = coords[j, paste(chr, start, end, strand, sep="_")])
  }
  message("Recomputing beta values")
  dataObj$beta <- getBeta(dataObj)

  message("Done.")
  dataObj
}

