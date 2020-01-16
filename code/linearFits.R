computeSVs <- function(Y, model, model0, n.sv = NULL, 
                      devel = FALSE, 
                      controls = NULL,
                      svaMethod = c("sva", "SmartSVA")) {
  require(sva)
  require(SmartSVA)

  svaMethod = match.arg(svaMethod)
  if (svaMethod == "sva")
    svafoo <- sva
  else
    svafoo <- smartsva.cpp

  # filter out loci where variance is 0
  v <- matrixStats::rowVars(Y, na.rm=TRUE)
  lidx <- v > 0
  Y <- Y[lidx, ]
  v <- v[lidx]

  if (devel) {
    message("DEVEL: computeSVs chosing most variable 5k loci")
    i <- order(v, decreasing = TRUE)[1:5000]
    Y <- Y[i,]
  }
  if (!is.numeric(n.sv)) {
    n.sv <- num.sv(na.omit(Y), model)
    message("computeSVs: using num.sv to guess number of SVs:", n.sv)
  } else {
    message("computeSVs: fixed number of n.sv:", n.sv)
  }
  if (is.null(controls)) {
    svobj <- svafoo(na.omit(Y), model, model0, n.sv = n.sv)  
  } else {
    controls <- controls[lidx]
    svobj <- sva(Y, model, model0, n.sv = n.sv, controls = controls)
  }
  
  res <- svobj$sv
  n <- ncol(res)
  colnames(res) <- paste("SVA", 1:n, sep="")
  res
}

permuteModelData <- function(modelData, permute, i = NULL) 
{
  modelData <- as.data.frame(modelData)
  if (is.null(i)) {
    i <- sample(nrow(modelData))
  }
  if (is.logical(permute)) {
    if (permute == TRUE) {
      message("permuteModelData: permuting all variables")
      modelData <- modelData[i,, drop=FALSE]
    } else {
      message("permuteModelData: not permuting")
    }
  } else {
    permute <- permute[ permute %in% colnames(modelData)]
    message("permuteModelData: shuffling variables ", 
      paste(permute, collapse=", "))
    modelData[, permute] <- modelData[i, permute, drop=FALSE]
  }
  return(modelData)
}


distributedChunkApply <-function(
  matrix, chunkSize, FUN, export = character(), ...) 
{
   require(itertools)
   require(foreach)
   require(bigmemory)

   sharedData <- as.big.matrix(matrix)
   desc <- describe(sharedData)
   iter <- isplitIndices(nrow(sharedData), chunkSize=chunkSize)

   res <- 
        foreach(i=iter,
                .packages=c("bigmemory"), .export=export) %dopar% 
        { 
                sharedData <- attach.big.matrix(desc)
                FUN(sharedData[i,], ...)
        }
   rm(sharedData)
   return (res)
}


# I have to reimplement mrlm to be resilient 
mymrlm <- function(M, 
                  design = NULL, 
                  ndups = 1, 
                  spacing = 1, 
                  weights = NULL, 
                  maxit=20) 
{
  if (!requireNamespace("MASS", quietly = TRUE)) 
    stop("MASS package required but is not available")
  M <- as.matrix(M)
  narrays <- ncol(M)
  if (is.null(design)) 
    design <- matrix(1, narrays, 1)
  design <- as.matrix(design)
  coef.names <- colnames(design)
  nbeta <- ncol(design)
  if (!is.null(weights)) {
    weights <- asMatrixWeights(weights, dim(M))
    weights[weights <= 0] <- NA
    M[!is.finite(weights)] <- NA
  }
  if (ndups > 1) {
    M <- unwrapdups(M, ndups = ndups, spacing = spacing)
    design <- design %x% rep(1, ndups)
    if (!is.null(weights)) 
        weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
  }
  ngenes <- nrow(M)
  stdev.unscaled <- beta <- matrix(NA, ngenes, nbeta, dimnames = list(rownames(M), 
    coef.names))
  sigma <- rep(NA, ngenes)
  df.residual <- rep(0, ngenes)
  for (i in 1:ngenes) {
    tryCatch({
      y <- as.vector(M[i, ])
      obs <- is.finite(y)
      X <- design[obs, , drop = FALSE]
      y <- y[obs]
      if (is.null(weights)) 
          w <- rep(1, length(y))
      else w <- as.vector(weights[i, obs])
      if (length(y) > nbeta) {
          out <- MASS::rlm(x = X, y = y, weights = w, maxit=maxit)
          beta[i, ] <- coef(out)
          stdev.unscaled[i, ] <- sqrt(diag(chol2inv(out$qr$qr)))
          df.residual[i] <- length(y) - out$rank
          if (df.residual[i] > 0) 
              sigma[i] <- out$s
      }
    }, error = function(e) {
      message("Error caught. Continuing...")
    })
  }
  QR <- qr(design)
  cov.coef <- chol2inv(QR$qr, size = QR$rank)
  est <- QR$pivot[1:QR$rank]
  dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
  list(coefficients = beta, stdev.unscaled = stdev.unscaled, 
      sigma = sigma, df.residual = df.residual, cov.coefficients = cov.coef, 
      pivot = QR$pivot, rank = QR$rank)
}

mylmFit <- function(object, 
                    design = NULL, 
                    ndups = 1, 
                    spacing = 1, 
                    block = NULL, 
                    correlation, 
                    weights = NULL, 
                    method = "ls", 
                    maxit = 20)
{
  require(limma)
    y <- getEAWP(object)
    if (is.null(design)) 
        design <- y$design
    if (is.null(design)) 
        design <- matrix(1, ncol(y$exprs), 1)
    else {
        design <- as.matrix(design)
        if (mode(design) != "numeric") 
            stop("design must be a numeric matrix")
        if (nrow(design) != ncol(y$exprs)) 
            stop("row dimension of design doesn't match column dimension of data object")
    }
    ne <- nonEstimable(design)
    if (!is.null(ne)) 
        cat("Coefficients not estimable:", paste(ne, collapse = " "), 
            "\n")
    if (missing(ndups) && !is.null(y$printer$ndups)) 
        ndups <- y$printer$ndups
    if (missing(spacing) && !is.null(y$printer$spacing)) 
        spacing <- y$printer$spacing
    if (missing(weights) && !is.null(y$weights)) 
        weights <- y$weights
    method <- match.arg(method, c("ls", "robust"))
    if (ndups > 1) {
        if (!is.null(y$probes)) 
            y$probes <- uniquegenelist(y$probes, ndups = ndups, 
                spacing = spacing)
        if (!is.null(y$Amean)) 
            y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean), 
                ndups = ndups, spacing = spacing), na.rm = TRUE)
    }
    if (method == "robust") 
        fit <- mymrlm(y$exprs, design = design, ndups = ndups, 
            spacing = spacing, weights = weights, maxit=maxit)
    else if (ndups < 2 && is.null(block)) 
        fit <- lm.series(y$exprs, design = design, ndups = ndups, 
            spacing = spacing, weights = weights)
    else {
        if (missing(correlation)) 
            stop("the correlation must be set, see duplicateCorrelation")
        fit <- gls.series(y$exprs, design = design, ndups = ndups, 
            spacing = spacing, block = block, correlation = correlation, 
            weights = weights, ...)
    }
    if (NCOL(fit$coef) > 1) {
        n <- rowSums(is.na(fit$coef))
        n <- sum(n > 0 & n < NCOL(fit$coef))
        if (n > 0) 
            warning("Partial NA coefficients for ", n, " probe(s)", 
                call. = FALSE)
    }
    fit$genes <- y$probes
    fit$Amean <- y$Amean
    fit$method <- method
    fit$design <- design
    new("MArrayLM", fit)
}


parallel_mylmFit <- function(M, 
                            model, 
                            chunkSize=50000, 
                            nNodes=4, ...) 
{
  require(limma)
  # compute many fits as a list
  fit <- distributedChunkApply(matrix=M, chunkSize=chunkSize,
    FUN=mylmFit, export=c("mymrlm"), design=model, ...)
  # combine into one fit object
  res <- list()
  for( item in names(fit[[1]])) {
    res[[item]] <- DistributedLmFit$mergeListItem(fit, item)
  }
  class(res) <- class(fit[[1]])
  return(res)
}

getLimmaFit <- function(M, 
                        md, 
                        formula, 
                        formula0, 
                        contrasts = NULL,
                        nsv, 
                        iteration, 
                        distributed=TRUE,
                        ...) 
{  
  message(glue("############# n.sv {nsv} Iteration {iteration} ##############"))
  if (iteration > 0)
    md <- permuteModelData(md, TRUE)
  m <- model.matrix(formula %>% as.formula, md)
  m0 <- model.matrix(formula0 %>% as.formula, md)
  coef <- which( !colnames(m) %in% colnames(m0))

  # Remove colons from design matrices
  if (any(grepl(":", colnames(m)))) {
    colnames(m) <- make.names(colnames(m))
    colnames(m0) <- make.names(colnames(m0))
  }

  if (nsv > 0) {
    svs <- computeSVs(M, m, m0, n.sv = nsv, ...)
    m <- cbind(m, svs)
    m0 <- cbind(m0, svs)
  }

  # Now run lmFit with eBayes
  if (distributed) {
    fit <- parallel_mylmFit(
      M = M,
      model = m, 
      chunkSize = ceiling(nrow(M)/params$workers),
      nNodes = params$workers,
      method = "robust",
      maxit = 100)
  } else {
    fit <- mylmFit(M, design = m, method = "robust", maxit = 100)
  }

  if (!is.null(contrasts)) {
    cont <- makeContrasts( 
      contrasts = contrasts,
      levels    = m
    )
    fit <- contrasts.fit(fit, cont) %>% eBayes 
    testF <-  
      topTable(fit, sort.by="none", number = Inf) %>% 
        setDT %>% 
        .[, list(P.Value, adj.P.Val)]
    n <- ncol(fit$coef)
    fit <- 
      foreach (i = 1:n, .combine = cbind) %do% {
        contrast <- colnames(fit$coef)[i]
        dt <- topTable(fit, coef = i, 
          sort.by = "none", number = Inf)
        dt <- data.table(C = dt$logFC, P = dt$P.Value)
        setnames(dt, 
          c("C", "P"), 
          paste0(c("C.", "P."), contrast))
        dt
      } %>% 
      data.table(ID = rownames(fit$coef), .) %>%
      cbind(., testF  )
  } else {
    fit <- fit %>% eBayes %>% 
      topTable(coef=coef, number=Inf, sort.by="none") %>% 
      data.table(ID = rownames(.), .)
  }
  return(fit)
}

permutationAnalysis <- function(M, 
                                md, 
                                formula, 
                                formula0, 
                                contrasts = NULL,
                                nsv, 
                                nPerm,
                                ...) 
{
  # Perform permutation analysis
  export = c("getLimmaFit", "parallel_mylmFit",
             "mylmFit", "mymrlm", 
             "computeSVs", "permuteModelData")
  permuted <- foreach( iteration = 0:nPerm, .combine = c, .export=export) %dopar% {
    tryCatch({
      source("../code/common.R")
      fit <- getLimmaFit(M, md, formula, formula0, contrasts, nsv, iteration, 
        distributed=FALSE, ...)
      fit[, sum(adj.P.Val < 0.05, na.rm=TRUE)]
    }, error = function(e) {
      e
    })
  }
  permuted
}
