catln <- function(...) { cat(...,"\n")}
assert <- function(...) { stopifnot(...)}
spaste <- function(...) { paste(..., sep="")}
less <- function( ...) { page(..., method="print")}

newDev <- function (width=6, height=6) 
{
	if ( length(dev.list()) > 0)
		dev.off()
	dev.new(width = width, height = height)
}

### subset tables inside the list based on the number of rows/columns

shrinkList <- function(myData, inds, ref=dim(myData$qntl)[1], verbose=1) {

	untouchedTables <- character()
	
	for(i in 1:length(myData)) {
	
		dims <- dim(myData[[i]])
		len <- length(myData[[i]])
			
		if(length(dims)!=0 && dims[1]==ref && dims[2]==ref) {
			stop("Cannot shrink - there is a table with dim_1 = dim_2 = reference")
		}
		
		if(length(dims)!=0 && dims[1]==ref) {
			myData[[i]] <- myData[[i]][inds,]
		} else if(length(dims)!=0 && dims[2] == ref) {
			myData[[i]] <- myData[[i]][,inds]
		} else if(len==ref)	{
			myData[[i]] <- myData[[i]][inds]
		} else {
			untouchedTables <- c(untouchedTables, names(myData)[i])
		}
	}
	if(verbose == 1) {
		print(paste("untouched tables:", untouchedTables))
	}
	return(myData)
}


### apply-type function for 2 matrices.
### returns function applied for each related row pair of mat1 and mat2

pairapply <- function(mat1, mat2, fun) {

	if(dim(mat1)[1]!=dim(mat2)[1] | dim(mat1)[2]!=dim(mat2)[2]) {
		stop("dimentions of matrices are not equal")
	}
	
	res <- sapply(1:dim(mat1)[1], function(i) fun(mat1[i,], mat2[i,]))
	return(res)
}


### function to get column and row indeces based on matrix index

matrixToInds <- function(matrixInds, matrixDim) {
	rowNo <- ((matrixInds-1) %/% matrixDim[2])+1
	colNo <- (matrixInds - (matrixDim[2]*rowNo) + matrixDim[2])
	return(cbind(rowNo, colNo))
}



### Garbage collection from Horvath
collectGarbage <- function () 
{
	while (gc()[2, 4] != gc()[2, 4] | gc()[1, 4] != gc()[1, 4]) 
	{
	}
}



getObject <- function(str) {
	# Given a variable name as a string, 
	# return the variable
	eval(parse(text =str))
}


loadObject <- function(f)
{
	# Given a filename, load it's contents into a variable
    env <- new.env()
    nm <- load(f, env, verbose=TRUE)[1]
    env[[nm]]
}



## Creates a parallel Cluster. Very useful on scinet
CLUSTER <- NULL
makeMyCluster <- function(outfile="", nNodes=0) {
  require(parallel)
  require(doSNOW)
  if (nNodes == 0) {
      nodefile <- Sys.getenv("PBS_NODEFILE")
      hosts <- readLines(nodefile)
  } else {
      hosts <- rep("localhost", nNodes)
  }
  message(sprintf("Starting cluster on %s", paste(hosts, collapse=", ")))
  cluster <- makeCluster(hosts, type="SOCK", outfile=outfile)    
  registerDoSNOW(cluster)
  clusterSetupRNG(cluster)
  return(cluster)
}
stopMyCluster <- function(cluster) {
  require(parallel)
  require(doSNOW)
  message("Stopping cluster")
  registerDoSEQ()
  stopCluster(cluster)
}
withCluster <- function(action, ...) {
  CLUSTER <<- makeMyCluster(...)  
  tryCatch(action, finally={
    stopMyCluster(CLUSTER)
    CLUSTER <<- NULL
  })
}


# This function will split vector x into groups of size n
chunks <- function(x, n, force.number.of.groups = FALSE, len = length(x), groups = trunc(len/n), overflow = len%%n) { 
  if(force.number.of.groups) {
    f1 <- as.character(sort(rep(1:n, groups)))
    f <- as.character(c(f1, rep(n, overflow)))
  } else {
    f1 <- as.character(sort(rep(1:groups, n)))
    f <- as.character(c(f1, rep("overflow", overflow)))
  }

  g <- split(x, f)

  if(force.number.of.groups) {
    g.names <- names(g)
    g.names.ordered <- as.character(sort(as.numeric(g.names)))
  } else {
    g.names <- names(g[-length(g)])
    g.names.ordered <- as.character(sort(as.numeric(g.names)))
    g.names.ordered <- c(g.names.ordered, "overflow")
  }

  return(g[g.names.ordered])
}


# Split a given matrix into chunks by rows and apply a function
# on each chunk in parallel
# Params:
# - matrix - input matrix
# - chunkSize - split by this number 
# - nNodes - number of nodes for distribution
# - FUN - the function to apply
# - ... - other params passed on to the function
distributedChunkApply <-function(
  matrix, chunkSize, nNodes, FUN, export = character(), ...) 
{
   require(itertools)
   require(foreach)
   require(bigmemory)

   sharedData <- as.big.matrix(matrix)
   desc <- describe(sharedData)
   iter <- isplitIndices(nrow(sharedData), chunkSize=chunkSize)

   res <- withCluster(
        foreach(i=iter,
                .packages=c("bigmemory"), .export=export) %dopar% 
        { 
                sharedData <- attach.big.matrix(desc)
                FUN(sharedData[i,], ...)
        }, nNodes=nNodes)
   rm(sharedData)
   return (res)
}


# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)

  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")

  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }

  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)

  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin

  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

# Prints and returns the same object, good for piping
tee <- function(obj, foo=print) { foo(obj); obj }
