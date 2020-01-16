source("utils.R")

DistributedLmFit <- list()

DistributedLmFit$mergeListItem <- function(mylist, item) {
	c_items <- c("sigma", "df.residual", "Amean", "pivot")
	rbind_items <- c("coefficients", "stdev.unscaled")
	first_items <- c("design", "method", "cov.coefficients")
	
	if (item %in% c_items)
		do.call("c", lapply(mylist, '[[', item))
	else if (item %in% rbind_items)
		do.call("rbind", lapply(mylist, '[[', item))
	else if (item %in% first_items) 
		mylist[[1]][[item]]
}

DistributedLmFit$lmFit <- function(M, model, chunkSize=50000, nNodes=4, ...) {
	require(limma)
	# compute many fits as a list
	fit <- distributedChunkApply(matrix=M, chunkSize=chunkSize, nNodes=nNodes,
		FUN=lmFit, export=character(), model, ...)
	# combine into one fit object
	res <- list()
	for( item in names(fit[[1]])) {
		res[[item]] <- DistributedLmFit$mergeListItem(fit, item)
	}
	class(res) <- class(fit[[1]])
	return(res)
}
