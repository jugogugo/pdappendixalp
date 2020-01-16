
require("latex2exp")

# Takes a p value and formats with Latex in scientific notation
# but only if the value is less than 0.01
texP <- function(p, digits=3, scthreshold = 0.01) {
  p[ p == 0 ] <- 2.2e-16
  scientific <- p < scthreshold
  scientific[is.na(p)] <- FALSE
  output <- character(length(p))

  output[!scientific] <- 
  	format(p[!scientific], digits=digits)

	#Transforms the number into scientific notation even if small
  output[scientific] <- 
  	format(p[scientific], digits = digits, scientific = TRUE) 
  output[scientific] <- strsplit(output[scientific], split = "e")
  output[scientific] <- 
  	sapply(output[scientific], 
  	       function(X) {
  	       	sprintf("$%s\\,\\mathrm{x}\\,10^{%s}$", 
  	       	        X[1], as.integer(X[2]))
  	       })
  output[!scientific] <- 
  	sapply(output[!scientific],
  				 function(X) {
  				 	sprintf("$%s$", X)
  				 })
  output[ is.na(p) ] <- "NA"
  return(unlist(output))
}


# texP(c(0.1, 0.05, 0.01, NA, 0.0001), digits=1)


texPstars <- function(pstars) {
	factor(pstars, 
         levels = c("***", "**", "*", ".", " "),
         labels = c("\\*\\*\\*", "\\*\\*", "\\*", "$\\cdot$", " ")
       )
}


# ps  <- c(0.2, 0.1, 0.05, 0.01, NA, 0.0001)
# texPstars(gtools::stars.pval(ps))