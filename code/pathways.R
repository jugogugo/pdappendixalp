# Human
getGoGenesHuman <- function(term, localhost=FALSE) {
	if (localhost) {
	  url <- glue("http://localhost:8080/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity_name,bioentity,bioentity_label&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=regulates_closure:%22{term}%22&fq=document_category:%22bioentity%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&facet.field=source&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=panther_family_label&facet.field=annotation_class_list_label&facet.field=regulates_closure_label&q=*:*")
	} else {
	  url <- glue("http://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity_name,bioentity,bioentity_label&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=regulates_closure:%22{term}%22&fq=document_category:%22bioentity%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&facet.field=source&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=panther_family_label&facet.field=annotation_class_list_label&facet.field=regulates_closure_label&q=*:*")
	}
  genes <- fread(url, header = FALSE)
  return(genes$V3)
}

# Mouse
getGoGenesMouse <- function(term) {
  require(glue)
  url <- "http://www.informatics.jax.org/go/report.txt?goID={term}&results=1000&startIndex=0&sort=term&dir="
  x <- fread(glue(url))
  x <- x[,2, with=FALSE]
  setnames(x, "Gene")
  unique(x$Gene)
}


ALP <- list(
  list(term="GO:0061912", name="Selective autophagy"),
  list(term="GO:0016236", name="Macroautophagy"),
  list(term="GO:0005764", name="Lysosomes"),
  list(term="GO:0061684", name="Chaperone-mediated autophagy"),
  list(term="GO:0006914", name="Autophagy")
  )

getHumanPathways <- function(pathways=ALP, localhost=FALSE) {
	require(foreach)
	foreach(p = pathways) %do% {
  	p$genes <- getGoGenesHuman(p$term, localhost=localhost)
  	p
	}
} 

getMousePathways <- function(pathways=ALP) {
	require(foreach)
	foreach(p = pathways) %do% {
	  p$genes <- getGoGenesMouse(p$term)
  	p
	}
}