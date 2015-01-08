########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Define the GO terms using biomaRt and create a GMT file for GSEA
##
## inputs:	
##      - 
##
## outputs, a list of 2 items:
##			- 
##
## Notes:	duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
#################################################

`formatGMT` <-
function (infile, outfile, replace=FALSE, verbose=TRUE) {
	## this function read a gmt file and replace all weird characters in gene set names by "_"

  if(file.exists(outfile)) {
    if(!replace) { stop("Output file already exists!") }
    file.remove(outfile)
  }
  
  badchars <- "[,][:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
	if(verbose) { message(sprintf("reading %s", basename(infile))) }
  tempff <- file.path(dirname(outfile), basename(tempfile(pattern="formatGMT_", tmpdir="", fileext=".tmp")))
  sink(file=tempff, type="output")
  rr <- GSA::GSA.read.gmt(filename=infile)
  sink()
  unlink(x=tempff, force=TRUE)
  dupln <- sum(duplicated(rr$geneset.names))
  rr$geneset.names <- gsub(pattern=badchars, replacement="_", x=toupper(rr$geneset.names))
  dupln <- c(dupln, sum(duplicated(rr$geneset.names)))
  if(dupln[1] > 0) { warning("Duplicated geneset names in the original gmt file") }
  if(any(dupln[-1] > dupln[1])) { warning("duplicated geneset names due to formatting") }
    names(rr$genesets) <- names(rr$geneset.descriptions) <- rr$geneset.names
  golist <- mapply(c, rr$geneset.names, rr$geneset.descriptions, rr$genesets)  
  
	## write gmt file
	if(verbose) { message(sprintf("writing %s to %s", basename(outfile), dirname(outfile))) }
	rr2 <- lapply(golist, function(x, file) {
    write(paste(c(x[1], x[2], unique(x[3:length(x)])), collapse="\t"), file=file, append=TRUE)
  }, file=path.expand(outfile))
	invisible(list("geneset"=rr$genesets, "description"=rr$geneset.descriptions))
}

## End
