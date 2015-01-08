########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Define the GO terms using biomaRt and create a GMT file for GSEA
##
## inputs:	
##      - gid is a vector of gene ids
##      - value is a string specifying what are the gene ids (for example "ensembl_gene_id" or "ensembl_transcript_id")
##      - gokeep is the catergory to keep.  Only one GO category
##      - mart.db is mart database such as mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
##      - outfile is the path to write the output file
## verbose=TRUE will display some messages
##      - ... are additional parameters to be passed to createGMT
##
## outputs, a list of 2 items:
##			- 
##
## Notes:	duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
#################################################

`createGMT` <-
function (gid, value, gokeep=c("biological_process", "molecular_function", "cellular_component"), mart.db, outfile, replace=FALSE, verbose=TRUE, ...) {

	## create GO terms for ensembl ids
	gokeep <- match.arg(gokeep)
	if(missing(gid)) { gid <- getBM(attributes=value, filters="", values="", mart.db, ...)[ ,1] }
  
  if(missing(outfile)) { outfile <- "GO_TERM.gmt" }
  outdir <- path.expand(dirname(outfile))
  if(!file.exists(outdir)) { dir.create(outdir, recursive=TRUE, showWarnings=FALSE) }
  if(file.exists(outfile)) {
    if(!replace) { stop("Output file already exists!") }
    file.remove(outfile)
  }
  
  badchars <- "[,][:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
	#for all genes, extract GO terms
	gene.an <- getBM(attributes=c(value, "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003"), filters=value, values=gid, mart=mart.db, ...)
	gene.an[gene.an == "" | gene.an == " "] <- NA
	gene.an <- gene.an[!is.na(gene.an[ ,"namespace_1003"]) & is.element(gene.an[ ,"namespace_1003"], gokeep), ,drop=FALSE]
	gene.an <- data.frame(gene.an, "GONAME"=gsub(pattern=badchars, replacement="_", x=toupper(gene.an[ , "name_1006"])))

	goo <- sort(unique(gene.an[ ,"go_id"]))
	names(goo) <- as.character(gene.an[match(goo, gene.an[ ,"go_id"]), "GONAME"])
	goo2 <- cbind(names(goo), goo)
	rownames(goo2) <- names(goo)
	golist <- apply(goo2, 1, function(x, z) {
	res <- c(x[1], x[2], unique(z[is.element(z[ ,"go_id"], x[2]), value]))
	names(res) <- NULL
	return(res)
	}, z=gene.an)
	names(golist) <- rownames(goo2)

	## write gmt file
	if(verbose) { message(sprintf("writing %s to %s", outfile, outdir)) }
	rr <- lapply(golist, function(x, file) { write(sprintf("%s\thttp://www.ebi.ac.uk/QuickGO/GTerm?id=%s\t%s", x[1], x[2], paste(unique(x[3:length(x)]), collapse="\t")), file=file, append=TRUE) }, file=outfile)
	invisible(golist)
}

## End
