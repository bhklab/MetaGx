########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## gene set enrichment analysis
##
## inputs:	
##      - 
##			- 
##
## outputs, a list of items:
##			- 
##
## Notes:	duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
#################################################

`prerankGSEA` <- 
function (exe.path, gmt.path, rank.path, chip.path, gsea.collapse=FALSE, nperm=1000, scoring.scheme=c("weighted", "weighted_p2", "weighted_p1.5", "classic"), make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=500, set.min=15, zip.report=FALSE, gsea.report, gsea.out, replace.res=FALSE, gsea.seed=987654321) {
	exe.path <- path.expand(exe.path)
	gmt.path <- path.expand(gmt.path)
	rank.path <- path.expand(rank.path)
	gsea.seed <- as.integer(gsea.seed)
	nperm <- as.integer(nperm)
	plot.top.x <- as.integer(plot.top.x)
	set.max <- as.integer(set.max)
	set.min <- as.integer(set.min)
	if(missing(gsea.out)) { gsea.out <- "." }
	if(missing(chip.path)) { chip.path <- "" } else { chip.path <- path.expand(chip.path) }
	if(!gsea.collapse) { gsea.collapse <- "false" } else { gsea.collapse <- "true" }
	if(!make.sets) { make.sets <- "false" } else { make.sets <- "true" }
	if(!include.only.symbols) { include.only.symbols <- "false" } else { include.only.symbols <- "true" }
	if(!zip.report) { zip.report <- "false" } else { zip.report <- "true" }
	if(missing(gsea.report)) { gsea.report <- paste("gsea_report", gsub("[.]", "_", gsub("[.]rnk", "", basename(rank.path))), sep="_") }
	scoring.scheme <- match.arg(scoring.scheme)
	rest <- dir(gsea.out)
	rest <- rest[grep(pattern=sprintf("^%s", gsea.report), x=rest)[1]]
	if(!replace.res && (length(rest) > 0 && !is.na(rest))) { 
	   #warning("output directory already exists!")
	   dirn <- file.path(gsea.out, rest)
     dirfn <- dir(dirn)
     if(length(grep("gsea_report_for_na_pos_", dirfn)) == 0 || length(grep("gsea_report_for_na_pos_", dirfn)) == 0) { stop("Output directory exists but results are not available!") }
	   tt <- rbind(read.csv(file.path(gsea.out, rest, dirfn[grep("gsea_report_for_na_pos_", dirfn)[2]]), stringsAsFactors=FALSE, sep="\t", header=TRUE), read.csv(file.path(gsea.out, rest, dirfn[grep("gsea_report_for_na_neg_", dirfn)[2]]), stringsAsFactors=FALSE, sep="\t", header=TRUE))
		rownames(tt) <- as.character(tt[ ,"NAME"])
	} else {
    outlog <- file.path(gsea.out, basename(tempfile(pattern="output_log_", tmpdir="", fileext=".txt")))
		gsea.cmd <- sprintf("java -Xmx10240m -cp %s xtools.gsea.GseaPreranked -gmx %s -chip %s -collapse %s -nperm %i -rnk %s -scoring_scheme %s -rpt_label %s -include_only_symbols %s -make_sets %s -plot_top_x %i -rnd_seed %i -set_max %i -set_min %i -zip_report %s -out %s -gui false &> %s", exe.path, gmt.path, chip.path, gsea.collapse, nperm, rank.path, scoring.scheme, gsea.report, include.only.symbols, make.sets, plot.top.x, gsea.seed, set.max, set.min, zip.report, gsea.out, outlog)
		system(gsea.cmd)
		## read results
		rest <- dir(gsea.out)
		rest <- rest[grep(pattern=gsea.report, x=rest)[1]]
		restn <- sapply(strsplit(rest, "[.]"), function(x) { return(x[length(x)]) })
		tt <- rbind(read.csv(file.path(gsea.out, rest, sprintf("gsea_report_for_na_pos_%s.xls",restn)), stringsAsFactors=FALSE, sep="\t", header=TRUE), read.csv(file.path(gsea.out, rest, sprintf("gsea_report_for_na_neg_%s.xls",restn)), stringsAsFactors=FALSE, sep="\t", header=TRUE))
		rownames(tt) <- as.character(tt[ ,"NAME"])
		## rename results directory
    dirn <- file.path(gsea.out, gsub(sprintf("[.]GseaPreranked[.]%s", restn), "", rest))
		file.rename(from=file.path(gsea.out, rest), to=dirn)
    file.rename(from=outlog, to=file.path(dirn, "output_log.txt"))
	}
  ## TODO: links should be conserved when transformed into excel
  ## =Hyperlink("http://www.techonthenet.com","Tech on the Net") for instance
  
  ## for each geneset, extract the geens in the leading edge
  core.genes <- lapply(rownames(tt), function(x, dirn) {
    xx <- NULL
    if(file.exists(file.path(dirn, sprintf("%s.xls", x)))) {
      leg <- read.csv(file.path(dirn, sprintf("%s.xls", x)), stringsAsFactors=FALSE, sep="\t", header=TRUE)
      xx <- leg[complete.cases(leg[ ,c("PROBE", "CORE.ENRICHMENT")]) & leg[ ,"CORE.ENRICHMENT"] == "Yes", "PROBE"]
    }
    return(xx)
  }, dirn=dirn)
  names(core.genes) <- rownames(tt)
	return (list("geneset.results"=tt, "geneset.core"=core.genes))
}
