########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

## run this code from the working directory
## source(file.path("R", "bcmeta_isdb_data"))

## remove all existing objects from the workspace
rm(list=ls(all=TRUE))

require(InSilicoDb)
require(genefu)

## this version of jetset allow for simple probe-gene mapping for affymetrix chips hgu95av2 (GPL8300), hgu133a (GPL96), hgu133b (GPL97), hgu133plus2 (GPL570), u133x3p (GPL1352)

# source(file.path("R", "bcmeta_foo.R"))
