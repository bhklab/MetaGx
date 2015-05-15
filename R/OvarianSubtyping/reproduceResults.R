library(gdata)
library(Biobase)
library(GSVA)

load("test/eset.scaled.rda")
source("getVerhaakSubtypes.R")

### Verhaak et al., 2013

spreadsheet.data <- read.xls("input/JCI65833sd1.xls", skip=1)

## Examine the correlation between normalized ssGSEA scores in the Tothill dataset, as predicted by our implementation and that of Verhaak et al.
tothill.eset <- esets$GSE9891_eset
implemented.verhaak.output <- getVerhaakSubtypes(tothill.eset)

# extract the official CLOVAR subtype scores and compare with our implementation
tothill.spreadsheet.data <- spreadsheet.data[grep("Tothill", spreadsheet.data$DATASET),]

# Include patients in both data sets
tothill.spreadsheet.data <- tothill.spreadsheet.data[c(1,grep(".ssGSEA.normalized.score", colnames(tothill.spreadsheet.data)))]
rownames(tothill.spreadsheet.data) <- tothill.spreadsheet.data$ID
tothill.spreadsheet.data <- tothill.spreadsheet.data[,-1]
merged <- merge(implemented.verhaak.output$gsva.out, tothill.spreadsheet.data, by="row.names")

# Plot correlations between ssGSEA scores
plot(merged$DIF, merged$Differentiated.ssGSEA.normalized.score)
plot(merged$IMR, merged$Immunoreactive.ssGSEA.normalized.score)
plot(merged$MES, merged$Mesenchymal.ssGSEA.normalized.score)
plot(merged$PRO, merged$Proliferative.ssGSEA.normalized.score)