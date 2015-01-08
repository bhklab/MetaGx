########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

`blueYellow` <-
function (n) {
  x <- (1:n)-1
  rgb(x, x, rev(x), maxColorValue=n)
}