########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################


`intersectList` <-
function(...) {
   args <- list(...)
   nargs <- length(args)
   if (nargs <= 1) {
     if (nargs == 1 && is.list(args[[1]])) {
       do.call("intersectList", args[[1]])
     } else {
       return (args[[1]])
     }
   } else if (nargs == 2) {
     return (intersect(args[[1]], args[[2]]))
   } else {
     return (intersect(args[[1]], intersectList(args[-1])))
   }
}

## End

