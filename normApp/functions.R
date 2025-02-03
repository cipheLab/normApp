
inverseArcsinhCIPHE <- function(flow.frame, marker = NULL, args) {
  raw <- flow.frame@exprs
  if (is.null(marker) || length(marker) < 1) {
    marker <- colnames(flow.frame)
  }
  for (i in seq_along(marker)) {
    mat <- raw[, marker[i]]
    mat <- sinh(mat) * args[i]
    raw[, marker[i]] <- mat
  }
  flow.frame@exprs <- raw
  return(flow.frame)
}

arcsinhCIPHE <- function(flow.frame, marker = NULL, args) {
  
  raw <- flow.frame@exprs
  if (is.null(marker) || length(marker) < 1) {
    marker <- colnames(flow.frame)
  }
  for (i in seq_along(marker)) {
    mat <- raw[, marker[i]]
    mat <- asinh(mat / args[i])
    raw[, marker[i]] <- mat
  }
  flow.frame@exprs <- raw
  return(flow.frame)
}
