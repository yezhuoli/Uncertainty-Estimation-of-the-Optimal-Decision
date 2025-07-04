#############################################################
##                          Scale                          ##
#############################################################

# Scale to unit
scale_to_uni <- function(xF, lb, ub) {
  # If x is a matrix, apply the operation to each column
  if (is.matrix(xF) || is.data.frame(xF)) {
    t(apply(xF, 1, function(col) (col - lb) / (ub - lb) ))
  } else {
    # If x is not a matrix (i.e., a vector/array), scale directly
    (xF - lb) / (ub - lb)
  }
}

# Scale to original
scale_to_org <- function(xF_uni, lb, ub) {
  # If x is a matrix, apply the operation to each column
  if (is.matrix(xF_uni) || is.data.frame(xF_uni)) {
    t(apply(xF_uni, 1, function(col) col * (ub - lb) + lb ))
  } else {
    # If x is not a matrix (i.e., a vector/array), scale directly
    lb + xF_uni * (ub - lb)
  }
}

#############################################################