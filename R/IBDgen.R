IBDgen <- function(K,n.trt,t,n.vec = NULL, L = NULL, l = NULL,W = NULL,balanced = T){
  stopifnot(!is.null(K),!is.null(n.trt),!is.null(t))
  stopifnot(n.trt>t)
  if (balanced) {
    if (is.null(L)) {
      L <-  K*t/n.trt
    }
    if (is.null(l)) {
      l <- L*(t-1)/(n.trt-1)
    }
    stopifnot(round(K*t/n.trt)==(K*t/n.trt) & round( L*(t-1)/(n.trt-1))==(L*(t-1)/(n.trt-1)))

    W <- crossdes::find.BIB(trt = n.trt, b = K, k = t) # W = \{R_k\}_{k=1}^K
    crossdes::isGYD(W)
  } else{
    stopifnot(!is.null(W))
  }

  stopifnot(is.matrix(W),K %% dim(W)[1] == 0, t == dim(W)[2])
  if (is.null(n.vec)){
    n.vec <- rep(t,K) # Block size
  }
  stopifnot(length(n.vec) == K & all(n.vec %% t ==0))

  # Enlarging W if nrow(W)<K #
  if (dim(W)[1] < K){
    Rk <-  kronecker(matrix(1, K/(dim(W)[1]), 1), W)
  } else {
    Rk <- W
  }

  # First stage randomization #
  Rk <- Rk[sample(1:K),]

  # Second stage randomization #
  # Initialize an empty list to store the results
  df_list <- list()

  # Loop over blocks
  for (k in 1:K) {
    Zi <- sample(rep(Rk[k,], each = n.vec[k] / t))  # Randomly sample treatment assignments
    df_list[[k]] <- data.frame(blk_id = rep(k, n.vec[k]), assign = Zi)  # Store each block's data
  }

  # Combine the list into a single data frame
  df <- do.call(rbind, df_list)


  return(list(W = W, Rk = Rk, blk_assign = df))
}
