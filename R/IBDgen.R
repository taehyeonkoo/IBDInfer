#' @title Generating Incomplete Block Designs
#' @description Generate incomplete block designs.
#'
#'
#' @param K The number of blocks.
#' @param n.trt The number of whole treatments.
#' @param t The number of treatments to be assigned to each block.
#' @param n.vec The vector of block sizes.
#' @param L The vector of the number of blocks having each treatment.
#' @param l The matrix of the number of blocks having each pair of treatments.
#' @param W The set of treatment subsets used in the design.
#' @param balanced Whether the design is balanced or not. If \code{TRUE}, generate a balanced design.
#'
#' @returns A list containing the following components:
#' \item{W}{The set of treatment subsets used in the design.}
#' \item{W.uniq}{The unique set of treatment subsets used in the design with proportion in \code{W}.}
#' \item{Rk}{The block assignment matrix.}
#' \item{blk_assign}{The block assignment data frame.}
#'
#' @importFrom crossdes find.BIB isGYD
#' @export
#'
#' @examples
#' K <- 6
#' n.trt <- 3
#' t <- 2
#' n.vec <- rep(4, K)
#' IBDgen(K = K, n.trt = n.trt, t = t, n.vec = n.vec)
#'
#' @references {
#' Sailer, M. O., & Bornkamp, M. B. (2022). Package ‘crossdes’: Construction of Crossover Designs.
#' }
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

  # Sort each row so permutations match
  sorted_mat <- t(apply(W, 1, sort))

  # Create string keys for row comparison
  row_keys <- apply(sorted_mat, 1, paste, collapse = ",")

  # Count frequencies and proportions
  freq_table <- table(row_keys)
  proportions <- as.numeric(freq_table) / nrow(W)

  # Convert row_keys back to numeric matrix
  unique_rows <- do.call(rbind, strsplit(names(freq_table), ","))
  unique_rows <- apply(unique_rows, 2, as.numeric)

  # Build data frame
  result_df <- as.data.frame(unique_rows)

  # Add proportion column
  result_df$prop <- proportions

  W.uniq <- result_df

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


  return(list(W = W, W.uniq = W.uniq, Rk = Rk, blk_assign = df))
}
