modselect <- function(est, X, l1, l2){

  K <- length(X)
  n_lambda <- nrow(expand.grid(l1, l2))
  thetas <- vector("list", n_lambda)
  A <- vector("list", n_lambda)

  K <- length(X)
  p <- ncol(X[[1]])
  n <- sapply(X, nrow)

  nu <- matrix(NA, n_lambda, K)
  lik_mat <- matrix(NA, n_lambda, K)
  lik <- numeric(n_lambda)
  aic <- numeric(n_lambda)
  ebic <- numeric(n_lambda)
  gamma <- 0.5

  lognu <- matrix(nrow = n_lambda, ncol = K)
  rhs <- matrix(nrow = n_lambda, ncol = K)

  for(i in 1:n_lambda)
  {
    for(k in 1:K){
      thetas[[i]][[k]] <- as.matrix(est[[i]]$Theta[[k]])
      A[[i]][[k]] <- thetas[[i]][[k]]
      # adjacency matrix
      A[[i]][[k]][which(A[[i]][[k]] != 0)] <- 1
      # nonzero elements penalty
      nu[i,k] <- sum(A[[i]][[k]][upper.tri(A[[i]][[k]])])
      # likelihood
      lik_mat[i,k] <- n[k]/2 *(determinant(as.matrix(thetas[[i]][[k]]), logarithm = TRUE)$modulus - sum(diag(as.matrix(est[[1]]$ES[[k]]) %*%as.matrix(thetas[[i]][[k]]))))

      # IC penalties
      lognu[i,k] <- log(n[k]) * nu[i,k]
      rhs[i,k] <- (4 * gamma * log(p) * nu[i,k])
    }

    lik[i] <- sum(lik_mat[i,])

    # information criteria
    aic[i] <- ( - 2 * lik[i] ) + ( 2 * sum(nu[i,]) )
    ebic[i] <- - 2 * lik[i] + sum(lognu[i,]) + sum(rhs[i,])
  }
  aic_idx <- which.min(aic)
  ebic_idx <- which.min(ebic)

  l1_aic <- as.numeric(expand.grid(l1, l2)[aic_idx,][1])
  l2_aic <- as.numeric(expand.grid(l1, l2)[aic_idx,][2])
  l1_ebic <- as.numeric(expand.grid(l1, l2)[ebic_idx,][1])
  l2_ebic <- as.numeric(expand.grid(l1, l2)[ebic_idx,][2])

  returnlist <- list("aic_idx", "ebic_idx", "l1_aic" = l1_aic, "l2_aic" = l2_aic, "l1_ebic" = l1_ebic, "l2_ebic" = l2_ebic)
  return(returnlist)
}
