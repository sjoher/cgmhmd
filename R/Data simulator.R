data_sim <- function(network, n, p, K, ncat, rho, gamma_g = NULL, gamma_o, gamma_b = NULL, gamma_p = NULL, prob = NULL, nclass = NULL){
  if(network == "Circle"){
    theta_base <- matrix(0, p, p)
    diag(theta_base) <- 1
    diag(theta_base[-nrow(theta_base),-1]) <- 0.5
    diag(theta_base[-1,-ncol(theta_base)]) <- 0.5
    theta_base[ 1, p ] <- 0.4
    theta_base[ p, 1 ] <- 0.4
  }

  if(network == "Random"){
    theta_base <- matrix(0, p, p)
    theta_base[upper.tri(theta_base)] <- rbinom(p * (p - 1)/2, 1, prob)
    theta_base <- theta_base + t(theta_base)
    diag(theta_base) <- 1
  }

  if(network == "Cluster"){
    theta_base <- matrix(0, p, p)
    g.large <- p%%nclass
    g.small <- nclass - g.large
    n.small <- floor(p/nclass)
    n.large <- n.small + 1
    vp <- c(rep(n.small, g.small), rep(n.large, g.large))
      if (length(prob) != nclass)
        prob = rep(prob, nclass)
      for (i in 1:nclass) {
        tmp <- if (i == 1)
          (1:vp[1])
        else ((sum(vp[1:(i - 1)]) + 1):sum(vp[1:i]))
        gg <- matrix(0, vp[i], vp[i])
        gg[upper.tri(gg)] <- rbinom(vp[i] * (vp[i] - 1)/2, 1, prob[i])
        theta_base[tmp, tmp] <- gg
      }
    nsamp <- length(theta_base[which(theta_base == 1)])
    theta_base[which(theta_base == 1)] <- sample(c(runif(nsamp, -1, -0.5),runif(nsamp, 0.5, 1)), nsamp)
    theta_base[lower.tri(theta_base, diag = TRUE)] <- 0
    theta_base <- theta_base + t(theta_base)
    diag(theta_base) <- 1
  }

  if(network == "Scale-free"){
    theta_base <- matrix(as.numeric(graph.sim(p, graph = "scale-free")),nrow = p, ncol = p)
    nsamp <- length(theta_base[which(theta_base == 1)])
    theta_base[which(theta_base == 1)] <- sample(c(runif(nsamp, -1, -0.5),runif(nsamp, 0.5, 1)), nsamp)
    theta_base[lower.tri(theta_base, diag = TRUE)] <- 0
    theta_base <- theta_base + t(theta_base)
    diag(theta_base) <- 1
  }

  if(network == "AR1")
  {
    sigma = matrix( 0, p, p )

    for( i in 1 : ( p - 1 ) )
      for( j in ( i + 1 ) : p )
        sigma[ i, j ] = ( 0.7 ) ^ abs( i - j )

    sigma = sigma + t( sigma ) + diag( p )
    theta_base     = solve( sigma )
  }

  if(network == "AR2")
  {
    theta_base = stats::toeplitz( c( 1, 0.5, 0.25, rep( 0, p - 3 ) ) )
  }

  Theta <- vector(mode = "list", K)
  for(k in 1:K){
    Theta[[k]] <- theta_base
  }

    # Make theta different
    theta_m <- Theta[[1]]
    theta_m[upper.tri(theta_m, diag = TRUE)] <- NA
    if(network == "AR1"){
      M <- length(which(abs(theta_m) > 0.0001))
      zeros <- which(abs(theta_m) < 0.0001 ,arr.ind = T)
    } else{
      M <- length(which(theta_m != 0))
      zeros <- which(theta_m ==0,arr.ind = T)
    }

    # If zeros > M, we randomly select rho*M elements and set those to nonzero to make precision matrices different
    if(nrow(zeros) > M){
      zero_samp <- vector("list", K)
      for(k in 1:K){
        zero_samp[[k]] <- zeros[sample(nrow(zeros), floor(rho*M)), ]
        for(i in 1:floor((rho*M))){
          if(floor((rho*M)) == 1){
            Theta[[k]][as.numeric(zero_samp[[k]][1]),as.numeric(zero_samp[[k]][2])] <- rep(0.5, rho*M)[i]
          } else{
            Theta[[k]][zero_samp[[k]][i,1],zero_samp[[k]][i,2]] <- rep(0.5, rho*M)[i]
          }
        }
        Theta[[k]][upper.tri(Theta[[k]], diag = FALSE)] <- NA
        Theta[[k]][upper.tri(Theta[[k]])] <- t(Theta[[k]])[upper.tri(Theta[[k]])]
      } # Otherwise, we randomly select (M*rho)/round(sqrt(p)) nonzero edges and set those to 0
    } else{
        non_zero_samp <- vector("list", K)
        nonzeros <- which(abs(theta_m) > 0.0001 ,arr.ind = T)
        for(k in 1:K){
          non_zero_samp[[k]] <- nonzeros[sample(nrow(nonzeros), (M*rho)/round(sqrt(p))), ]
          for(i in 1:((M*rho)/round(sqrt(p)))){
            Theta[[k]][non_zero_samp[[k]][i,1],non_zero_samp[[k]][i,2]] <- rep(0.5, (M*rho)/round(sqrt(p)))[i]
          }
          Theta[[k]][upper.tri(Theta[[k]], diag = FALSE)] <- NA
          Theta[[k]][upper.tri(Theta[[k]])] <- t(Theta[[k]])[upper.tri(Theta[[k]])]
        }
    }

    Sigma <- vector("list", K)
    #Generate latent correlation matrix
    for(k in 1:K){
      diag(Theta[[k]]) <- 0
      Theta[[k]] <- Theta[[k]]*0.5
      diag(Theta[[k]]) <- abs(min(eigen(Theta[[k]])$values)) + 0.5
      Sigma[[k]] <- cov2cor(solve(Theta[[k]]))
      Theta[[k]] <- solve(Sigma[[k]])
    }
    # Simulate N(0,Sigma) latent variables Z
    z <- vector("list", K)
    for(k in 1:K){
      z[[k]] <- mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma[[k]])
    }

    n_ordinal <- floor(gamma_o*p)
    n_binary <- floor(gamma_b*p)
    n_pois <- floor(gamma_p*p)
    ord_cols <- sample(1:p, n_ordinal)
    bin_cols <- sample((1:p)[-ord_cols], n_binary)
    pois_cols <- sample((1:p)[-c(ord_cols,bin_cols)], n_pois)

    # Simulate ordinal data
    if(gamma_o > 0){
      cutoff   <- matrix( runif( ncat * n_ordinal ), nrow = n_ordinal, ncol = ncat )
      marginals <- apply( cutoff, 1, function(x) { qnorm( cumsum( x / sum(x) )[-length(x)] ) } )
      for(k in 1:K){
        for(j in 1:n_ordinal)
        {
          breaks <- c( min( z[[k]][,ord_cols[j]] ) - 1, marginals[,j], max(z[[k]][,ord_cols[j]]) + 1 )
          z[[k]][,ord_cols[j]]  <- as.integer( cut( z[[k]][,ord_cols[j]], breaks = breaks, right = FALSE ) )
        }
      }
    }

    # Simulate binary data
    if(gamma_b > 0){
      for(k in 1:K){
        for(j in 1:n_binary){
          bin_prob <- pnorm(z[[k]][, bin_cols[j]])
          z[[k]][ ,bin_cols[j]] <- qbinom(p = bin_prob, size = 1, prob = 0.5 )
        }
      }
    }

    # Simulate poisson data
    if(gamma_p > 0){
      for(k in 1:K){
        for(j in 1:n_pois){
          pois_prob <- pnorm(z[[k]][, pois_cols[j]])
          z[[k]][ ,pois_cols[j]] <- qpois(p = pois_prob, lambda = 10)
        }
      }
    }

  return(list(z = z, theta = Theta))
}

