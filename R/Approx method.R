#Approx method 2

initialize = function(y, ncores)
{
  p = ncol(y[[1]])
  K <- length(y)
  n_k <- numeric(K)
  lower_upper <- vector("list", K)

  for(k in 1:K){
    lower_upper[[k]] <- lower.upper(y[[k]])
    n_k[k] <- nrow(y[[k]])
  }

  ## Initialize S
  Z <- vector("list", K)
  diag_element <- vector("list", K)
  for(k in 1:K){
    Z[[k]] <- matrix(0, n_k[k], p)
    diag_element[[k]] <- rep(0, p)
  }


  tmp2 <- vector("list", K)
  ES <- vector("list", K)
  for(k in 1:K){
    if(ncores > 1){
      cl <- makeCluster(ncores)
      clusterExport(cl, varlist = c("element_S_j", "lower_upper"), envir=environment())
      tmp2[[k]] <- parLapply(cl = cl, 1:p, function(i){ element_S_j(i, lower_upper[[k]] );})
      stopCluster(cl)
    } else{
      tmp2[[k]] <- lapply(1:p, function(i){ element_S_j(i, lower_upper[[k]] );})
    }
    # So for every p for every k, we get one EX and one EXX
    Z[[k]] <-  do.call(cbind, lapply(1:p, function(x) tmp2[[k]][[x]]$EX ))
    diag_element[[k]] <- unlist(lapply(1:p, function(x)mean(tmp2[[k]][[x]]$EXX)))
    ES[[k]] <- t(Z[[k]]) %*% Z[[k]] / n_k[k]
    diag(ES[[k]]) <- diag_element[[k]]
  }

  return(list(Z=Z, ES = ES, lower_upper=lower_upper))
}

cutoffs = function(y){
  p<-ncol(y)
  n<-nrow(y)
  k<-unique(sort(unlist(y)))
  n.levels<-length(k)
  q<-matrix(nrow=p,ncol=n.levels)
  for(i in 1:p){
    X=factor(y[,i],levels=k)
    No<-tabulate(X, nbins=n.levels)
    q[i,]<-qnorm(cumsum(No)/n)
  }
  q[ ,n.levels] <- Inf
  q<-cbind(-Inf,q)
  return(q)
}

lower.upper = function(y){
  cutoffs <- cutoffs(y)
  levels	<- unique(sort(unlist(y)))
  n.levels<- length(levels)
  n <- nrow(y)
  p <- ncol(y)
  lower = matrix(nrow=n,ncol=p)
  upper = matrix(nrow=n,ncol=p)
  for (i in 1:n){
    sel <- match(y[i,],levels)
    lower[i,]<-apply(cbind(sel,cutoffs),1,function(x){x[x[1]+1]})
    upper[i,]<-apply(cbind(sel,cutoffs),1,function(x){x[x[1]+2]})
  }
  lower[is.na(lower)] <- -Inf
  upper[is.na(upper)] <- Inf

  return(list(lower=lower,upper=upper))
}

element_S_j = function(j, lower_upper, mu=0, sigma=1)
{
  delta1 <- (lower_upper$lower[ ,j] - mu) / sigma
  delta2 <- (lower_upper$upper[ ,j] - mu) / sigma
  tmp1 <- (dnorm(delta1) - dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
  EX <- mu + tmp1 * sigma

  delta1[delta1 < -1e+10] <- -1e+10
  delta2[delta2 > 1e+10] <- 1e+10
  tmp2 <- (delta1*dnorm(delta1) - delta2*dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
  EXX <- sigma^2 + mu^2 + sigma^2 * tmp2 + 2 * mu * sigma * tmp1
  rm(delta1, delta2, tmp1, tmp2)
  gc()

  return(list(EX=EX, EXX=EXX))
}




Approx_method <- function(y, Z, ES=NULL, lambdas, lower_upper=NULL, combination = 1, em_tol=.001, em_iter=10, ncores){

  element_S <- function( lower, upper, mu=0, sigma=1)
  {
    delta1 <- (lower - mu) / sigma
    delta2 <- (upper - mu) / sigma
    tmp1 <- (dnorm(delta1) - dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
    EX <- mu + tmp1 * sigma

    delta1[delta1 < -1e+10] <- -1e+10
    delta2[delta2 > 1e+10] <- 1e+10
    tmp2 <- (delta1*dnorm(delta1) - delta2*dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
    EXX <- sigma^2 + mu^2 + sigma^2 * tmp2 + 2*mu*sigma*tmp1
    rm(delta1, delta2, tmp1, tmp2 )
    gc()

    return(list(EX=EX, EXX=EXX))
  }

  cond_N <- function(j, Sigma, Z , Z_new, diag_element, lower_upper)
  {
    p <- ncol(Sigma)
    tmp <- matrix(Sigma[j, -j], 1, p-1)
    tmp1 <- solve(Sigma[-j, -j])

    mu <- tmp %*% tmp1 %*% t(Z[, -j])
    mu <- as.vector(mu)
    sigma <- Sigma[j, j] - tmp %*% tmp1 %*% t(tmp)
    sigma <- sqrt(sigma)

    obj <- element_S( lower= lower_upper$lower[ ,j], upper= lower_upper$upper[ ,j], mu=mu, sigma=sigma)

    Z_new <- obj$EX
    diag_element <- mean(obj$EXX)

    rm(tmp, tmp1, mu, sigma, obj )
    gc()

    return(list(Z_new=Z_new, diag_element=diag_element))
  }

  calculate.R.approx.internal <- function(y, Z, lower_upper, Sigma, ncores = NULL)
  {
    if(is.null(ncores)) ncores = detectCores() - 1
    p <- ncol(y)
    n <- nrow(y)
    Z_new <- rep(0, n)

    if(ncores > 1){
      cl <- makeCluster(ncores)
      clusterExport(cl, varlist = c("cond_N", "lower_upper"), envir=environment())
      Sigma <- Sigma
      Z <- Z
      cond_norm <- parLapply(cl = cl, 1:p, function(j) {
        cond_N(j, Sigma=Sigma, Z=Z , Z_new=Z_new, diag_element=diag_element, lower_upper=lower_upper );
      })
      stopCluster(cl)
    }else{
      cond_norm <- lapply(1:p, function(j){ cond_N(j, Sigma, Z , Z_new, diag_element, lower_upper); } )
    }


    Z_new <- do.call(cbind, lapply(1:p, function(x) cond_norm[[x]]$Z_new ))
    diag_element <- sapply(1:p, function(x) cond_norm[[x]]$diag_element)

    ES <- t(Z_new) %*% Z_new / n
    diag(ES) <- diag_element

    output <- list()
    output$ES <- ES
    output$Z <- Z_new

    rm(lower_upper, ES, Z_new )
    return(output)
  }

  setup_approx = function(y, Z, ES=NULL, lambdas, lower_upper=NULL, combination = 1, em_tol=.001, em_iter=10, ncores)
  {
    K <- length(y)
    if(class(lambdas) == "numeric")
    {
      obj <- fused_gl_copula(S = ES, Y = y, lambda1 = lambdas[1], lambda2 = lambdas[2], maxiter=1000, penalize.diagonal=FALSE)
    }
    else
    {
      obj <- fused_gl_copula(S = ES, Y = y, lambda1 = lambdas[combination,1], lambda2 = lambdas[combination,2], maxiter=1000, penalize.diagonal=FALSE)
    }
    Theta <- vector("list", K)
    Sigma <- vector("list", K)
    sd_marginal <- vector("list", K)
    for(k in 1:K){
      Theta[[k]] <- Matrix((t(obj$theta[[k]]) + obj$theta[[k]]) / 2, sparse = TRUE)
      Sigma[[k]] <- (t(solve(obj$theta[[k]])) + solve(obj$theta[[k]])) / 2
      sd_marginal[[k]] <- sqrt(diag(Sigma[[k]]))
      sd_marginal[[k]][abs(sd_marginal[[k]]) < 1e-10] <- 1e-10
      Sigma[[k]] <- diag(1/sd_marginal[[k]]) %*% Sigma[[k]] %*% diag(1/sd_marginal[[k]])
      Theta[[k]] <- diag(sd_marginal[[k]]) %*% Theta[[k]] %*% diag(sd_marginal[[k]])
    }
    if(class(lambdas) == "numeric")
    {
      approx.method <- calculate_EM_approx(combination = combination, y=y, Z=Z, lambda1 = lambdas[1], lambda2 = lambdas[2], Sigma=Sigma, Theta= Theta, lower_upper=lower_upper, em_tol=em_tol, em_iter = em_iter, ncores = ncores)
    }
    else
    {
      approx.method <- calculate_EM_approx(combination = combination, y=y, Z=Z, lambda1 = lambdas[combination,1], lambda2 = lambdas[combination,2], Sigma=Sigma, Theta= Theta, lower_upper=lower_upper, em_tol=em_tol, em_iter = em_iter, ncores = ncores)
    }


    invisible(return(approx.method))
  }

  calculate_EM_approx = function(combination, y, Z, lambda1, lambda2, Sigma, Theta, lower_upper, em_tol, em_iter, verbose=FALSE, ncores)
  {
    c_em_iter = 1
    p = ncol(y[[1]])
    K <- length(y)
    difsum <- K*100
    S_obj <- vector(mode = "list", length = K)
    S_obj_ES <- vector(mode = "list", length = K)
    Sigma_new <- vector(mode = "list", length = K)
    sd_marginal <- vector(mode = "list", length = K)
    dif <- numeric(K)
    n_k <- numeric(K)
    for(k in 1:K){
      n_k[k] <- nrow(y[[k]])
    }

    while((c_em_iter <= em_iter) && (difsum >= (K*em_tol)))
    {
      for(k in 1:K){
        S_obj[[k]]	<- calculate.R.approx.internal(y=y[[k]], Z=Z[[k]], lower_upper=lower_upper[[k]], Sigma= Sigma[[k]], ncores = ncores)
        Z[[k]]		<- S_obj[[k]]$Z
        S_obj_ES[[k]] <- S_obj[[k]]$ES
      }
      S_fgl	<- fused_gl_copula(S = S_obj_ES, Y = y, lambda1 = lambda1, lambda2 = lambda2, maxit=1000, penalize.diagonal=FALSE)
      for(k in 1:K){
        Theta[[k]] <- (t(S_fgl$theta[[k]]) + S_fgl$theta[[k]]) / 2
        Sigma_new[[k]] <- (t(solve(S_fgl$theta[[k]])) + solve(S_fgl$theta[[k]])) / 2
        sd_marginal[[k]] <- sqrt(diag(Sigma_new[[k]]))
        sd_marginal[abs(sd_marginal[[k]]) < 1e-10][[k]] <- 1e-10
        Sigma_new[[k]]	<- diag(1/sd_marginal[[k]]) %*% Sigma_new[[k]] %*% diag(1/sd_marginal[[k]])
        Theta[[k]]	<- matrix(diag(sd_marginal[[k]]) %*% Theta[[k]] %*% diag(sd_marginal[[k]]), ncol = p)
        dif[k]		<-	sum(abs(Sigma_new[[k]] - Sigma[[k]])/p^2)
        Sigma[[k]]	<-  Sigma_new[[k]]
      }
      difsum <- sum(dif)
      c_em_iter <- c_em_iter + 1
    }

    loglik <- numeric(K)
    for(k in 1:K){
      #loglik[k] <- n_k[k]/2 *(log(det(Theta[[k]])) - sum(diag(S_obj_ES[[k]] %*%Theta[[k]])))
      loglik[k] <- n_k[k]/2 *(determinant(Theta[[k]], logarithm = TRUE)$modulus - sum(diag(S_obj_ES[[k]] %*%Theta[[k]])))
    }
    loglik <- sum(loglik)

    results <- list(Z=Z, ES = S_obj_ES, Sigma=Sigma, Theta = Theta, loglik = loglik)

    return(results)
  }
  Approx_method <- setup_approx(y, Z, ES=ES, lambdas=lambdas, lower_upper=lower_upper, combination = combination, em_tol=0.001, em_iter=5, ncores = ncores)
  return(Approx_method)
}
