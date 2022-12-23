library(Matrix)
library(tmvtnorm)

Gibbs_method <- function(y, lambdas=lambdas, n_lambda=n_lambda, Theta = Theta, combination = combination, max.elongation = em.iter, em.tol=em.tol, ncores){

  p <- dim(y[[1]])[2]
  K <- length(y)

  cutoffs = function(Y){
    p<-ncol(Y)
    n<-nrow(Y)
    k<-unique(sort(unlist(Y)))
    n.levels<-length(k)
    q<-matrix(nrow=p,ncol=n.levels)
    for(i in 1:p){
      X=factor(Y[,i],levels=k)
      No<-tabulate(X, nbins=n.levels)
      q[i,]<-qnorm(cumsum(No)/n)
    }
    q[ ,n.levels] <- Inf
    q<-cbind(-Inf,q)
    return(q)
  }

  lower.upper = function(Y){
    cutoffs <- cutoffs(Y)
    levels	<- unique(sort(unlist(Y)))
    n.levels<- length(levels)
    n <- nrow(Y)
    p <- ncol(Y)
    lower = matrix(nrow=n,ncol=p)
    upper = matrix(nrow=n,ncol=p)
    for (i in 1:n){
      sel <- match(Y[i,],levels)
      lower[i,]<-apply(cbind(sel,cutoffs),1,function(x){x[x[1]+1]})
      upper[i,]<-apply(cbind(sel,cutoffs),1,function(x){x[x[1]+2]})
    }
    lower[is.na(lower)] <- -Inf
    upper[is.na(upper)] <- Inf
    return(list(lower=lower,upper=upper))
  }

  lower_upper <- vector("list", K)
  for(k in 1:K){
    lower_upper[[k]] <- lower.upper(y[[k]])
  }

  calcRcore <- function(i, n, p, gibbs.iter, mc.iter, theta, lower_upper, verbose = TRUE){
    z.gibbs <- tmvtnorm::rtmvnorm.sparseMatrix(n=mc.iter, H= theta , lower=lower_upper$lower[i,], upper=lower_upper$upper[i,], burn.in=gibbs.iter)
    summed <- matrix(0, p, p)
    for(iteration in 1: mc.iter)
    {
      summed <- summed + (z.gibbs[iteration,] %*% t(z.gibbs[iteration,]))
    }
    R <- summed /  mc.iter
    return(R)
  }

  calculate.R.internal = function(y, theta=NULL, lower_upper=NULL, gibbs.iter=1000, mc.iter=1000, verbose = TRUE, ncores)
  {
    if(missing(y)) stop("argument \"y\" is missing, with no default")
    S <- 0
    p <- ncol(y)
    n <- nrow(y)
    if(missing(lower_upper)) lower_upper <- lower.upper(y)

    if(ncores > 1){
      cl <- makeCluster(ncores)
      clusterExport(cl, varlist = c("calcRcore", "Theta"), envir=environment())
      scovs <- parLapply(cl = cl, 1:n, function(i) {
        calcRcore(i, n, p, gibbs.iter, mc.iter, theta, lower_upper);
      })
      stopCluster(cl)
    }else{
      scovs <- lapply(1:n, function(i){ calcRcore(i, n, p, gibbs.iter, mc.iter, theta, lower_upper); })
    }

    for(i in 1:n){
      S <- S + scovs[[i]]
    }
    ES  <- cov2cor(S/n)
    rm(S)
    return(ES)
  }

  calculate_EM_Gibbs = function(combination, y, lambda1 = lambdas[1,1], lambda2 = lambdas[1,2], Theta=NULL, lower_upper, em.tol=.001, em.iter=10, c.em.iter=1, gibbs.iter=1000, mc.iter=1000, ncores)
  {
    p = ncol(y[[1]])
    c.em.iter = 1
    K <- length(y)
    n_k <- numeric(K)
    for(k in 1:K){
      n_k[k] <- nrow(y[[k]])
    }

    difsum <- K*100
    R <- vector(mode = "list", length = K)
    Sigma <- vector(mode = "list", length = K)
    dif <- numeric(K)
    while((c.em.iter <= em.iter) && (difsum >= (K*em.tol)))
    {
      for(k in 1:K){
        R[[k]]   <- calculate.R.internal(y[[k]], theta = Theta[[k]], lower_upper = lower_upper[[k]], gibbs.iter = gibbs.iter, mc.iter = mc.iter, ncores = ncores)
      }
      R.fgl <- fused_gl_copula(S = R, Y = y, lambda1, lambda2)
      for(k in 1:K){
        Sigma[[k]] <- solve(R.fgl$theta[[k]])
        if(det(Sigma[[k]]) <= 0){
          Sigma[[k]] <- nearPD(Sigma[[k]],keepDiag=TRUE)$mat
        }
        dif[k] <- sum(abs(Theta[[k]] - R.fgl$theta[[k]])/p^2)
      }
      difsum <- sum(dif)
      for(k in 1:K){
        Theta[[k]] = as(R.fgl$theta[[k]], "dgTMatrix")
        Theta[[k]] = as(Theta[[k]], "sparseMatrix")
      }
      c.em.iter <- c.em.iter + 1
    }

    results <- list()
    thetalist <- vector("list", K)
    sigmalist <- vector("list", K)
    eslist <- vector("list", K)
    for(k in 1:K){
      thetalist[[k]] <- Matrix(Theta[[k]], sparse=TRUE)
      sigmalist[[k]] <- Matrix(Sigma[[k]], sparse=TRUE)
      eslist[[k]] <- Matrix(R[[k]], sparse=TRUE)
    }

    loglik <- numeric(K)
    for(k in 1:K){
      loglik[k] <- n_k[k]/2 *(determinant(Theta[[k]], logarithm = TRUE)$modulus - sum(diag(R[[k]] %*%Theta[[k]])))
    }
    loglik <- sum(loglik)


    results$Theta	<- thetalist
    results$Sigma	<- sigmalist
    results$ES		<- sigmalist
    results$loglik <- loglik

    return(results)
  }
  if(class(lambdas) == "numeric")
  {
    Gibbs.method <- calculate_EM_Gibbs(combination, y, lambda1 = lambdas[1], lambda2 = lambdas[2], Theta=Theta, lower_upper = lower_upper, em.tol = em.tol, em.iter = max.elongation, gibbs.iter = 500, mc.iter = 1000, ncores = ncores)
  }
  else
  {
    Gibbs.method <- calculate_EM_Gibbs(combination, y, lambda1 = lambdas[combination,1], lambda2 = lambdas[combination,2], Theta=Theta, lower_upper = lower_upper, em.tol = em.tol, em.iter = max.elongation, gibbs.iter = 500, mc.iter = 1000, ncores = ncores)
  }
  invisible(return(Gibbs.method))
}


