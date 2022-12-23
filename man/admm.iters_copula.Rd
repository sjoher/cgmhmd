\name{admm.iters_copula}
\alias{admm.iters_copula}

\title{admm.iters_copula}

\description{This is a function that optimizes the \eqn{k} precision matrices.}

\usage{
admm.iters_copula(S, lambda1,lambda2,rho=1,rho.increment=1,weights,penalize.diagonal,maxiter = 1000,tol=1e-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{List containing \eqn{k} weighted version of \egn{Theta}.}
  \item{lambda1}{Value for lambda1.}
  \item{lambda2}{Value for lambda2.}
  \item{rho}{A step size parameter.  Large values decrease step size.}
  \item{rho.increment}{Factor to adjust rho.}
  \item{weights}{Determines the putative sample size of each class's data.  Allowed values: a vector with 
  length equal to the number of classes; "equal", giving each class weight 1; "sample.size",
  giving each class weight corresponding to its sample size.}
  \item{penalize.diagonal}{Choice whether to penalize the diagonal (TRUE) or not (FALSE)}
  \item{maxiter}{Maximum number of iterations}
  \item{tol}{Determines convergence criterion.}
}

\value{
\item{theta}{A list of the estimated inverse covariance matrices.}
\item{Z}{Optimized version of the \eqn{k} \egn{Theta}.}
\item{diff}{Loss value.}
\item{iter}{Number of iterations.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2022). Copula graphical models for heterogeneous mixed data. arXiv preprint, arXiv:2210.13140. \cr
}

\author{
  Sjoerd Hermes, Joost van Heerwaarden and Pariya Behrouzi\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
A <- list(matrix(runif(25), 5, 5),matrix(runif(25), 5, 5),matrix(runif(25), 5, 5))
S <- list()
for(k in 1:3){
  S[[k]] <- cor(A[[k]])
}
lam1 <- 0.1
lam2 <- 0.5
lam1.bl = penalty.as.matrix_copula(lam1,dim(S[[1]])[2],penalize.diagonal=penalize.diagonal)
lam2.bl = penalty.as.matrix_copula(lam2,dim(S[[1]])[2],penalize.diagonal=TRUE)
Thetabl = admm.iters_copula(S, lam1.bl,lam2.bl,rho=1,weights="equal",penalize.diagonal=FALSE,maxiter=1000,tol=1e-5)
}
