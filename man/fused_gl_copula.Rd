\name{fused_gl_copula}
\alias{fused_gl_copula}

\title{fused_gl_copula}

\description{This is a function that computes estimates for the \eqn{K} precision matrices using the fused group lasso penalty.}

\usage{
fused_gl_copula(S, Y, lambda1,lambda2,rho=1,weights="equal",penalize.diagonal=FALSE,maxiter=1000,tol=1e-5,warm=NULL,return.whole.theta=TRUE,truncate=1e-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{List of \eqn{K} estimated covariance matrices.}
  \item{Y}{List of \eqn{K} data matrices.}
  \item{lambda1}{Value for lambda1.}
  \item{lambda2}{Value for lambda2.}
  \item{rho}{A step size parameter.  Large values decrease step size.}
  \item{weights}{Determines the putative sample size of each class's data.  Allowed values: a vector with 
  length equal to the number of classes; "equal", giving each class weight 1; "sample.size",
  giving each class weight corresponding to its sample size.}
  \item{penalize.diagonal}{Choice whether to penalize the diagonal (TRUE) or not (FALSE)}
  \item{maxiter}{Maximum number of iterations}
  \item{tol}{Determines convergence criterion.}
  \item{warm}{Input a warm start to theta in the form of a K-length list of pXp matrices.}
  \item{return.whole.theta}{If TRUE, each class's inverse covariance matrix is returned whole.  If FALSE, the inverse covariance
  matrix is only returned over the connected nodes, and only the diagonal of the matrix is returned over
  the unconnected nodes.}
}

\value{
\item{theta}{A list of the estimated inverse covariance matrices.}
\item{connected}{Identities of the connected nodes.}
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
  S[[k]] <- cov(A[[k]])
}
lam1 <- 0.1
lam2 <- 0.5
R.fgl <- fused_gl_copula(S = S, Y = A, lam1, lam2)
}