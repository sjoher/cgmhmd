\name{data_sim}
\alias{data_sim}

\title{data_sim}

\description{
Simulate mixed multi-group data.
}

\usage{
data_sim(network, n, p, K, ncat, rho, gamma_g = NULL, gamma_o, gamma_b = NULL,
gamma_p = NULL, prob = NULL, nclass = NULL)
}

\arguments{
  \item{network}{Type of network, either "circle", "Random", "Cluster",
  "Scale-free", "AR1" or "AR2".}
  \item{n}{Number of observations.}
  \item{p}{Number of variables.}
  \item{K}{Number of groups.}
  \item{ncat}{Number of categories for ordinal variables.}
  \item{rho}{Dissimilarity parameter inducing dissimilarity between the K datasets.}
  \item{gamma_g}{Proportion of Gaussian variables in the data.}
  \item{gamma_o}{Proportion of ordinal variables in the data.}
  \item{gamma_b}{Proportion of binomial variables in the data.}
  \item{gamma_p}{Proportion of Poisson variables in the data..}
  \item{prob}{Edge occurency probability in random graph.}
  \item{nclass}{Number of clusters in cluster graph.}

}

\value{
\item{z}{A list of \eqn{K} \eqn{n} by \eqn{p} matrices representing the latent
Gaussian transformed (observed) data.}
\item{theta}{A list of \eqn{K} \eqn{n} by \eqn{p} matrices representing the
precision matrices corresponding to the latent Gaussian (unobserved) data.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2022). Copula graphical models for heterogeneous mixed data. arXiv preprint, arXiv:2210.13140. \cr
}

\author{
  Sjoerd Hermes, Joost van Heerwaarden and Pariya Behrouzi\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
data_sim(network = "Random", n = 10, p = 50, K = 3, ncat = 6, rho = 0.25, gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05)
}
