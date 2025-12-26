\name{NNRPanel-package}
\alias{NNRPanel-package}
\alias{NNRPanel}
\docType{package}

\title{
NNRPanel: Nonlinear Panel Estimation with Interactive Fixed Effects
}

\description{
The \pkg{NNRPanel} package implements nonlinear panel estimators with 
interactive fixed effects. The proposed method involves two steps. 
In the first step, we convexify the optimization problem using nuclear norm regularization (NNR) 
and obtain preliminary NNR estimators of the parameters, including the fixed effects. 
Then, we find the global solution of the original optimization problem using a 
standard gradient descent method initialized at these preliminary estimates.
The proposed estimator avoids solving a high-dimensional optimization problem 
and can be feasibly computed in large nonlinear panels. 
We also provide analytical and sample-splitting bias corrections. 
}

\details{
\pkg{NNRPanel} is designed for panel data, and support for network data will be 
added in the near future.

The package supports logistic, Probit, and Poisson outcome models. It also provides
a data-driven tuning-parameter selection method and a simple procedure for 
determining the number of factors.

Key functionality includes:
\itemize{
  \item Data-driven tuning-parameter selection used in the first step estimation. 
  \item Nuclear-norm–regularized first-step estimation of nonlinear panel models.
  \item Data-driven selection of the number of latent factors.
  \item Second-step local MLE estimation using first-step estimators as initial value. 
  \item Analytical and sample-splitting bias corrections for high-dimensional panel models.
}

\strong{Main exported functions:}
\itemize{
  \item \code{NNRPanel_estimate()}: Main estimator pipeline (NNR \eqn{\rightarrow} local MLE \eqn{\rightarrow} bias correction).
  \item \code{calculate_tuning_parameter()}: Data-driven tuning parameter \eqn{\varphi}.
  \item \code{determine_num_factor()}: Estimate the factor dimension.
}

This package is intended for econometric applications involving nonlinear
panel models with interactive fixed effects  when the number of
individuals and/or time periods is large.
}

\author{
Weisheng Zhang

Maintainer: Weisheng Zhang <wszhang.econ@gmail.com>
}



\keyword{package}



\examples{
\dontrun{
# Simulated example (Logit)
data_list <- gen_data_logit(N = 200, T = 100, num_factor = 2, num_Z = 1, beta = c(0.2))
data_frame <- list_to_data_frame(data_list)

# Main estimator with data-driven factor selection
res <- NNRPanel_estimate(data_frame, func = "logit", R_max = 5, R_true = 2, delta = 0.05)
}
}

\references{
Zeleneev, Andrei, and Weisheng Zhang. (2025).  
\emph{Tractable Estimation of Nonlinear Panels with Interactive Fixed Effects}.  
}