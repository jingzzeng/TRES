\name{PMSE}
\alias{PMSE}

\title{
Predictions and Mean squared error for tensor predictor regression
}
\description{
Evaluate tensor predictor regression through prediction mean squared error.
}
\usage{
PMSE(Xn, Yn, Bhat)
}

\arguments{
  \item{Xn}{A predictor tensor.
}
\item{Yn}{A response vector.}
\item{Bhat}{An estimation of coefficient tensor.}
}

\value{
\item{mse}{Mean squared error. Defined as \eqn{trace(\sum_{i=1}^n(\mathbf{Y}_i-\hat{\mathbf{Y}}_i)(\mathbf{Y}_i-\hat{\mathbf{Y}}_i)'/n)}, where \eqn{\hat{\mathbf{Y}}_i} is the predictions. }
\item{Yhat}{The predictions of tensor predictor regression.}
}
