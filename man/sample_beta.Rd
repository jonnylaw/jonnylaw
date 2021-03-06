% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_dlm.R
\name{sample_beta}
\alias{sample_beta}
\title{Sample beta}
\usage{
sample_beta(ys, xs, theta, mod, beta_mu, beta_sigma)
}
\arguments{
\item{ys}{p x n matrix of observation}

\item{xs}{p x m matrix of time-independent covariates}

\item{theta}{latent state}

\item{mod}{dlm model}

\item{beta_mu}{prior mean of the coefficients}

\item{beta_sigma}{prior standard deviation of the coefficients}
}
\value{

}
\description{
Sample beta
}
