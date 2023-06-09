% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/igMobil.R
\name{Mobility}
\alias{Mobility}
\title{Calculates 14 Intergenerational Mobility Indices including single stage and transition matrix indices}
\usage{
Mobility(X, Y, grps = 4, showtmat = FALSE)
}
\arguments{
\item{X}{Parents SES data (Preferably income)}

\item{Y}{Childs SES data (Preferably income)}

\item{grps}{Number of percentiles to divide the data into for calculating the transiton matrix}

\item{showtmat}{Determines whether to print the transition matrix or not}
}
\description{
Takes the SES data for Parents and their children, and prints out the value of 14 different mobility indices. Transition Matrix is also constructed by dividing the input data into several percentiles
}
\references{
Savegnago, Marco. (2016). Igmobil: A Command for Intergenerational Mobility Analysis in Stata. Stata Journal. 16. 386-402.
}
\author{
Chitti Ankith Reddy,Bheemeshwar Reddy
}
