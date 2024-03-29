% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_utilities.R
\name{sparse2DT}
\alias{sparse2DT}
\title{Convert a sparse matrix to a data.table}
\usage{
sparse2DT(matrix, with.names = FALSE, single.cell = FALSE, key = TRUE)
}
\arguments{
\item{matrix}{sparse matrix to be used}

\item{with.names}{return the row and column names from the sparse matrix \link{default:FALSE}}

\item{single.cell}{returns the row names as Genes and column names as cells if from a single cell dataset \link{default: FALSE}}

\item{key}{logical for whether the data.table key should be set \link{default:FALSE}}
}
\description{
Convert a sparse matrix to a data.table
}
