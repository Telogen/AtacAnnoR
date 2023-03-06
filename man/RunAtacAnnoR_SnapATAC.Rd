% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_Main.R
\name{RunAtacAnnoR_SnapATAC}
\alias{RunAtacAnnoR_SnapATAC}
\title{Run AtacAnnoR on a snap object in SnapATAC pipeline}
\usage{
RunAtacAnnoR_SnapATAC(
  ref_mtx,
  ref_celltype,
  ref_type = "sc",
  query_snapObj,
  threads = 10,
  verbose = TRUE
)
}
\arguments{
\item{ref_mtx}{a reference gene expression matrix whose rows are genes and columns are cells/samples}

\item{ref_celltype}{a vector of cell types whose orders are same as the columns of sc.mtx}

\item{ref_type}{which type of the reference matrix is
\itemize{
\item{'sc': }{the reference matrix is a single cell RNS-seq matrix}
\item{'bulk': }{the reference matrix is a bulk RNA-seq matrix}
}}

\item{query_snapObj}{the query snap object}

\item{threads}{the number of threads, default is 10}

\item{verbose}{whether to display messages and plots, default is TRUE}
}
\value{
Returns a new query snap object with cell types predicted by AtacAnnoR restored in
\code{query_snapObj@metaData$final_pred}.
}
\description{
Run AtacAnnoR on a snap object in SnapATAC pipeline. Note that before using this
function, the \code{gmat} and \code{pmat} should be stored in
\code{query_snapObj}.
}