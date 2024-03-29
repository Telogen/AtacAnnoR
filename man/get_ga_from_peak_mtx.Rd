% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Round1Anno.R
\name{get_ga_from_peak_mtx}
\alias{get_ga_from_peak_mtx}
\title{Get gene activity matrix from the peak counts matrix}
\usage{
get_ga_from_peak_mtx(peak_counts, gene_gr, upstream = 2000, threads = 10)
}
\arguments{
\item{peak_counts}{The peak counts matrix whose rows are peaks and columns are cells}

\item{gene_gr}{The GRanges objects of all coding genes genome, can be got from \code{get_gene_gr()}}

\item{upstream}{Number of bases to extend upstream of the TSS, default is 2000}

\item{threads}{Number of threads}
}
\value{
Return a gene activity matrix whose rows are genes and columns are cells.
}
\description{
Get gene activity matrix from the peak counts matrix
}
