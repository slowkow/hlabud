% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hla-frequencies.R
\name{hla_frequencies}
\alias{hla_frequencies}
\title{Get HLA frequences from Allele Frequency Net Database (AFND)}
\usage{
hla_frequencies(verbose = FALSE)
}
\arguments{
\item{verbose}{If TRUE, print messages along the way.}
}
\value{
A dataframe with HLA allele frequencies for all genes.
}
\description{
Download and read a table of HLA allele frequencies from the \href{http://www.allelefrequencies.net/}{Allele Frequency Net Database (AFND)}.
}
\details{
If you use this data, please cite the latest manuscript about Allele Frequency Net Database:
\itemize{
\item Gonzalez-Galarza FF, McCabe A, Santos EJMD, Jones J, Takeshita L, Ortega-Rivera ND, et al. \href{https://pubmed.ncbi.nlm.nih.gov/31722398}{Allele frequency net database (AFND) 2020 update: gold-standard data classification, open access genotype data and new query tools.} Nucleic Acids Res. 2020;48: D783–D788. doi:10.1093/nar/gkz1029
}
}
\examples{
\donttest{
hla_frequencies()
}
}
