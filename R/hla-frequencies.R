#' Get HLA frequences from Allele Frequency Net Database (AFND)
#'
#' Download and read a table of HLA allele frequencies from the [Allele Frequency Net Database (AFND)](http://www.allelefrequencies.net/).
#'
#' If you use this data, please cite the latest manuscript about Allele Frequency Net Database:
#'
#' - Gonzalez-Galarza FF, McCabe A, Santos EJMD, Jones J, Takeshita L, Ortega-Rivera ND, et al. [Allele frequency net database (AFND) 2020 update: gold-standard data classification, open access genotype data and new query tools.](https://pubmed.ncbi.nlm.nih.gov/31722398) Nucleic Acids Res. 2020;48: D783–D788. doi:10.1093/nar/gkz1029
#'
#' @param verbose If TRUE, print messages along the way.
#' @return A dataframe with HLA allele frequencies for all genes.
#' @examples
#' \donttest{
#' hla_frequencies()
#' }
#' @export
hla_frequencies <- function(verbose = FALSE) {
  hlabud_dir <- get_hlabud_dir()
  my_file <- file.path(hlabud_dir, "afnd.tsv")
  my_url <- "https://github.com/slowkow/allelefrequencies/raw/main/afnd.tsv"
  if (!file.exists(my_file)) {
    if (verbose) { message(glue("Downloading {my_url}")) }
    lines <- readLines(my_url)
    if (verbose) { message(glue("Writing {my_file}")) }
    writeLines(lines, my_file)
  }
  if (verbose) { message(glue("Reading {my_file}")) }
  d <- as_tibble(read.delim(my_file))
  d$indivs_over_n <- as.numeric(str_remove(d$indivs_over_n, "\\(\\*\\)"))
  d$alleles_over_2n <- as.numeric(str_remove(d$alleles_over_2n, "\\(\\*\\)"))
  d$n <- as.numeric(str_remove_all(d$n, ","))
  return(d)
}

