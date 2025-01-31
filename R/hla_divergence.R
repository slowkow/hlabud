
#' Get a pairwise 20x20 distance matrix for all pairs of amino acids
#'
#' By default, we return the amino acid distance matrix by Grantham 1974 (doi:10.1126/science.185.4154.862).
#'
#' @param method `"grantham"` for the Grantham 1974 matrix or `"uniform"` for a matrix with ones on the non-diagonal.
#' @return A 20x20 symmetric matrix with positive numbers and zeros on the diagonal.
#' @examples
#' # By default, the Grantham 1974 matrix
#' amino_distance_matrix("grantham")
#'
#' # All ones, and zeros on the diagonal
#' amino_distance_matrix("uniform")
#' @export
amino_distance_matrix <- function(method = "grantham") {

  mat_grantham <- read.delim(text = "X	A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V
A	0	112	111	126	195	91	107	60	86	94	96	106	84	113	27	99	58	148	112	64
R	112	0	86	96	180	43	54	125	29	97	102	26	91	97	103	110	71	101	77	96
N	111	86	0	23	139	46	42	80	68	149	153	94	142	158	91	46	65	174	143	133
D	126	96	23	0	154	61	45	94	81	168	172	101	160	177	108	65	85	181	160	152
C	195	180	139	154	0	154	170	159	174	198	198	202	196	205	169	112	149	215	194	192
Q	91	43	46	61	154	0	29	87	24	109	113	53	101	116	76	68	42	130	99	96
E	107	54	42	45	170	29	0	98	40	134	138	56	126	140	93	80	65	152	122	121
G	60	125	80	94	159	87	98	0	98	135	138	127	127	153	42	56	59	184	147	109
H	86	29	68	81	174	24	40	98	0	94	99	32	87	100	77	89	47	115	83	84
I	94	97	149	168	198	109	134	135	94	0	5	102	10	21	95	142	89	61	33	29
L	96	102	153	172	198	113	138	138	99	5	0	107	15	22	98	145	92	61	36	32
K	106	26	94	101	202	53	56	127	32	102	107	0	95	102	103	121	78	110	85	97
M	84	91	142	160	196	101	126	127	87	10	15	95	0	28	87	135	81	67	36	21
F	113	97	158	177	205	116	140	153	100	21	22	102	28	0	114	155	103	40	22	50
P	27	103	91	108	169	76	93	42	77	95	98	103	87	114	0	74	38	147	110	68
S	99	110	46	65	112	68	80	56	89	142	145	121	135	155	74	0	58	177	144	124
T	58	71	65	85	149	42	65	59	47	89	92	78	81	103	38	58	0	128	92	69
W	148	101	174	181	215	130	152	184	115	61	61	110	67	40	147	177	128	0	37	88
Y	112	77	143	160	194	99	122	147	83	33	36	85	36	22	110	144	92	37	0	55
V	64	96	133	152	192	96	121	109	84	29	32	97	21	50	68	124	69	88	55	0")
  mat_grantham_rows <- mat_grantham[,1]
  mat_grantham <- mat_grantham[,2:ncol(mat_grantham)]
  rownames(mat_grantham) <- mat_grantham_rows
  mat_grantham <- as.matrix(mat_grantham)

  if (method == "grantham") {
    return(mat_grantham)
  }

  if (method == "uniform") {
    mat_uniform <- mat_grantham
    mat_uniform[mat_uniform > 0] <- 1
    return(mat_uniform)
  }

  stop(glue("method is '{method}', but it should be one of: 'grantham', 'uniform'"))

}

#' Calculate HLA divergence for each individual
#'
#' First, we convert each allele name (e.g. `A*01:01`) to an amino acid sequence.
#' The divergence is the sum of the distances between each pair of amino acids at each position, divided by the total sequence length.
#' The amino acid distance matrix we use is the one published by Grantham 1974 (doi:10.1126/science.185.4154.862), based on three physical properties of amino acids (composition, polarity, and molecular volume) that are correlated with an estimate of relative substitution frequency.
#' 
#' The code in this function is a translation of the [original Perl code](https://sourceforge.net/projects/granthamdist) by [Tobias Lenz](https://orcid.org/0000-0002-7203-0044), which was published in Pierini & Lenz 2018 MolBiolEvol (https://doi.org/10.1093/molbev/msy116).
#'
#' When comparing two amino acid sequences, only characters that are one of the 20 amino acids are considered in the divergence calculation, so gaps (and any other characters) do not count.
#'
#' @param alleles A character vector of comma-delimited alleles for each individual. We usually expect two alleles per individual, but it is possible to have more (or fewer) copies due to copy number alterations. This function still works when each individual has a different number of alleles.
#' @param method A pairwise amino acid matrix, or a method name: `"grantham"` or `"uniform"` to indicate which pairwise amino acid distance matrix to use. If you choose to pass a matrix, then it should be a 20x20 symmetric matrix with zeros on the diagonal, and the rownames and colnames should be the one-letter amino acid codes `A R N D C Q E G H I L K M F P S T W Y V`.
#' @param release Default is "latest". Should be a release name like "3.51.0".
#' @param verbose If TRUE, print messages along the way.
#' @seealso [hla_releases()] to get a complete list of all release names.
#' @seealso [amino_distance_matrix()] to get a amino acid distance matrix that you can use with `hla_divergence()`.
#' @return A dataframe with divergence for each individual.
#' @examples
#' my_genos <- c("A*23:01:12,A*24:550", "A*25:12N,A*11:27",
#'               "A*24:381,A*33:85", "A*01:01:,A*01:01,A*02:01")
#' hla_divergence(my_genos, method = "grantham")
#'
#' # This is equivalent
#' hla_divergence(my_genos, method = amino_distance_matrix("grantham"))
#' @export
hla_divergence <- function(
  alleles = c("A*01:01,A*02:01"), method = "grantham", release = "latest", verbose = FALSE,
  positions = NULL
) {

  aminos <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
    "F", "P", "S", "T", "W", "Y", "V")

  # alleles = c("HLA-A*01:01,HLA-A*03:01", "HLA-A*02:01,HLA-A*03:01")
  alleles <- str_replace_all(alleles, "HLA-", "")

  if (is.character(method) && method %in% c("uniform", "grantham")) {
    method_mat <- amino_distance_matrix(method)
  } else if (
    is.matrix(method) && all(dim(method) == c(20, 20)) &&
    all(aminos %in% colnames(method)) && all(aminos %in% rownames(method))
  ) {
    method_mat <- method
  } else {
    stop("method is invalid")
  }

  my_gene <- unique(str_split_fixed(alleles, "\\*", 2)[,1])
  stopifnot(length(my_gene) == 1)
  a <- hla_alignments(my_gene, release = release)

  retval <- sapply(alleles, function(my_ind) {

    my_alleles <- str_split(my_ind, ",")[[1]]
    my_alleles <- sapply(my_alleles, function(my_a) {
      my_regex <- str_replace_all(my_a, "\\*", "\\\\\\*")
      ix <- head(which(str_starts(rownames(a$alleles), my_regex)), 1)
      if (length(ix) == 0) {
        stop(glue("Input '{my_a}' not found"), call. = FALSE)
      }
      my_real <- rownames(a$alleles)[ix]
      if (verbose && my_real != my_a) {
        message(glue("Matching input '{my_a}' to '{my_real}'"))
      }
      my_real
    })
    my_mat <- a$alleles[my_alleles,]
    sequenceLength <- ncol(my_mat)
    allelePairs <- list()

    if (is.null(positions)) {
      positions <- 1:sequenceLength
    } else {
      positions <- as.character(positions)
      if (!all(positions %in% colnames(my_mat))) {
        stop(
          glue("positions must be a subset of column names from the hla_alignments()$alleles matrix:\n{paste(head(colnames(my_mat), 10), collapse=',')} ... {paste(tail(colnames(my_mat), 10), collapse=',')}"),
          call. = FALSE
        )
      }
    }

    for (i in 1:nrow(my_mat)) {
      for (j in 1:nrow(my_mat)) {
        distanceSum <- 0
        if (i != j) { # Calculate sum of pairwise aa distance between ith and jth allele
          # for (id in 1:sequenceLength) {
          for (id in positions) {
            x <- my_mat[i,id]
            y <- my_mat[j,id]
            if (x %in% aminos && y %in% aminos) {
              my_dist <- method_mat[x,y]
              if (length(my_dist) && !is.na(my_dist)) {
                distanceSum <- distanceSum + my_dist
              }
            }
          }
          distance <- distanceSum / sequenceLength # Average over sequence length
          x <- my_alleles[i]
          y <- my_alleles[j]
          if (!(x %in% allelePairs)) {
            allelePairs[[x]] <- list()
          }
          if (!(y %in% allelePairs)) {
            allelePairs[[y]] <- list()
          }
          allelePairs[[x]][[y]] <- distance
          allelePairs[[y]][[x]] <- distance
        }
      }
    }
    allelePairs

    my_combn <- combn(my_alleles, 2)
    my_dist <- 0
    for (i in 1:ncol(my_combn)) {
      x <- my_combn[1,i]
      y <- my_combn[2,i]
      d <- allelePairs[[x]][[y]]
      if (!is.null(d)) {
        my_dist <- my_dist + d
      }
    }
    my_dist <- my_dist / ncol(my_combn)

  })

  retval
}
