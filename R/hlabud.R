
#' @keywords internal
mkdir <- function(path) {
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

#' @keywords internal
get_release <- function(release = "latest") {
  releases <- hla_releases()
  if (release == "latest") {
    # If the user didn't set it, then set it to the latest release
    release <- releases[1]
  }
  # don't allow unrecognized releases
  if (!release %in% releases) {
    stop("Unrecognized release '{release}' not in releases: {paste(releases, ' ')}")
  }
  return(release)
}

#' Get the name of the folder for caching downloaded IMGTHLA files
#'
#' Get the folder name from `getOption("hlabud_dir")` or else use the
#' [rappdirs](https://github.com/r-lib/rappdirs) package to choose an
#' appropriate folder for your operating system. The folder will be created
#' automatically if it does not exist. And the `hlabud_dir` option will bet set
#' to that folder.
#'
#' Linux:
#'
#'     ~/.local/share/hlabud
#'
#' Mac:
#'
#'     ~/Library/Application Support/hlabud
#'
#' Windows:
#'
#'     C:\Documents and Settings\<User>\Application Data\slowkow\hlabud
#'
#' @examples
#' \dontrun{
#' hlabud_dir <- get_hlabud_dir()
#' }
#' @return The name of the folder.
#' @export
get_hlabud_dir <- function() {
  hlabud_dir <- getOption("hlabud_dir")
  if (!is.null(hlabud_dir) && nchar(hlabud_dir) > 1) {
    mkdir(hlabud_dir)
  } else {
    appname <- "hlabud"
    appauthor <- "slowkow"
    hlabud_dir <- path.expand(user_data_dir(appname, appauthor))
    mkdir(hlabud_dir)
    options(hlabud_dir = hlabud_dir)
  }
  return(hlabud_dir)
}

#' Download and unpack a tarball release of IMGTHLA from GitHub
#' 
#' The release tarball from Github is unpacked into the
#' `getOption("hlabud_dir")` folder.
#'
#' Note that the latest releases are more than 100 MB in size, so the download
#' might take a while on slow connections.
#'
#' @param release Default is "latest". Should be a release name like "3.51.0".
#' @param overwrite If TRUE, overwrite existing files in the release folder.
#' @param verbose If TRUE, print messages along the way.
#' @examples
#' \dontrun{
#' install_hla()
#' install_hla("3.51.0")
#' install_hla("3.51.0", verbose = TRUE)
#' # Change the install directory
#' options(hlabud_dir = "path/to/my/dir")
#' install_hla()
#' }
#' @return NULL
#' @export
install_hla <- function(release = "latest", overwrite = FALSE, verbose = FALSE) {

  hlabud_dir <- get_hlabud_dir()

  if (verbose) { message(glue("Fetching releases from GitHub")) }
  tags <- read_json("https://api.github.com/repos/ANHIG/IMGTHLA/tags")
  writeLines(toJSON(tags, pretty = TRUE, auto_unbox = TRUE), file.path(hlabud_dir, "tags.json"))

  release <- get_release(release)
  ix <- which(str_detect(sapply(tags, "[[", "name"), glue("v{release}-")))

  url <- tags[[ix]]$tarball_url
  tar_file <- file.path(hlabud_dir, glue("{basename(url)}.tar.gz"))
  release_dir <- file.path(hlabud_dir, release)

  if (!file.exists(tar_file) && (overwrite || !file.exists(release_dir))) {
    options(timeout = max(300, getOption("timeout")))
    message(glue("Downloading {url}"))
    # download.file(url, destfile = tar_file)
    curl_download(url, tar_file, quiet = !verbose)
  }

  if (file.exists(tar_file) && (overwrite || !file.exists(release_dir))) {
    dir.create(release_dir)
    message(glue("Unpacking {tar_file}"))
    untar(tar_file, exdir = release_dir, extras = '--strip-components=1')
  } else {
    message(glue("Using existing installation: {release_dir}"))
  }

}

#' Get IMGTHLA gene names
#'
#' Retrieve the contents of `ANHIG/IMGTHLA/alignments` and return a list of gene
#' names derived from the txt files in that folder.
#'
#' See the files at: https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
#'
#' @param release Default is "latest". Should be a release name like "3.51.0".
#' @param overwrite Overwrite the existing `genes.json` file with a new one from GitHub
#' @param verbose If TRUE, print messages along the way.
#' @return A character vector of HLA gene names like "DRB1"
#' @examples
#' \donttest{
#' hla_genes() 
#' }
#' @export
hla_genes <- function(release = "latest", overwrite = FALSE, verbose = FALSE) {
  hlabud_dir <- get_hlabud_dir()
  tags_file <- file.path(hlabud_dir, "tags.json")
  release <- get_release(release)
  genes_file <- file.path(hlabud_dir, release, "genes.json")
  # get the shasum
  tags <- read_json(tags_file)
  ix <- which(str_detect(sapply(tags, "[[", "name"), glue("v{release}-")))
  sha <- tags[[ix]]$commit$sha
  # get the genes.json
  if (overwrite || !file.exists(genes_file)) {
    j <- read_json(glue("https://api.github.com/repos/ANHIG/IMGTHLA/contents/alignments?ref={sha}"))
    if (verbose) { message(glue("Writing {genes_file}")) }
    mkdir(dirname(genes_file))
    writeLines(toJSON(j, pretty = TRUE, auto_unbox = TRUE), genes_file)
  }
  j <- read_json(genes_file)
  genes <- sapply(j, "[[", "name")
  genes <- unique(str_extract(genes, "^[^_]+"))
  return(genes)
}

#' Get a table of allele names for a particular IMGTHLA release
#'
#' Download a list of all allele names for all HLA genes for a particular
#' IMGTHLA release.
#'
#' @param release Default is "latest". Should be a release name like "3.51.0".
#' @param overwrite Overwrite the existing `alleles.json` file and `Allelelist.{version}.txt` file
#' @param verbose If TRUE, print messages along the way.
#' @return A data frame with HLA allele ids and names
#' @examples
#' \donttest{
#' head(hla_alleles())
#' }
#' @export
hla_alleles <- function(release = "latest", overwrite = FALSE, verbose = FALSE) {
  hlabud_dir <- get_hlabud_dir()
  alleles_file <- file.path(hlabud_dir, "alleles.json")
  release <- get_release(release)
  if (overwrite || !file.exists(alleles_file)) {
    j <- read_json("https://api.github.com/repos/ANHIG/IMGTHLA/contents/allelelist")
    if (verbose) { message(glue("Writing {alleles_file}")) }
    writeLines(toJSON(j, pretty = TRUE, auto_unbox = TRUE), alleles_file)
  }
  j <- read_json(alleles_file)
  #
  # Convert a release like "3.51.0" to a slug like "3510"
  my_slug <- str_replace_all(release, "\\.", "")
  my_name <- glue("Allelelist.{my_slug}.txt")
  #
  my_names <- sapply(j, function(x) x$name)
  if (!(my_name %in% my_names)) {
    warning(glue("unrecognized release name '{my_name}'"))
    return(NULL)
  }
  #
  out_dir <- file.path(hlabud_dir, "allelelist")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, my_name)
  if (overwrite || !file.exists(out_file)) {
    curl_download(
      glue("https://github.com/ANHIG/IMGTHLA/raw/Latest/allelelist/{my_name}"),
      out_file
    )
  }
  alleles <- read.table(out_file, comment.char = "#", sep = ",", header = TRUE)
  return(alleles)
}

#' Get IMGTHLA releases from GitHub
#'
#' Retrieve the tags from the ANHIG/IMGTHLA repo and get the associated release
#' names.
#'
#' @param overwrite Overwrite the existing tags.json file in `file.path(getOption("hlabud_dir"), "tags.json")`with a new one from the [IMGTHLA](https://github.com/ANHIG/IMGTHLA) GitHub repo
#' @return A character vector of release names like "3.51.0"
#' @examples
#' \donttest{
#' hla_releases() 
#' }
#' @export
hla_releases <- function(overwrite = FALSE) {
  hlabud_dir <- get_hlabud_dir()
  tags_file <- file.path(hlabud_dir, "tags.json")
  if (overwrite || !file.exists(tags_file)) {
    j <- read_json("https://api.github.com/repos/ANHIG/IMGTHLA/tags")
    writeLines(toJSON(j, pretty = TRUE, auto_unbox = TRUE), tags_file)
  }
  j <- read_json(tags_file)
  release_names <- sapply(j, function(x) x$name)
  releases <- str_extract(release_names, "[\\d.]+")
  return(releases)
}

#' Get aligned sequences in a dataframe
#'
#' Here are the conventions used for alignments ([EBI IMGT-HLA help page](https://www.ebi.ac.uk/ipd/imgt/hla/alignment/help/)):
#' * The entry for each allele is displayed in respect to the reference sequences.
#' * Where identity to the reference sequence is present the base will be displayed as a hyphen (-).
#' * Non-identity to the reference sequence is shown by displaying the appropriate base at that position.
#' * Where an insertion or deletion has occurred this will be represented by a period (.).
#' * If the sequence is unknown at any point in the alignment, this will be represented by an asterisk (*).
#' * In protein alignments for null alleles, the 'Stop' codons will be represented by a hash (X).
#' * In protein alignments, sequence following the termination codon, will not be marked and will appear blank.
#' * These conventions are used for both nucleotide and protein alignments.
#'
#' @param gene The name of a gene like "DRB1"
#' @param type The type of sequence, one of "prot", "nuc", "gen"
#' @param release Default is "latest". Should be a release name like "3.51.0".
#' @param verbose If TRUE, print messages along the way.
#' @return A dataframe.
#' @examples
#' \donttest{
#' a <- hla_alignments("DRB1")
#' head(a$sequences)
#' a$alleles[1:6,1:6]
#' a$onehot[1:6,1:6]
#' }
#' @export
hla_alignments <- function(gene = "DRB1", type = "prot", release = "latest", verbose = FALSE) {
  hlabud_dir <- get_hlabud_dir()
  tags_file <- file.path(hlabud_dir, "tags.json")
  release <- get_release(release)
  if (!type %in% c("nuc", "gen", "prot")) {
    stop("Unrecognized type '{type}' not in: nuc gen prot")
  }
  if (verbose) { message(glue("hlabud is using IMGTHLA release {release}")) }
  # 
  genes_file <- file.path(hlabud_dir, release, "genes.json")
  if (!file.exists(genes_file)) {
    hla_genes()
  }
  genes <- read_json(genes_file)
  genes <- genes[which(str_detect(sapply(genes, "[[", "name"), gene))]
  #
  my_names <- sapply(genes, "[[", "name")
  my_urls <- sapply(genes, "[[", "download_url")
  i <- which(str_detect(my_names, type))
  if (length(i) == 0) {
    message(glue("IMGT/HLA does not provide {type} for {gene}"))
    return()
  }
  #
  my_name <- my_names[i]
  my_type <- str_replace(my_name, "^.+_(\\w+)\\..+", "\\1")
  my_file <- file.path(hlabud_dir, release, "alignments", my_name)
  my_url <- my_urls[i]
  if (!file.exists(my_file)) {
    if (verbose) { message(glue("Downloading {my_url}")) }
    lines <- readLines(my_url)
    if (verbose) { message(glue("Writing {my_file}")) }
    mkdir(dirname(my_file))
    writeLines(lines, my_file)
  }
  if (verbose) { message(glue("Reading {my_file}")) }
  return(read_alignment(my_file))
}

#' Read an alignment file `*_prot.txt` file from IMGTHLA.
#'
#' The prot file has the amino acid sequence for each HLA allele.
#'
#' @return A list with a dataframe called `sequences` and two matrices `alleles` and `onehot` The data frame has two columns:
#' * allele: the name of the allele, e.g., `DQB*01:01`
#' * seq: the amino acid sequence
#' The matrix `alleles` has one row for each allele, and one column for each position, with the values representing the residues at each position in each allele.
#' The matrix `onehot` has a one-hot encoding of the variants that distinguish the alleles, with one row for each allele and one column for each amino acid at each position.
#' @param my_file File name for a txt file from IMGTHLA like "DQB1_prot.txt"
#' @examples
#' my_file <- file.path(
#'   "https://github.com/ANHIG/IMGTHLA/raw",
#'   "5f2c562056f8ffa89aeea0631f2a52300ee0de17",
#'   "alignments/DRB1_prot.txt"
#' )
#' a <- read_alignment(my_file)
#' head(a$sequences)
#' a$alleles[1:5,1:5]
#' a$onehot[1:5,1:5]
#' @export
read_alignment <- function(my_file) {
  my_type <- str_remove(str_split_fixed(basename(my_file), "_", 2)[,2], ".txt")
  my_gene <- str_split_fixed(basename(my_file), "_", 2)[,1]
  lines <- readLines(my_file)
  #
  if (my_type == "prot") {
    # Many amino acids are located before the position labeled as "1"
    al_i <- which(str_detect(lines, "Prot.+ 1$"))
    pre <- lines[al_i]
    pre_i <- str_locate(pre, "-")[1,1]
    pre_j <- str_locate(pre, "1")[1,1]
    n_pre <- nchar(str_replace_all(substr(lines[al_i + 2], pre_i, pre_j - 1), " ", ""))
  } else if (my_type == "nuc") {
    al_i <- which(str_detect(lines, "AA codon"))
    n_pre <- 3 * abs(as.numeric(str_remove(lines[al_i][1], " *AA codon +")))
  } else if (my_type == "gen") {
    al_i <- which(str_detect(lines, "gDNA"))
    n_pre <- abs(as.numeric(str_remove(lines[al_i][1], " *gDNA +")))
  }
  # Convert lines to a simple data frame
  my_regex <- glue("^ {my_gene}\\\\*")
  al <- lines[str_detect(lines, my_regex)]
  al <- str_replace_all(al, "\\|", "")
  al <- str_split_fixed(al, " +", 3)
  al <- al[,2:3]
  colnames(al) <- c("allele", "seq")
  al <- as.data.frame(al)
  al <- al %>% group_by(.data$allele) %>% mutate(id = sprintf("V%s", seq(n()))) %>% ungroup()
  al <- al %>% pivot_wider(names_from = "id", values_from = "seq", values_fill = "")
  al <- al %>% unite("seq", starts_with("V"), sep = "")
  al$seq <- str_replace_all(al$seq, " ", "")
  #
  oh <- get_onehot(al, n_pre)
  return(list(
    sequences = as_tibble(al),
    alleles = oh$alleles,
    onehot = oh$onehot
  ))
}

#' Make a one-hot encoded matrix from a dataframe of amino acid
#' sequences.
#'
#' @param al A dataframe with columns allele, seq
#' @param n_pre The number of amino acid sequences before position 1.
#' @keywords internal
get_onehot <- function(al, n_pre) {
  seq_chars <- str_split(al$seq, "")
  ref_chars <- seq_chars[[1]]
  max_chars <- max(sapply(seq_chars, length))
  min_chars <- min(sapply(seq_chars, length))
  for (i in 2:length(seq_chars)) {
    ix <- which(seq_chars[[i]] == "-")
    ix <- ix[ix < length(ref_chars)]
    seq_chars[[i]][ix] <- ref_chars[ix]
    seq_len <- length(seq_chars[[i]])
    if (seq_len < max_chars) {
      seq_chars[[i]] <- c(seq_chars[[i]], rep("*", max_chars - seq_len))
    }
  }
  if (length(ref_chars) < max_chars) {
    ref_chars <- c(ref_chars, rep("*", max_chars - length(ref_chars)))
    seq_chars[[1]] <- ref_chars
  }
  # Create a one-hot-encoded matrix
  # with allele names (rows) and positions (columns)
  # 
  #         P1_* P1_M P2_* P2_A P3_*
  # A*01:01    0    1    0    1    0
  # A*01:02    0    1    0    1    0
  # A*02:01    0    1    0    1    0
  # A*02:02    0    1    0    1    0
  # A*02:03    0    1    0    1    0
  # 
  #######################################################################
  alleles <- do.call(rbind, seq_chars)
  rownames(alleles) <- al$allele
  colnames(alleles) <- str_replace_all(
    sprintf("pos%s", c(-n_pre:-1, 1:(ncol(alleles) - n_pre))), "-", "n"
  )
  # colnames(alleles) <- sprintf("P%s", seq(ncol(alleles)))
  # keep positions with more than 1 allele
  alleles <- alleles[,apply(alleles, 2, function(x) length(unique(x))) > 1, drop = FALSE]
  alleles <- as.data.frame(alleles)
  for (i in seq(ncol(alleles))) {
    alleles[,i] <- as.factor(alleles[,i])
  }
  p <- onehot(alleles, max_levels = 20)
  retval <- predict(p, alleles)
  rownames(retval) <- al$allele
  # Discard positions where we don't know the allele
  retval <- retval[,!str_detect(colnames(retval), "\\*")]
  colnames(retval) <- str_replace(colnames(retval), "=", "_")
  return(list(alleles = as.matrix(alleles), onehot = retval))
}

#' Dosage
#'
#' For each genotype, return the the dosage for each amino acid (or nucleotide)
#' at each position.
#'
#' Each genotype should be represented like `"HLA-A*01:01,HLA-A*01:01"`
#'
#' By default, the returned data frame is filtered to exclude:
#' * positions where all input genotypes have the same allele
#' * positions that are identical to previous positions
#'
#' @param genotypes Input character vector with one genotype for each individual.
#' @param alleles A one-hot encoded matrix with one row per allele and one
#' column per amino acid position (or nucleotide position).
#' @param drop_constants Filter out constant amino acid positions by default.
#' @param drop_duplicates Filter out duplicate amino acid positions by default.
#' @returns A data frame with one row for each input genotype.
#' @examples
#' DRB1_file <- file.path(
#'   "https://github.com/ANHIG/IMGTHLA/raw",
#'   "5f2c562056f8ffa89aeea0631f2a52300ee0de17",
#'   "alignments/DRB1_prot.txt"
#' )
#' a <- read_alignment(DRB1_file)
#' genotypes <- c(
#'   "DRB1*12:02:02:03,DRB1*12:02:02:03,DRB1*14:54:02",
#'   "DRB1*04:174,DRB1*15:152",
#'   "DRB1*04:56:02,DRB1*15:01:48",
#'   "DRB1*14:172,DRB1*04:160",
#'   "DRB1*04:359,DRB1*04:284:02"
#' )
#' dosage <- dosage(genotypes, a$onehot)
#' dosage[,1:5]
#' @export
dosage <- function(genotypes, alleles, drop_constants = TRUE, drop_duplicates = TRUE) {
  dosages <- matrix(0, ncol = ncol(alleles), nrow = length(genotypes))
  for (i in seq_along(genotypes)) {
    # Split a string of genotypes like "HLA-A*01:01,HLA-A*01:01"
    a <- str_split(genotypes[i], ",")[[1]]
    for (my_a in a) {
      # Find the first row in alleles where the prefix matches our genotype
      ix <- which(str_starts(rownames(alleles), fixed(my_a)))
      if (length(ix) > 0) {
        dosages[i,] <- dosages[i,] + as.numeric(alleles[ix[1],])
      } else {
        warning(glue("'{my_a}' not found in rownames"))
      }
    }
  }
  rownames(dosages) <- genotypes
  colnames(dosages) <- colnames(alleles)
  if (drop_constants) {
    # Select positions where we observe > 1 possible dosage [0, 1, 2]
    ix <- apply(dosages, 2, function(x) length(unique(x)))
    dosages <- dosages[, ix > 1, drop = FALSE]
    if (!any(ix)) {
      return(dosages)
    }
  }
  if (drop_duplicates) {
    # Discard positions that are identical to previous positions
    ix <- !duplicated(apply(dosages, 2, function(x) paste(x, collapse = ",")))
    dosages <- dosages[,which(ix)]
  }
  return(dosages)
}

