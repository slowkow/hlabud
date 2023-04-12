
#' @keywords internal
mkdir <- function(path) {
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

#' @keywords internal
setup_hlabud_dir <- function() {
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
#' Fetch a list of releases from Github. By default, install the latest
#' release to the folder appropriate for your operating system:
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
#' A typical release is about 120 MB in size, so the download can take a few
#' minutes.
#' 
#' We can get or set the folder with the `hlabud_dir` option:
#'
#'     my_dir <- getOption("hlabud_dir")
#'     options(hlabud_dir = my_dir)
#'
#' The release tarball from Github is unpacked into the folder.
#'
#' Other functions in the hlabud package will use the unpacked data, or else
#' they will automatically download the minimum necessary files.
#'
#' @param release Either "latest" or else a release name like "3.51.0"
#' @param quiet If FALSE, print messages along the way.
#' @examples
#' \dontrun{
#' install_hla()
#' install_hla("3.51.0")
#' install_hla("3.51.0", quiet = TRUE)
#' # Change the install directory
#' options(hlabud_dir = "path/to/my/dir")
#' install_hla()
#' }
#' @return NULL
#' @export
install_hla <- function(release = "latest", quiet = FALSE) {

  hlabud_dir <- setup_hlabud_dir()

  message(glue("Fetching releases from GitHub"))
  j <- read_json("https://api.github.com/repos/ANHIG/IMGTHLA/tags")
  writeLines(toJSON(j, pretty = TRUE, auto_unbox = TRUE), file.path(hlabud_dir, "tags.json"))
  release_names <- sapply(j, function(x) x$name)
  releases <- str_extract(release_names, "[\\d.]+")
  message(glue("{length(releases)} releases, latest release: {releases[1]}"))

  if (release == "latest") {
    my_j <- j[[1]]
    release <- releases[1]
  } else {
    ix <- which(str_detect(releases, release))
    if (length(ix) != 1) {
      stop(glue("Could not find '{release}' in: {paste(releases, collapse=' ')}"))
    }
    my_j <- j[[ix]]
  }

  hlabud_release <- getOption("hlabud_release")
  if (is.null(hlabud_release) || !hlabud_release %in% releases) {
    options(hlabud_release = release) 
  }

  url <- my_j$tarball_url
  tar_file <- file.path(hlabud_dir, glue("{basename(url)}.tar.gz"))

  if (!file.exists(tar_file)) {
    options(timeout = max(300, getOption("timeout")))
    message(glue("Downloading {url}"))
    # download.file(url, destfile = tar_file)
    curl_download(url, tar_file, quiet = quiet)
  }

  release_dir <- file.path(hlabud_dir, release)
  if (!file.exists(release_dir)) {
    dir.create(release_dir)
    message(glue("Unpacking {tar_file}"))
    untar(tar_file, exdir = release_dir, extras = '--strip-components=1')
  } else {
    message(glue("Using existing installation: {release_dir}"))
  }

  # prot_files <- Sys.glob(glue("{dirname(tar_file)}/*/alignments/*_prot.txt"))
}

#' Get IMGTHLA releases from GitHub
#'
#' @return A character vector of release names like "3.51.0"
#' @param overwrite Overwrite the existing tags.json file with a new one from GitHub
#' @examples
#' \donttest{
#' hla_releases() 
#' }
#' @export
hla_releases <- function(overwrite = FALSE) {
  hlabud_dir <- setup_hlabud_dir()
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
#' Here are the conventions used for alignments:
#' * The entry for each allele is displayed in respect to the reference sequences.
#' * Where identity to the reference sequence is present the base will be displayed as a hyphen (-).
#' * Non-identity to the reference sequence is shown by displaying the appropriate base at that position.
#' * Where an insertion or deletion has occurred this will be represented by a period (.).
#' * If the sequence is unknown at any point in the alignment, this will be represented by an asterisk (*).
#' * In protein alignments for null alleles, the 'Stop' codons will be represented by a hash (X).
#' * In protein alignments, sequence following the termination codon, will not be marked and will appear blank.
#' * These conventions are used for both nucleotide and protein alignments.
#'
#' This information about conventions was copied from the following URL:
#'
#' https://www.ebi.ac.uk/ipd/imgt/hla/alignment/help/
#'
#' @return A dataframe.
#' @param gene The name of a gene like "DRB1"
#' @param type The type of sequence, one of "prot", "nuc", "gen"
#' @param release The name of a release like "3.51.0"
#' @param quiet If FALSE, print messages along the way.
#' @examples
#' \donttest{
#' a <- hla_alignments("DRB1")
#' head(a$sequences)
#' head(a$aminos)
#' head(a$onehot)
#' }
#' @export
hla_alignments <- function(gene = "DRB1", type = "prot", release = NULL, quiet = TRUE) {
  hlabud_dir <- setup_hlabud_dir()
  tags_file <- file.path(hlabud_dir, "tags.json")
  releases <- hla_releases()
  if (is.null(release)) {
    release <- getOption("hlabud_release")
  }
  if (is.null(release)) {
    release <- releases[1]
    options(hlabud_release = release)
  }
  if (!release %in% releases) {
    stop("Unrecognized release '{release}' not in releases: {paste(releases, ' ')}")
  }
  if (!type %in% c("nuc", "gen", "prot")) {
    stop("Unrecognized type '{type}' not in: nuc gen prot")
  }
  if (!quiet) {
    message(glue("hlabud is using IMGTHLA release {release}"))
  }
  if (type == "prot") {
    prot_file <- file.path(hlabud_dir, release, "alignments", glue("{gene}_{type}.txt"))
    if (!file.exists(prot_file)) {
      ix <- which(releases == release)
      j <- read_json(tags_file)
      sha <- j[[ix]]$commit$sha
      repo_url <- "https://github.com/ANHIG/IMGTHLA"
      prot_url <- glue("{repo_url}/raw/{sha}/alignments/{gene}_{type}.txt")
      if (!quiet) {
        message(glue("Downloading {prot_url}"))
      }
      lines <- readLines(prot_url)
      mkdir(dirname(prot_file))
      if (!quiet) {
        message(glue("Writing {prot_file}"))
      }
      writeLines(lines, prot_file)
    }
    if (!quiet) {
      message(glue("Reading {prot_file}"))
    }
    return(read_prot(prot_file))
  }
  stop("not implemented yet")
}

#' Read a `*_prot.txt` file from IMGTHLA.
#'
#' @return A list with a dataframe and a matrix. The dataframe has two columns:
#' * allele: the name of the allele, e.g., `DQB*01:01`
#' * seq: the amino acid sequence
#' The matrix has a one-hot encoding of the variants among the alleles, with
#' one row for each allele and one column for each amino acid at each position.
#' @param prot_file File name for a txt file from IMGTHLA like "DQB1_prot.txt"
#' @examples
#' DRB1_file <- file.path(
#'   "https://github.com/ANHIG/IMGTHLA/raw",
#'   "5f2c562056f8ffa89aeea0631f2a52300ee0de17",
#'   "alignments/DRB1_prot.txt"
#' )
#' a <- read_prot(DRB1_file)
#' head(a$sequences)
#' head(a$aminos)
#' head(a$onehot)
#' @export
read_prot <- function(prot_file) {
  my_gene <- str_split_fixed(basename(prot_file), "_", 2)[,1]
  al <- readLines(prot_file)
  # Many amino acids are located before the position labeled as "1"
  #
  al_i <- which(str_detect(al, "Prot.+ 1$"))
  pre <- al[al_i]
  pre_i <- str_locate(pre, "-")[1,1]
  pre_j <- str_locate(pre, "1")[1,1]
  n_pre <- nchar(str_replace_all(substr(al[al_i + 2], pre_i, pre_j - 1), " ", ""))
  #
  my_regex <- glue("^ {my_gene}\\\\*")
  al <- al[str_detect(al, my_regex)]
  al <- str_split_fixed(al, " +", 3)
  al <- al[,2:3]
  colnames(al) <- c("allele", "seq")
  al <- as.data.frame(al)
  al <- al %>% group_by(.data$allele) %>% mutate(id = sprintf("V%s", seq(n()))) %>% ungroup()
  al <- al %>% pivot_wider(names_from = "id", values_from = "seq")
  al <- al %>% unite("seq", starts_with("V"), sep = "")
  al$seq <- str_replace_all(al$seq, " ", "")
  #
  # if (digits == 4) {
  #   # keep 4 digits
  #   my_allele_regex <- "[^: ]+:[^: ]+"
  # } else if (digits == 6) {
  #   # keep 6 digits
  #   my_allele_regex <- "[^: ]+:[^: ]+(:[^: ]+){0,1}"
  # } else {
  #   stop("digits must be 4 or 6")
  # }
  # al$d4 <- str_extract(al$allele, my_allele_regex)
  # al <- unique(al[,c("d4", "seq")])
  # al <- al %>% group_by(d4) %>% filter(row_number() == 1) %>%
  #   select(d4, allele, seq)
  oh <- get_onehot(al, n_pre)
  return(list(sequences = as.data.frame(al), aminos = oh$aminos, onehot = oh$onehot))
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
  #         P1=* P1=M P2=* P2=A P3=*
  # A*01:01    0    1    0    1    0
  # A*01:02    0    1    0    1    0
  # A*02:01    0    1    0    1    0
  # A*02:02    0    1    0    1    0
  # A*02:03    0    1    0    1    0
  # 
  #######################################################################
  aminos <- do.call(rbind, seq_chars)
  rownames(aminos) <- al$allele
  colnames(aminos) <- str_replace_all(
    sprintf("P%s", c(-n_pre:-1, 1:(ncol(aminos) - n_pre))), "-", "n"
  )
  # colnames(aminos) <- sprintf("P%s", seq(ncol(aminos)))
  # keep positions with more than 1 allele
  aminos <- aminos[,apply(aminos, 2, function(x) length(unique(x))) > 1, drop = FALSE]
  aminos <- as.data.frame(aminos)
  for (i in seq(ncol(aminos))) {
    aminos[,i] <- as.factor(aminos[,i])
  }
  p <- onehot(aminos, max_levels = 20)
  retval <- predict(p, aminos)
  rownames(retval) <- al$allele
  # Discard positions where we don't know the allele
  retval <- retval[,!str_detect(colnames(retval), "\\*")]
  colnames(retval) <- str_replace(colnames(retval), "=", "_")
  return(list(aminos = aminos, onehot = retval))
}

#' Amino acid dosage
#'
#' For each genotype, return the the dosage for each amino acid at each
#' position.
#'
#' Each genotype should be represented like `"HLA-A*01:01+HLA-A*01:01"`
#'
#' By default, the returned data frame is filtered to exclude:
#' * amino acid positions where all input genotypes have the same allele
#' * amino acid positions that are identical to previous positions
#'
#' @param genotypes Input character vector with one genotype for each individual.
#' @param aminos A one-hot encoded matrix with one row per allele and one
#' column per amino acid position.
#' @param drop_constants Filter out constant amino acid positions by default.
#' @param drop_duplicates Filter out duplicate amino acid positions by default.
#' @returns A data frame with one row for each input genotype.
#' @examples
#' DRB1_file <- file.path(
#'   "https://github.com/ANHIG/IMGTHLA/raw",
#'   "5f2c562056f8ffa89aeea0631f2a52300ee0de17",
#'   "alignments/DRB1_prot.txt"
#' )
#' a <- read_prot(DRB1_file)
#' genotypes <- c(
#'   "DRB1*12:02:02:03+DRB1*12:02:02:03+DRB1*14:54:02",
#'   "DRB1*04:174+DRB1*15:152",
#'   "DRB1*04:56:02+DRB1*15:01:48",
#'   "DRB1*14:172+DRB1*04:160",
#'   "DRB1*04:359+DRB1*04:284:02"
#' )
#' dosage <- amino_dosage(genotypes, a$aminos)
#' dosage[,1:5]
#' @export
amino_dosage <- function(genotypes, aminos, drop_constants = TRUE, drop_duplicates = TRUE) {
  dosages <- matrix(0, ncol = ncol(aminos), nrow = length(genotypes))
  for (i in seq_along(genotypes)) {
    # Split a string of genotypes like "HLA-A*01:01+HLA-A*01:01"
    a <- str_split(genotypes[i], "\\+")[[1]]
    for (my_a in a) {
      # Find the first row in aminos where the prefix matches our genotype
      ix <- which(str_starts(rownames(aminos), fixed(my_a)))
      if (length(ix) > 0) {
        dosages[i,] <- dosages[i,] + as.numeric(aminos[ix[1],])
      }
    }
  }
  rownames(dosages) <- genotypes
  colnames(dosages) <- colnames(aminos)
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

