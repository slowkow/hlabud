
#' Convert one letter amino acid codes to three letter amino acid codes
#'
#' @keywords internal
one_to_three <- function(aminos) {
  dict <- c(
    "A" = "Ala", "C" = "Cys", "D" = "Asp", "E" = "Glu", "F" = "Phe", "G" = "Gly",
    "H" = "His", "I" = "Ile", "K" = "Lys", "L" = "Leu", "M" = "Met", "N" = "Asn",
    "P" = "Pro", "Q" = "Gln", "R" = "Arg", "S" = "Ser", "T" = "Thr", "V" = "Val",
    "W" = "Trp", "Y" = "Tyr", "*" = "Ter"
  )
  dict[aminos]
}

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
#' This function will:
#' - Get the folder name from `getOption("hlabud_dir")` or else automatically choose an appropriate folder for your operating system thanks to [rappdirs](https://github.com/r-lib/rappdirs).
#' - Create the folder automatically if it does not already exist.
#' - Set the the `hlabud_dir` option to that new folder.
#'
#' Here are the locations of the `hlabud_dir` folder on each operating system.
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
#'     C:\Documents and Settings\{User}\Application Data\slowkow\hlabud
#'
#' To set the `hlabud_dir` option, please use:
#' 
#'     options(hlabud_dir = "/my/favorite/path")
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

#' Download and unpack a tarball release from IMGTHLA
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
#' @seealso [hla_releases()] to get a complete list of all release names.
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

#' Get HLA gene names from IMGTHLA
#'
#' Retrieve the list of txt files in [`github.com/ANHIG/IMGTHLA/alignments`](https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments) and return a list of gene names derived from the file names.
#'
#' @param release Default is "latest". Should be a release name like "3.51.0".
#' @param overwrite Overwrite the existing `genes.json` file with a new one from GitHub
#' @param verbose If TRUE, print messages along the way.
#' @return A tibble with two columns: HLA gene names ("A", "DRB1") and types ("nuc", "gen", "prot").
#' @seealso [hla_releases()] to get a complete list of all release names.
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
  genes <- unique(str_extract(genes, "^[^.]+"))
  genes <- str_split_fixed(genes, "_", 2) %>% as.data.frame
  colnames(genes) <- c("gene", "type")
  # if (type) {
  #   genes <- unique(str_extract(genes, "^[^.]+"))
  # } else {
  #   genes <- unique(str_extract(genes, "^[^_]+"))
  # }
  # Exclude "ClassI", because it is an alignment of 3 HLA genes A, B, C
  # genes <- genes[!str_detect(genes, "ClassI")]
  genes <- genes[genes$gene != "ClassI",] %>% as_tibble
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
#' @seealso [hla_releases()] to get a complete list of all release names.
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

#' Get the names of releases from IMGTHLA
#'
#' Get tags from [github.com/ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA), save them to a file called `tags.json` in `getOption("hlabud_dir")`,  and return the release names from that file.
#'
#' The `tags.json` file will be automatically overwritten if it is older than 24 hours.
#'
#' @param overwrite Overwrite the existing tags.json file in `getOption("hlabud_dir")`
#' @return A character vector of release names like "3.51.0"
#' @examples
#' \donttest{
#' hla_releases() 
#' }
#' @export
hla_releases <- function(overwrite = FALSE) {
  hlabud_dir <- get_hlabud_dir()
  tags_file <- file.path(hlabud_dir, "tags.json")
  # Overwrite the tags.json file if it is older than a day.
  tags_old <- difftime(Sys.time(), file.info(tags_file)$mtime, units = "hours") > 24
  if (tags_old || overwrite || !file.exists(tags_file)) {
    j <- read_json("https://api.github.com/repos/ANHIG/IMGTHLA/tags")
    writeLines(toJSON(j, pretty = TRUE, auto_unbox = TRUE), tags_file)
  }
  j <- read_json(tags_file)
  release_names <- sapply(j, function(x) x$name)
  releases <- str_extract(release_names, "[\\d.]+")
  return(releases)
}

#' Get sequence alignments from IMGTHLA
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
#' @return A list with a dataframe called `sequences` and two matrices `alleles` and `onehot`.
#'
#' The dataframe has two columns:
#' * `allele`: the name of the allele, e.g., `DQB*01:01`
#' * `seq`: the amino acid sequence
#'
#' The matrix `alleles` has one row for each allele, and one column for each position, with the values representing the residues at each position in each allele.
#' The matrix `onehot` has a one-hot encoding of the variants that distinguish the alleles, with one row for each allele and one column for each amino acid at each position.
#' @seealso [hla_releases()] to get a complete list of all release names.
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
  # 
  genes_file <- file.path(hlabud_dir, release, "genes.json")
  if (!file.exists(genes_file)) {
    hla_genes(release = release)
  }
  # Get the file names from GitHub
  genes <- read_json(genes_file)
  file_names <- sapply(genes, "[[", "name")
  # Look for an exact match to the file we want
  i <- which(file_names == glue("{gene}_{type}.txt"))
  if (length(i) != 1) {
    message(glue("IMGT/HLA {release} does not provide {gene}_{type}.txt"))
    return()
  }
  genes <- genes[[i]]
  #
  my_name <- genes$name
  my_url <- genes$download_url
  my_type <- str_replace(my_name, "^.+_(\\w+)\\..+", "\\1")
  my_file <- file.path(hlabud_dir, release, "alignments", my_name)
  if (!file.exists(my_file)) {
    if (verbose) { message(glue("Downloading {my_url}")) }
    lines <- readLines(my_url)
    if (verbose) { message(glue("Writing {my_file}")) }
    mkdir(dirname(my_file))
    writeLines(lines, my_file)
  }
  if (verbose) { message(glue("Reading {my_file}")) }
  retval <- read_alignments(my_file)
  retval$file <- my_file
  retval$release <- release
  return(retval)
}

#' Read an alignment file `*_(nuc|gen|prot).txt` from IMGTHLA
#'
#' This function reads the txt files that are provided by IMGTHLA.
#'
#' Consider using [`hla_alignments()`] instead of this function. If you already have your own txt file that you want to read, then you can read it with `read_alignments("myfile.txt")`.
#'
#' These are the sequences contained in each file:
#' - `{gene}_prot.txt` has the amino acid sequence for each HLA allele.
#' - `{gene}_nuc.txt` has the nucleotide sequence for the exons.
#' - `{gene}_gen.txt` has the genomic sequence for the exons and introns.
#'
#' @return A list with a dataframe called `sequences` and two matrices `alleles` and `onehot`.
#'
#' The dataframe has two columns:
#' * `allele`: the name of the allele, e.g., `DQB*01:01`
#' * `seq`: the amino acid sequence
#'
#' The matrix `alleles` has one row for each allele, and one column for each position, with the values representing the residues at each position in each allele.
#' The matrix `onehot` has a one-hot encoding of the variants that distinguish the alleles, with one row for each allele and one column for each amino acid at each position.
#' @param file File name for a txt file from IMGTHLA like "DQB1_prot.txt"
#' @examples
#' my_file <- file.path(
#'   "https://github.com/ANHIG/IMGTHLA/raw",
#'   "5f2c562056f8ffa89aeea0631f2a52300ee0de17",
#'   "alignments/DRB1_prot.txt"
#' )
#' a <- read_alignments(my_file)
#' head(a$sequences)
#' a$alleles[1:5,1:5]
#' a$onehot[1:5,1:5]
#' @export
read_alignments <- function(file) {
  # Here is some ugly parsing code
  # Maybe one day we can switch to using a BNF parser, e.g. the {rly} R package
  my_type <- str_remove(str_split_fixed(basename(file), "_", 2)[,2], ".txt")
  my_gene <- str_split_fixed(basename(file), "_", 2)[,1]
  lines <- readLines(file)
  #
  if (my_type == "prot") {
    # Many amino acids are located before the position labeled as "1"
    al_i <- which(str_detect(lines, "Prot.+ 1$"))
    pre <- lines[al_i]
    pre_i <- str_locate(pre, "-")[1,1]
    if (is.na(pre_i)) {
      n_pre <- 0
    } else {
      pre_j <- str_locate(pre, "1")[1,1]
      n_pre <- nchar(str_replace_all(substr(lines[al_i + 2], pre_i, pre_j - 1), " ", ""))
    }
  } else if (my_type == "nuc") {
    al_i <- which(str_detect(lines, "cDNA"))
    n_pre <- as.numeric(str_remove(lines[al_i][1], " *cDNA +"))
    if (n_pre > 0) {
      n_pre <- 0
    } else {
      n_pre <- 3 * abs(n_pre)
    }
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
  # Convert to a matrix
  sequences <- as.matrix(al$seq)
  rownames(sequences) <- al$allele
  # Create a one-hot encoding
  oh <- get_onehot(al, n_pre)
  return(list(
    sequences = sequences,
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
  # Split the sequences into character vectors
  seq_chars <- str_split(al$seq, "")
  # The first sequence is the reference
  ref_chars <- seq_chars[[1]]
  # Not all sequences are the same length, so we need to know the range
  max_chars <- max(sapply(seq_chars, length))
  min_chars <- min(sapply(seq_chars, length))
  # Skip the reference (first sequence), and loop through the remaining sequences
  for (i in 2:length(seq_chars)) {
    # The "-" character indicates identity to the reference
    ix <- which(seq_chars[[i]] == "-")
    ix <- ix[ix <= length(ref_chars)]
    seq_chars[[i]][ix] <- ref_chars[ix]
    seq_len <- length(seq_chars[[i]])
    # If necessary, add "*" characters to pad the length of this sequence
    if (seq_len < max_chars) {
      seq_chars[[i]] <- c(seq_chars[[i]], rep("*", max_chars - seq_len))
    }
  }
  # If necessary, add "*" characters to pad the length of the reference sequence
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
  #
  # gaps "." do not increment the position
  gap <- alleles[1,] == "."
  #
  # deal with gaps of various lengths
  gap_i <- which(gap)
  nums <- cumsum(!gap) - n_pre
  lt1 <- nums < 1
  nums[lt1] <- nums[lt1] - 1
  nums[lt1] <- str_replace(nums[lt1], "-", "n")
  for (i in gap_i) {
    j <- i
    while (j < length(nums)) {
      j <- j + 1
      if (nums[i] != nums[j]) {
        nums[i] <- sprintf("%s_%s", nums[i], nums[j])
        break
      }
    }
  }
  nums <- sprintf("p%s", nums)
  # rbind(nums, alleles[1,])[,1:200]
  colnames(alleles) <- nums
  #
  # collapse the gap columns
  col_n <- table(colnames(alleles))
  col_n <- names(col_n[col_n > 1])
  for (my_col in col_n) {
    new_col <- apply(alleles[,colnames(alleles) == my_col], 1, paste, collapse = "")
    new_i <- which(colnames(alleles) == my_col)
    alleles[,new_i[1]] <- new_col
    alleles <- alleles[,-tail(new_i, length(new_i) - 1)]
  }
  #
  # keep positions with more than 1 allele
  keep_pos <- apply(alleles, 2, function(x) length(unique(x))) > 1
  if (!any(keep_pos)) {
    warning(glue("all positions have exactly 1 allele (there are no polymorphisms)"))
    return(list(alleles = as.matrix(alleles)))
  }
  alleles <- alleles[,keep_pos, drop = FALSE]
  alleles <- as.data.frame(alleles)
  for (i in seq(ncol(alleles))) {
    alleles[,i] <- as.factor(alleles[,i])
  }
  p <- onehot(alleles, max_levels = 20)
  retval <- predict(p, alleles)
  rownames(retval) <- al$allele
  colnames(retval) <- str_replace(colnames(retval), "=", "_")
  # Rename "*" to "unk" so we can use these names in formulas
  colnames(retval) <- str_replace(colnames(retval), "\\*", "unk")
  # Rename "." to "gap" so we can use these names in formulas
  colnames(retval) <- str_replace(colnames(retval), "\\.", "gap")
  #
  return(list(alleles = as.matrix(alleles), onehot = retval))
}

#' Convert a set of genotype names into a dosage matrix of each residue at each position
#'
#' For each genotype name, return the the dosage matrix for each residue (amino acid or nucleotide) at each position.
#'
#' Each genotype should be represented like this `"HLA-A*01:01,HLA-A*01:01"`
#'
#' By default, the returned matrix is filtered to exclude:
#' * positions where all input genotypes have the same allele
#'
#' @param mat A one-hot encoded matrix with one row per allele and one column for each residue (amino acid or nucleotide) at each position.
#' @param names Input character vector with one genotype for each individual. All entries must be present in `rownames(mat)`.
#' @param drop_constants Filter out constant amino acid positions. TRUE by default.
#' @param drop_duplicates Filter out duplicate amino acid positions. FALSE by default.
#' @return A matrix with one row for each input genotype, and one column for each residue at each position.
#' @examples
#' DRB1_file <- file.path(
#'   "https://github.com/ANHIG/IMGTHLA/raw",
#'   "5f2c562056f8ffa89aeea0631f2a52300ee0de17",
#'   "alignments/DRB1_prot.txt"
#' )
#' a <- read_alignments(DRB1_file)
#' genotypes <- c(
#'   "DRB1*12:02:02:03,DRB1*12:02:02:03,DRB1*14:54:02",
#'   "DRB1*04:174,DRB1*15:152",
#'   "DRB1*04:56:02,DRB1*15:01:48",
#'   "DRB1*14:172,DRB1*04:160",
#'   "DRB1*04:359,DRB1*04:284:02"
#' )
#' dosage <- dosage(a$onehot, genotypes)
#' dosage[,1:5]
#' @export
dosage <- function(mat, names, drop_constants = TRUE, drop_duplicates = FALSE) {
  dosages <- matrix(0, ncol = ncol(mat), nrow = length(names))
  real_names <- vector(mode = "character", length = length(names))
  for (i in seq_along(names)) {
    # Split a string of genotypes like "HLA-A*01:01,HLA-A*01:01"
    a <- str_split(names[i], ",")[[1]]
    for (my_a in a) {
      # Find the first row in mat where the prefix matches our genotype
      ix <- which(str_starts(rownames(mat), fixed(my_a)))
      if (length(ix) > 0) {
        real_names[i] <- rownames(mat)[ix[1]]
        dosages[i,] <- dosages[i,] + as.numeric(mat[ix[1],])
      } else {
        warning(glue("'{my_a}' not found in rownames"))
      }
    }
  }
  rownames(dosages) <- real_names
  colnames(dosages) <- colnames(mat)
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

