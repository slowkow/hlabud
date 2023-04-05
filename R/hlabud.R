
library(jsonlite)
library(rappdirs)
library(glue)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(magrittr)
library(tibble)
library(onehot)

mkdir <- function(path) {
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

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

# Download and unpack a tarball release of IMGTHLA from GitHub
#
install_hla <- function(release = "latest") {

  hlabud_dir <- setup_hlabud_dir()

  message(glue("Fetching releases from GitHub"))
  j <- read_json("https://api.github.com/repos/ANHIG/IMGTHLA/tags")
  writeLines(toJSON(j, pretty = TRUE, auto_unbox = TRUE), file.path(hlabud_dir, "tags.json"))

  release_names <- sapply(j, function(x) x$name)
  releases <- str_extract(release_names, "[\\d.]+")

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
    download.file(url, destfile = tar_file)
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

hla_alignments <- function(release = NULL, gene = "DRB", type = "prot") {
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
  if (type == "prot") {
    prot_file <- file.path(hlabud_dir, release, "alignments", glue("{gene}_{type}.txt"))
    if (!file.exists(prot_file)) {
      ix <- which(releases == release)
      j <- read_json(tags_file)
      sha <- j[[ix]]$commit$sha
      repo_url <- "https://github.com/ANHIG/IMGTHLA"
      lines <- readLines(glue("{repo_url}/raw/{sha}/alignments/{gene}_{type}.txt"))
      mkdir(dirname(prot_file))
      writeLines(lines, prot_file)
    }
    return(read_prot(prot_file))
  }
  stop("not implemented yet")
}

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
  al <- as.data.table(al)
  al$V1 <- NULL
  colnames(al) <- c("allele", "seq")
  al <- as_tibble(al)
  al <- al %>% group_by(allele) %>% mutate(id = sprintf("V%s", seq(n()))) %>% ungroup()
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
  return(list(sequences = al, onehot = get_onehot(al, n_pre)))
}

# al has two columns: allele, seq
get_onehot <- function(al, n_pre) {
  seq_chars <- str_split(al$seq, "")
  ref_chars <- seq_chars[[1]]
  max_chars <- max(sapply(seq_chars, length))
  for (i in 2:length(seq_chars)) {
    ix <- which(seq_chars[[i]] == "-")
    ix <- ix[ix < length(ref_chars)]
    seq_chars[[i]][ix] <- ref_chars[ix]
    if (length(seq_chars[[i]]) < max_chars) {
      seq_chars[[i]] <- c(
        seq_chars[[i]],
        rep(NA, max_chars - length(seq_chars[[i]]))
      )
    }
  }
  if (length(seq_chars[[1]]) < max_chars) {
    seq_chars[[1]] <- c(
      seq_chars[[1]],
      rep(NA, max_chars - length(seq_chars[[1]]))
    )
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
  aminos <- predict(p, aminos)
  rownames(aminos) <- al$allele
  # Discard positions where we don't know the allele
  aminos <- aminos[,!str_detect(colnames(aminos), "\\*")]
  colnames(aminos) <- str_replace(colnames(aminos), "=", "_")
  return(aminos)
}

amino_dosage <- function(alleles, aminos) {
  dosages <- matrix(NA, ncol = ncol(aminos), nrow = length(alleles))
  for (i in seq_along(alleles)) {
    a <- str_split(alleles[i], "\\+")[[1]]
    if (!all(a %in% rownames(aminos))) {
      stop(glue("allele not found in rownames(aminos): {a}"))
    }
    dosages[i,] <- colSums(aminos[a,], na.rm = TRUE)
  }
  rownames(dosages) <- alleles
  colnames(dosages) <- colnames(aminos)
  # Select positions where we observe more than 1 possible dosage [0, 1, 2]
  ix <- apply(dosages, 2, function(x) length(unique(x)))
  dosages <- dosages[, ix > 1, drop = FALSE]
  if (!any(ix)) {
    return(dosages)
  }
  # Discard positions that are identical to previous positions
  ix <- !duplicated(apply(dosages, 2, function(x) paste(x, collapse = ",")))
  dosages <- dosages[,which(ix)]
  return(dosages)
}

# # libraries, functions {{{
# 
# library(pacman)
# pacman::p_load(
#   conflicted,
#   onehot,
#   data.table,
#   ggbeeswarm,
#   ggforce,
#   ggplot2,
#   ggstance,
#   ggpmisc,
#   ggrepel,
#   ggtext,
#   dplyr,
#   ggforestplot,
#   glue,
#   janitor,
#   magrittr,
#   naturalsort,
#   pals,
#   patchwork,
#   parameters,
#   pbapply,
#   qs,
#   scales,
#   scico,
#   stringr,
#   tidyr
# )
# # devtools::install_github("NightingaleHealth/ggforestplot")
# # source("R/batch_1-4_5p/helpers.R")
# source("functions/helpers.R")
# source("functions/mpn65.R")
# source("functions/palettes.R")
# source("functions/theme-kamil.R")
# theme_set(theme_kamil)
# conflict_prefer("filter", "dplyr")
# conflict_prefer("select", "dplyr")
# conflict_prefer("mutate", "dplyr")
# 
# # }}}
# 
# # HLA genotypes
# ########################################################################
# hla <- fread("/mnt/covid-mad/output-neutrophils/arcasHLA/genotypes.tsv")
# 
# # hla-analysis/amino/positions {{{
# 
# 
# alleles_to_aminos <- function(gene = "DQB1", alleles) {
# 
#   # my_gene <- "DQB1"
#   my_gene <- gene
# 
#   # keep 4 digits or 6 digits?
#   my_allele_regex <- "[^: ]+:[^: ]+"
#   # my_allele_regex <- "[^: ]+:[^: ]+(:[^: ]+){0,1}"
# 
#   # Get the amino acid sequences for alleles of this gene
#   # my_gene <- "DQB1"
#   my_regex <- glue("^ {my_gene}\\\\*")
#   if (my_gene == "DRB1") {
#     al <- readLines(glue("data/alleles.org/DRB_prot.txt"))
#   } else {
#     al <- readLines(glue("data/alleles.org/{my_gene}_prot.txt"))
#   }
# 
#   # Many amino acids are located before the position labeled as "1"
#   #
#   al_i <- which(str_detect(al, "Prot.+ 1$"))
#   pre <- al[al_i]
#   pre_i <- str_locate(pre, "-")[1,1]
#   pre_j <- str_locate(pre, "1")[1,1]
#   n_pre <- nchar(str_replace_all(substr(al[al_i + 2], pre_i, pre_j - 1), " ", ""))
#   #
#   al <- al[str_detect(al, my_regex)]
#   al <- str_split_fixed(al, " +", 3)
#   al <- as.data.table(al)
#   al$V1 <- NULL
#   colnames(al) <- c("allele", "seq")
#   al <- as_tibble(al)
#   al <- al %>% group_by(allele) %>% mutate(id = sprintf("V%s", seq(n()))) %>% ungroup()
#   al <- al %>% pivot_wider(names_from = "id", values_from = "seq")
#   al <- al %>% unite("seq", starts_with("V"), sep = "")
#   al$seq <- str_replace_all(al$seq, " ", "")
#   #
#   # Keep 4 digits
#   al$d4 <- str_extract(al$allele, my_allele_regex)
#   al <- unique(al[,c("d4", "seq")])
#   al <- al %>% group_by(d4) %>% filter(row_number() == 1)
#   #
#   # Get the alleles for this gene
#   my_key <- sprintf("%s_", my_gene)
#   my_hla <- hla %>%
#     # filter(covid == 1, time_point == "D0") %>%
#     filter(time_point == "D0") %>%
#     select(
#       sample, donor, neut, crp_0_cat, crp_0_n, covid, acuity_maxn, acuity_max_cat,
#       cd4_clonality, cd8_clonality, bcr_clonality, aa, vl, il6, ctrc,
#       starts_with(!!my_key)
#     ) %>%
#     pivot_longer(-c(
#       sample, donor, neut, crp_0_cat, crp_0_n, covid, acuity_maxn, acuity_max_cat,
#       cd4_clonality, cd8_clonality, bcr_clonality, aa, vl, il6, ctrc
#     ))
#   #
#   keep_donors <- (
#     my_hla %>%
#     group_by(donor) %>%
#     summarize(n_alleles = length(unique(value)), .groups = "keep") %>%
#     arrange(n_alleles) %>%
#     filter(n_alleles <= 2)
#   )$donor
#   length(keep_donors)
#   #
#   my_hla <- my_hla %>%
#     filter(donor %in% keep_donors) %>%
#     select(
#       donor, neut, crp_0_cat, crp_0_n, covid, acuity_maxn, acuity_max_cat,
#       cd4_clonality, cd8_clonality, bcr_clonality, aa, vl, il6, ctrc,
#       name, value
#     ) %>%
#     filter(!is.na(value)) %>%
#     unique()
#   #
#   # Keep first few digits
#   my_hla$value <- str_extract(my_hla$value, my_allele_regex)
#   #
#   al <- al[al$d4 %in% my_hla$value,]
#   #
#   seq_chars <- str_split(al$seq, "")
#   ref_chars <- seq_chars[[1]]
#   max_chars <- max(sapply(seq_chars, length))
#   for (i in 2:length(seq_chars)) {
#     ix <- which(seq_chars[[i]] == "-")
#     seq_chars[[i]][ix] <- ref_chars[ix]
#     if (length(seq_chars[[i]]) < max_chars) {
#       seq_chars[[i]] <- c(
#         seq_chars[[i]],
#         rep(NA, max_chars - length(seq_chars[[i]]))
#       )
#     }
#   }
#   setdiff(my_hla$value, al$d4)
#   #
#   # Create a one-hot-encoded matrix
#   # with allele names (rows) and positions (columns)
#   # 
#   #         P1=* P1=M P2=* P2=A P3=*
#   # A*01:01    0    1    0    1    0
#   # A*01:02    0    1    0    1    0
#   # A*02:01    0    1    0    1    0
#   # A*02:02    0    1    0    1    0
#   # A*02:03    0    1    0    1    0
#   # 
#   #######################################################################
#   aminos <- do.call(rbind, seq_chars)
#   rownames(aminos) <- al$d4
#   colnames(aminos) <- str_replace_all(
#     sprintf("P%s", c(-n_pre:-1, 1:(ncol(aminos) - n_pre))), "-", "n"
#   )
#   # colnames(aminos) <- sprintf("P%s", seq(ncol(aminos)))
#   # keep positions with more than 1 allele
#   aminos <- aminos[,apply(aminos, 2, function(x) length(unique(x))) > 1, drop = FALSE]
#   aminos <- as.data.frame(aminos)
#   for (i in seq(ncol(aminos))) {
#     aminos[,i] <- as.factor(aminos[,i])
#   }
#   p <- onehot(aminos)
#   aminos <- predict(p, aminos)
#   rownames(aminos) <- al$d4
#   # Discard positions where we don't know the allele
#   aminos <- aminos[,!str_detect(colnames(aminos), "\\*")]
#   colnames(aminos) <- str_replace(colnames(aminos), "=", "_")
#   aminos_file <- glue("{out_dir}/amino/table/aminos/{my_gene}.tsv")
#   mkdir(dirname(aminos_file))
#   message(glue("Writing {aminos_file}"))
#   fwrite(
#     as.data.frame(aminos) %>% tibble::rownames_to_column("alelle"),
#     aminos_file,
#     sep = "\t"
#   )
#   #
#   d <- my_hla %>%
#       group_by(
#         cd4_clonality, cd8_clonality, bcr_clonality, aa, vl, il6, ctrc,
#         neut, crp_0_cat, crp_0_n, covid, acuity_maxn, acuity_max_cat, donor
#       ) %>%
#       summarize(
#           alleles = paste(sort(unique(value)), collapse = ","),
#           .groups = "drop"
#       )
#   dosages <- matrix(NA, ncol = ncol(aminos), nrow = nrow(d))
#   for (i in seq(nrow(d))) {
#     alleles <- str_split(d$alleles[i], ",")[[1]]
#     if (length(alleles) == 1) {
#       alleles <- c(alleles, alleles)
#     }
#     dosages[i,] <- colSums(aminos[alleles,], na.rm = TRUE)
#   }
#   colnames(dosages) <- colnames(aminos)
#   # Select positions where we have observed all possible dosages [0, 1, 2]
#   dosages <- dosages[,apply(dosages, 2, function(x) length(unique(x)) == 3),drop=FALSE]
#   if (ncol(dosages) == 0) {
#     next
#   }
#   #
#   # There are haplotype blocks, so some positions are identical to each other.
#   ix <- duplicated(apply(dosages, 2, function(x) paste(x, collapse = ",")))
#   my_is <- which(!ix)
#   my_alleles <- colnames(dosages)[my_is]
#   #
#   d <- cbind(d, dosages)
# 
#   d_file <- glue("{out_dir}/amino/table/dosage/dosage-{my_gene}.tsv")
#   mkdir(dirname(d_file))
#   message(glue("Writing {d_file}"))
#   fwrite(d, d_file, sep = "\t")
# 
#   # Number of donors with each classical allele, e.g., "DQB1*03:01"
#   # d %>% count(alleles, P13_A, acuity_max_cat) %>% arrange(-P13_A) %>% head(20)
#   d2_file <- glue("{out_dir}/amino/table/alleles/alleles-{my_gene}.tsv")
#   mkdir(dirname(d2_file))
#   message(glue("Writing {d2_file}"))
#   fwrite(
#     d %>% select(donor, alleles, acuity_max_cat),
#     d2_file, sep = "\t"
#   )
# 
#   # logistic regression with death (acuity_max=1) {{{
#   # plot these results
#   my_glm <- rbindlist(pblapply(unique(my_alleles), function(my_a) {
#     f <- sprintf("death ~ %s", my_a)
#     glm(
#       as.formula(f), family = "binomial",
#       # data = d %>% mutate(acuity_maxn = as.integer(rescale(acuity_maxn) >= 0.75))
#       data = d %>% mutate(death = as.integer(rescale(acuity_maxn) >= 1))
#     ) %>% parameters(exponentiate = TRUE) %>% mutate(allele = my_a)
#   })) %>% mutate(fdr = p.adjust(p, method = "fdr"))
#   my_glm <- clean_names(my_glm)
#   my_file <- glue("{out_dir}/amino/table/models/{my_coef}/{my_coef}-{my_gene}-glm.tsv")
#   dir.create(dirname(my_file), recursive = TRUE, showWarnings = FALSE)
#   message(glue("Writing {my_file}"))
#   fwrite(x = my_glm, file = my_file, sep = "\t")
# 
#   my_d <- my_glm %>% filter(fdr < 0.1) %>% filter(parameter != "(Intercept)")
#   p <- ggplot(my_d) +
#     aes(x = coefficient, y = parameter) +
#     ggforestplot::geom_stripes(odd = "grey95", even = "#00000000") +
#     geom_vline(xintercept = 1, linewidth = 0.3, alpha = 0.3) +
#     geom_errorbar(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.5, width = 0, alpha = 0.3) +
#     geom_point(aes(fill = fdr < 0.05), shape = 21, size = 3, stroke = 0.3) +
#     scale_fill_manual(
#       name = glue("FDR < 0.05"),
#       values = c("TRUE" = "grey20", "FALSE" = "grey95")
#     ) +
#     scale_x_continuous(
#       breaks = pretty_breaks(4),
#       # labels = function(x) signif(2^x, 2)
#     ) +
#     theme(
#       plot.subtitle   = element_textbox(width = unit(1, "npc")),
#       legend.position = "top",
#       strip.text.y    = element_blank(),
#       panel.spacing   = unit(1, "lines")
#     ) +
#     labs(
#       title = glue("{my_gene} association with death"),
#       caption = "death ~ dosage",
#       y = "Allele",
#       x = "Odds Ratio"
#     )
#   my_width <- 5
#   my_height <- 1 + length(unique(my_d$parameter)) * 0.4
#   my_ggsave(
#     slug = glue("death-errorbars-{my_gene}"),
#     out_dir = glue("{out_dir}/amino/glm"),
#     plot = p,
#     type = "pdf",
#     scale = 1, width = my_width, height = my_height, units = "in", dpi = 300
#   )
#   my_d <- my_glm %>% filter(parameter != "(Intercept)")
#   p <- ggplot(my_d) +
#     aes(x = coefficient, y = -log10(p)) +
#     # ggplot2::geom_errorbarh(
#     #   data = function(d) d %>% filter(fdr < 0.05),
#     #   aes(xmin = CI_low, xmax = CI_high)
#     # ) +
#     geom_vline(xintercept = 1, linewidth = 0.3, alpha = 0.3) +
#     geom_point(aes(fill = fdr < 0.05), shape = 21, size = 3, stroke = 0.3) +
#     geom_text_repel(data = function(d) d %>% filter(fdr < 0.05), mapping = aes(label = allele)) +
#     scale_fill_manual(
#       name = NULL,
#       values = c("TRUE" = "grey20", "FALSE" = "grey95"),
#       labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR >= 0.05")
#     ) +
#     scale_x_continuous(
#       breaks = pretty_breaks(4),
#       expand = expansion(mult = c(0.1, 0.1))
#     ) +
#     scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
#     labs(x = "OR", y = "-log10 P", title = glue("{my_gene} association with death"), caption = "death ~ dosage") +
#     theme(legend.position = "top")
#   my_width <- 5
#   my_height <- 5.5
#   my_ggsave(
#     slug = glue("death-volcano-{my_gene}"),
#     out_dir = glue("{out_dir}/amino/glm"),
#     plot = p,
#     type = "pdf",
#     scale = 1, width = my_width, height = my_height, units = "in", dpi = 300
#   )
#   # }}}
# 
#   # # ordinal logistic regression with MASS::polr() {{{
#   # # plot these results for top alleles
#   # m <- MASS::polr(acuity_max_cat ~ P13_A, data = d %>% mutate(acuity_max_cat = fct_rev(acuity_max_cat)))
#   # # m <- polr(acuity_max_cat ~ P13_A, data = d)
#   # # p <- predict(m, d)
#   # # table(d$acuity_max_cat, p)
#   # parameters(m, exponentiate = TRUE)
#   # # }}}
# 
#   # linear regression with lm() {{{
#   my_coefs <- c(
#     "acuity_maxn", "neut", "cd8_clonality", "cd4_clonality", "bcr_clonality",
#     "aa", "vl", "il6", "ctrc"
#   )
#   for (my_coef in my_coefs) {
# 
#     # null model {{{
#     if (FALSE) {
#       set.seed(42)
#       null_models <- pblapply(seq_len(2000), function(i) {
#           # lapply(my_alleles, function(my_allele) {
#           lapply(c("P45_A"), function(my_allele) {
#           if (my_coef == "neut") {
#             my_form <- as.formula(glue("{my_coef} ~ {my_allele} + crp_0_n"))
#           } else {
#             my_form <- as.formula(glue("{my_coef} ~ {my_allele}"))
#           }
#           d_null <- d
#           ix <- sample(seq_len(nrow(d_null)))
#           d_null[[my_coef]] <- d_null[[my_coef]][ix]
#           retval <- parameters(lm(my_form, d_null))
#           retval$Parameter[retval$Parameter == my_allele] <- "dosage"
#           retval$allele <- my_allele
#           retval
#         })
#       })
#       null_ps <- sapply(seq_along(null_models), \(i) null_models[[i]][[1]]$p[2])
# 
#       null_hist <- hist(null_points)
# 
#       # Looks good!
#       my_breaks <- c(0.005, 0.05, 0.1, 0.2, 0.5)
#       null_points <- sapply(
#         my_breaks,
#         \(x) (sum(null_ps <= x) + 1) / (length(null_ps) + 1)
#       )
#       # [1] 0.04795205 0.09690310 0.20479520 0.50649351
#       null_d <- data.frame(exp = -log10(my_breaks), obs = -log10(null_points))
#       ci <- 0.95
#       n <- nrow(null_d)
#       null_d$clower <- -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1))
#       null_d$cupper <- -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
#       p <- ggplot(null_d) +
#         aes(x = exp, y = obs) +
#         scale_x_continuous(
#           trans = "log10",
#           breaks = -log10(my_breaks),
#           labels = \(x) 10 ^ -x
#         ) +
#         scale_y_continuous(
#           trans = "log10",
#           breaks = -log10(my_breaks),
#           labels = \(x) 10 ^ -x
#         ) +
#         geom_abline(intercept = 0, slope = 1, size = 0.3, alpha = 0.3) +
#         # geom_ribbon(
#         #   mapping = aes(x = exp, ymin = clower, ymax = cupper),
#         #   alpha = 0.1
#         # ) +
#         geom_line() +
#         geom_point() +
#         theme(
#           panel.grid.major = element_line()
#         ) +
#       labs(
#         title = "Distribution of p-values",
#         subtitle = "(n = 1000 permutations)",
#         x = "Expected p-value", y = "Observed p-value"
#       )
#       my_ggsave(
#         slug = glue("{my_coef}__{my_gene}_null_pvalues"),
#         out_dir = file.path(out_dir, "amino", "positions"),
#         type = "pdf",
#         plot = p,
#         scale = 1, units = "in", dpi = 300,
#         width = 4,
#         height = 4
#       )
# 
#       null_d <- data.frame(p = null_ps)
#       p <- ggplot(null_d) +
#         aes(x = p) +
#         geom_histogram(bins = 5)
#       my_ggsave(
#         slug = glue("{my_coef}__{my_gene}_null_pvalues"),
#         out_dir = file.path(out_dir, "amino", "positions"),
#         type = "pdf",
#         plot = p,
#         scale = 1, units = "in", dpi = 300,
#         width = 5,
#         height = 4
#       )
# 
#       p <- gg_qqplot(null_ps)
#       my_ggsave(
#         slug = glue("{my_coef}__{my_gene}_null_pvalues"),
#         out_dir = file.path(out_dir, "amino", "positions"),
#         type = "pdf",
#         plot = p,
#         scale = 1, units = "in", dpi = 300,
#         width = 4,
#         height = 3.5
#       )
#     }
#     # }}}
# 
#     models <- lapply(my_alleles, function(my_allele) {
#       if (my_coef == "neut") {
#         my_form <- as.formula(glue("{my_coef} ~ {my_allele} + crp_0_n"))
#       } else {
#         my_form <- as.formula(glue("{my_coef} ~ {my_allele}"))
#       }
#       retval <- parameters(lm(my_form, d))
#       retval$Parameter[retval$Parameter == my_allele] <- "dosage"
#       retval$allele <- my_allele
#       retval
#     })
#     models <- as.data.table(do.call(rbind, models))
#     models <- clean_names(models)
#     models %<>%
#       group_by(parameter) %>%
#       mutate(fdr = p.adjust(p, method = "fdr")) %>%
#       ungroup()
#     models$gene <- my_gene
#     my_file <- glue("{out_dir}/amino/table/models/{my_coef}/{my_coef}-{my_gene}.tsv")
#     dir.create(dirname(my_file), recursive = TRUE, showWarnings = FALSE)
#     message(glue("Writing {my_file}"))
#     fwrite(x = models, file = my_file, sep = "\t")
# 
#     models %>%
#       select(allele, parameter, coefficient, ci_low, ci_high, p, fdr) %>%
#       filter(parameter == "dosage") %>%
#       arrange(p)
#       # filter(allele == "P72_A")
# 
#     # mean(d$P13_A[d$acuity_max_cat == "1 Death"]) - mean(d$P13_A[d$acuity_max_cat == "3 SupO2"])
#     # mean(d$P13_G[d$acuity_max_cat == "1 Death"]) - mean(d$P13_G[d$acuity_max_cat == "3 SupO2"])
# 
#     my_top_alleles <- (
#       models %>%
#         filter(parameter == "dosage") %>%
#         top_n(n = 3, wt = -p)
#     )$allele
#     for (my_a in my_top_alleles) {
#       formula <- y ~ as.numeric(x)
#       p <- d %>% select(!!sym(my_a), acuity_max_cat, acuity_maxn) %>%
#         mutate(acuity_maxn = as.integer(rescale(acuity_maxn) == 1)) %>%
#         ggplot() +
#         aes(x = !!sym(my_a), y = acuity_maxn) +
#         geom_point(aes(fill = acuity_max_cat),
#           size = 3, shape = 21, stroke = 0.3, position = position_jitter()) +
#         scale_x_continuous(breaks = c(0, 1, 2)) +
#         # geom_smooth(method = lm, formula = formula)
#         stat_poly_line(formula = formula) +
#         stat_poly_eq(mapping = use_label("eq"), formula = formula, label.x = "right", label.y = "bottom")
#       my_ggsave(
#         slug = glue("{my_coef}-{my_gene}-{my_a}"),
#         out_dir = file.path(out_dir, "amino", "lm", my_coef),
#         type = "pdf",
#         plot = p,
#         scale = 1, units = "in", dpi = 300,
#         width = 5,
#         height = 4
#       )
#     }
# 
#     d_p <- models %>%
#       filter(parameter == "dosage") %>%
#       arrange(p)
#     d_p$position <- as.integer(str_extract(d_p$allele, "\\d+"))
#     #
#     p <- ggplot(d_p) +
#         aes(y = -log10(p), x = position, fill = fdr < 0.05, size = fdr < 0.05) +
#         geom_hline(
#             linetype = 2,
#             yintercept = -log10(0.05 / nrow(d_p)), linewidth = 0.3, color = "grey50") +
#         geom_point(shape = 21, stroke = 0.3) +
#         scale_size_manual(values = c(1, 2)) +
#         scale_fill_manual(
#             values = c("grey80", "grey20"),
#             guide = guide_legend(override.aes = list(label = "", size = 3))
#         ) +
#         theme(
#             legend.position = "top",
#             legend.box.spacing = unit(0, "lines")
#         ) +
#         labs(
#             title = glue("{my_gene} amino acid associations with {my_coef}"),
#             y = "-log10 P", x = "Amino Acid Position"
#         )
#     my_ggsave(
#       slug = glue("{my_coef}-{my_gene}"),
#       out_dir = file.path(out_dir, "amino", "positions", my_coef),
#       type = "pdf",
#       plot = p,
#       scale = 1, units = "in", dpi = 300,
#       width = 5,
#       height = 4
#     )
# 
#     if (my_coef == "acuity_maxn") {
#       for (my_allele in my_top_alleles) {
#           my_d <- d %>%
#               select(donor, alleles, acuity_max_cat, allele = !!my_allele) %>%
#               count(acuity_max_cat, allele)
#           my_p <- d_p %>% filter(allele == my_allele)
#           if (my_p$fdr < 0.05) {
#             p <- ggplot(my_d) +
#                 aes(y = acuity_max_cat, x = n, fill = acuity_max_cat) +
#                 geom_colh() +
#                 scale_fill_manual(values = palettes$acuity_max_cat) +
#                 scale_x_continuous(breaks = pretty_breaks(2), labels = round) +
#                 scale_y_discrete(position = "right") +
#                 facet_wrap(~ allele, ncol = 3, scales = "free_x") +
#                 theme(
#                     legend.position = "none",
#                     panel.spacing = unit(1, "lines")
#                 ) +
#                 labs(
#                     title = glue("acuity_max by copies of {my_gene} {my_allele}"),
#                     subtitle = glue("P = {signif(my_p$p, 2)}, FDR = {signif(my_p$fdr, 2)}"),
#                     y = NULL, x = "Donors"
#                 )
#             my_ggsave(
#               slug = glue("{my_gene}__{my_allele}"),
#               out_dir = file.path(out_dir, "amino", "alleles", my_coef),
#               type = "pdf",
#               plot = p,
#               scale = 1, units = "in", dpi = 300,
#               width = 5,
#               height = 3
#             )
#         }
#       }
#     }
# 
#     if (my_coef == "neut") {
#       for (my_allele in my_top_alleles) {
#           my_d <- d %>%
#               select(
#                 donor, alleles, neut, crp_0_cat, crp_0_n, covid, acuity_max_cat,
#                 allele = !!my_allele
#               ) %>%
#               mutate(allele = naturalfactor(allele))
#           my_p <- d_p %>% filter(allele == my_allele)
#           if (my_p$fdr < 0.05) {
# 
#             n_colors <- length(unique(my_d$crp_0_cat))
#             p <- ggplot(my_d) +
#                 aes(x = allele, y = neut) +
#                 geom_point(
#                   # aes(color = as.factor(covid)),
#                   aes(fill = crp_0_cat),
#                   shape = 21,
#                   size = 2, stroke = 0.3,
#                   position = position_quasirandom(dodge.width = 0.7)
#                 ) +
#                 geom_boxplot(
#                   alpha = 0.0, coef = 0, outlier.shape = NA,
#                   size = 0.3
#                 ) +
#                 scale_fill_manual(
#                   values = scico::scico(
#                     n = n_colors + 1, palette = "tokyo",
#                     direction = -1
#                   )[2:(n_colors+1)],
#                   na.value = "white"
#                 ) +
#                 guides(
#                   fill = guide_legend(override.aes = list(size = 3), reverse = TRUE)
#                 ) +
#                 scale_y_continuous(breaks = pretty_breaks(2)) +
#                 theme(
#                   # legend.position = "none",
#                   # legend.position = "top",
#                   legend.position = "right",
#                   panel.spacing = unit(1, "lines"),
#                   legend.box.spacing = unit(0, "lines")
#                 ) +
#                 labs(
#                   title = glue("{my_gene} {my_allele}"),
#                   subtitle = glue("P = {signif(my_p$p, 2)}, FDR = {signif(my_p$fdr, 2)}"),
#                   y = "Neutralization", x = "Genotype"
#                 )
#             my_ggsave(
#               slug = glue("{my_gene}__{my_allele}"),
#               out_dir = file.path(out_dir, "amino", "alleles", my_coef),
#               type = "pdf",
#               plot = p,
#               scale = 1, units = "in", dpi = 300,
#               width = 5,
#               height = 4
#             )
# 
#         }
#       }
#     }
# 
#     if (my_coef == "bcr_clonality") {
#       for (my_allele in my_top_alleles) {
#           my_d <- d %>%
#               select(
#                 donor, alleles, bcr_clonality, covid, acuity_max_cat,
#                 allele = !!my_allele
#               ) %>%
#               mutate(allele = naturalfactor(allele))
#           my_p <- d_p %>% filter(allele == my_allele)
#           if (my_p$fdr < 0.05) {
# 
#             p <- ggplot(my_d) +
#                 aes(x = allele, y = bcr_clonality) +
#                 geom_point(
#                   # aes(color = as.factor(covid)),
#                   # aes(fill = crp_0_cat),
#                   shape = 21,
#                   size = 2, stroke = 0.3,
#                   position = position_quasirandom()
#                 ) +
#                 geom_boxplot(
#                   alpha = 0.0, coef = 0, outlier.shape = NA,
#                   size = 0.3
#                 ) +
#                 # scale_fill_manual(
#                 #   values = scico::scico(
#                 #     n = n_colors + 1, palette = "tokyo",
#                 #     direction = -1
#                 #   )[2:(n_colors+1)],
#                 #   na.value = "white"
#                 # ) +
#                 # guides(
#                 #   fill = guide_legend(override.aes = list(size = 3), reverse = TRUE)
#                 # ) +
#                 scale_y_continuous(breaks = pretty_breaks(2)) +
#                 theme(
#                   # legend.position = "none",
#                   # legend.position = "top",
#                   # legend.position = "right",
#                   panel.spacing = unit(1, "lines"),
#                   legend.box.spacing = unit(0, "lines")
#                 ) +
#                 labs(
#                   title = glue("{my_gene} {my_allele}"),
#                   subtitle = glue("P = {signif(my_p$p, 2)}, FDR = {signif(my_p$fdr, 2)}"),
#                   y = "BCR Clonality", x = "Genotype"
#                 )
#             my_ggsave(
#               slug = glue("{my_gene}__{my_allele}"),
#               out_dir = file.path(out_dir, "amino", "alleles", my_coef),
#               type = "pdf",
#               plot = p,
#               scale = 1, units = "in", dpi = 300,
#               width = 4,
#               height = 4
#             )
# 
#         }
#       }
#     }
# 
#   }
#   # }}}
# 
# }
# 
