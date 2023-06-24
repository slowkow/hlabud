library(stringr)

test_that("hla_alignment() works for all genes", {

  # load each of the hlabud dependencies
  # lapply(setdiff(unique(renv::dependencies()$Package), c("R")), function(x) library(x, character.only = TRUE))

  my_release <- "3.52.0"
  # my_gene <- "A"
  my_genes <- hla_genes(release = my_release)

  for (i in seq(nrow(my_genes))) {

    my_gene <- my_genes$gene[i]
    my_type <- my_genes$type[i]
    message(sprintf("gene %s  type %s", my_gene, my_type))

    expect_no_error({
      a <- hla_alignments(my_gene, type = my_type, release = my_release, verbose = TRUE)
    }, message = "^(Downloading|Writing|Reading)")
    expect_named(a, c("sequences", "alleles", "onehot"))
    expect_equal(nrow(a$sequences), nrow(a$alleles))
    expect_equal(nrow(a$sequences), nrow(a$onehot))

  }

})
