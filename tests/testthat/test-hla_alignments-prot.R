library(stringr)
library(glue)

my_release <- "3.52.0"
my_genes <- hla_genes(release = my_release)

for (i in seq(nrow(my_genes))) {
  my_gene <- my_genes$gene[i]
  my_type <- my_genes$type[i]
  if (my_type != "prot") {
    next
  }
  test_that(glue('hla_alignment("{my_gene}", type = "{my_type}", release = "{my_release}")'), {
    expect_no_error({
      a <- hla_alignments(my_gene, type = my_type, release = my_release, verbose = FALSE)
    }, message = "^(Downloading|Writing|Reading)")
    expect_named(a, c("sequences", "alleles", "onehot", "gene", "type", "release", "file"))
    expect_equal(length(a$sequences), nrow(a$alleles))
    expect_equal(length(a$sequences), nrow(a$onehot))
  })
}
