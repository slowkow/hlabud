library(stringr)
library(glue)

# load each of the hlabud dependencies
# lapply(setdiff(unique(renv::dependencies()$Package), c("R")), function(x) library(x, character.only = TRUE))

my_release <- "3.52.0"
# my_gene <- "A"
# my_type <- "prot"
my_genes <- hla_genes(release = my_release)

# my_gene <- "DRB1"
# my_type <- "gen"

for (i in seq(nrow(my_genes))) {
  my_gene <- my_genes$gene[i]
  my_type <- my_genes$type[i]
  test_that(glue("hla_alignment({my_gene}, type = {my_type}, release = {my_release})"), {
      if (my_gene == "N" && my_type == "nuc") {
        expect_warning({
          a <- hla_alignments(my_gene, type = my_type, release = my_release, verbose = TRUE)
        }, regexp = "all positions have exactly 1 allele")
      } else {
        expect_no_error({
          a <- hla_alignments(my_gene, type = my_type, release = my_release, verbose = TRUE)
        }, message = "^(Downloading|Writing|Reading)")
        expect_named(a, c("sequences", "alleles", "onehot", "file", "release"))
        expect_equal(nrow(a$sequences), nrow(a$alleles))
        expect_equal(nrow(a$sequences), nrow(a$onehot))
      }
  })
}
