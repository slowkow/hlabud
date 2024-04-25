
test_that("hla_divergence fails with a non-existant allele name", {
  my_alleles <- c("DRB1*01:01,DRB1*02:01")
  my_release <- "3.52.0"
  expect_error({
    hla_divergence(my_alleles, release = my_release)
  }, regexp = "Input 'DRB1\\*02:01' not found")
})

test_that("hla_divergence works with a valid allele names", {
  my_release <- "3.52.0"
  a <- hla_alignments("DRB1", release = my_release)
  set.seed(42)
  my_alleles <- replicate(n = 100, expr = paste(sample(names(a$sequences), size = 2), collapse = ","))
  # my_alleles <- c("DRB1*01:01,DRB1*03:01")
  my_div <- hla_divergence(my_alleles, release = my_release)
  expect_length(my_div, length(my_alleles))
  expect_vector(my_div, ptype = numeric(), size = length(my_alleles))
  expect_true(all(my_div >= 0))
})
