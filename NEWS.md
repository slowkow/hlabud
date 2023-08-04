# hlabud 1.0.0.9999

* Instead of discarding positions with `*`, we include them and label them as `unk`, for example `pos241_unk` indicates an unknown amino acid at position 241. Thanks to Sreekar Mantena for reporting this issue!

* Fix an off-by-one error. For example, HLA-A had `pos361_-` in the `colnames(a$onehot)` but it should have had the reference allele instead of `-`. This is now fixed. Thanks again to Sreekar Mantena for reporting this issue!

# hlabud 1.0.0

* Added a `NEWS.md` file to track changes to the package.
