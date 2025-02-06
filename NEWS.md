
hlabud 2.1.0
============

## Bug fixes

* Fix incorrect position numbering again, thanks to Vinicius Stelet for
  bringing this to my attention in issue #4.

In hlabud version 2.0.0, some genes were correctly numbered and some were not:

    A     incorrect
    B     incorrect
    C     ok
    DMA   incorrect
    DMB   incorrect
    DOA   incorrect
    DOB   incorrect
    DPA1  ok
    DPB1  ok
    DQA1  ok
    DQA2  incorrect
    DQB1  ok
    DQB2  ok
    DRA   incorrect
    DRB   ok
    DRB1  ok
    DRB3  ok
    DRB4  ok
    DRB5  ok
    E     ok
    F     incorrect
    G     incorrect
    HFE   incorrect
    MICA  incorrect
    MICB  incorrect
    TAP1  incorrect
    TAP2  ok

We now have additional tests to confirm that the IMGT files are being parsed
correctly, and the positions are numbered correctly.


hlabud 2.0.0
============

## Bug fixes

* Fix incorrect position numbering, accounting for insertions and deletions that are indicated with the "." character. Thanks to Vinicius Stelet for bringing this to my attention in issue #3.

* Instead of discarding positions with `*`, we include them and label them as `unk`, for example `pos241_unk` indicates an unknown amino acid at position 241. Thanks to Sreekar Mantena for reporting this issue!

* Fix an off-by-one error. For example, HLA-A had `pos361_-` in the `colnames(a$onehot)` but it should have had the reference allele instead of `-`. This is now fixed. Thanks again to Sreekar Mantena for reporting this issue!

## Changes

* Change position names from `pos21_D` to `D21`. When negative, `posn21_D` to `Dn21`.

* Change `dosage()` to take a one-hot matrix as the first argument.

* Change `dosage()` to return full allele names from IMGT when matching to partial allele names like `DRB1*03` or `DRB1*03:01`. And show messages indicating which alleles were matched when `verbose=TRUE`.

* Automatically overwrite `{hlabud_dir}/alleles.json` if it is older than 24 hours.

* Automatically overwrite `{hlabud_dir}/tags.json` if it is older than 24 hours.


hlabud 1.0.0
============

* Initial release.

* Added a `NEWS.md` file to track changes to the package.
