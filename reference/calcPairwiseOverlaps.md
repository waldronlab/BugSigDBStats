# Calculate pairwise overlaps between all signatures

Calculate pairwise overlaps between all signatures

## Usage

``` r
calcPairwiseOverlaps(sets, targetset = NULL)
```

## Arguments

- sets:

  a named list of signatures

- targetset:

  NULL to test all pairwise overlaps, or the name of one element of the
  \`sets\` list. If the name of one element of \`sets\` is provided,
  overlaps will be calculated only with that one signature.

## Value

a \`data.frame\` with one row per pairwise overlap calculated, and
colnames:

name1 name2 length_set1 length_set2 length_union length_intersection
overlap jaccard

## Details

The hard work in this code is is by John Blischak from the blog post on
efficient calculation of pairwise overlaps between list elements at
https://blog.jdblischak.com/posts/pairwise-overlaps/. See that page for
definitions of overlap and jaccard. This function adds a few columns, an
option to calculate pairwise overlaps with one list element only, and
documentation.

## Examples

``` r
testlist <- list(a = 1:3, b = 3, c = 3:4)
(all <- calcPairwiseOverlaps(testlist))
#>   name1 name2 length_set1 length_set2 length_union length_intersection overlap
#> 1     a     b           3           1            3                   1     1.0
#> 2     a     c           3           2            4                   1     0.5
#> 3     b     c           1           2            2                   1     1.0
#>     jaccard
#> 1 0.3333333
#> 2 0.2500000
#> 3 0.5000000
calcPairwiseOverlaps(testlist, targetset = "b")
#>   name1 name2 length_set1 length_set2 length_union length_intersection overlap
#> 1     b     a           1           3            3                   1       1
#> 2     b     c           1           2            2                   1       1
#>     jaccard
#> 1 0.3333333
#> 2 0.5000000
##
## Calculate overlaps between existing signatures with one additional signature
testlist <- c(testlist, d = list(4:5))
calcPairwiseOverlaps(testlist, targetset = "d")
#>   name1 name2 length_set1 length_set2 length_union length_intersection overlap
#> 1     d     c           2           2            3                   1     0.5
#>     jaccard
#> 1 0.3333333
```
