# Title Create a Jaccard or overlap dissimilarity matrix from the output of calcPairwiseOverlaps

Title Create a Jaccard or overlap dissimilarity matrix from the output
of calcPairwiseOverlaps

## Usage

``` r
makeDist(overlaps, use = "jaccard", names)
```

## Arguments

- overlaps:

  A data.frame output from the calcPairwiseOverlaps function

- use:

  Either "overlap" or "jaccard"

- names:

  A list of signature names

## Value

a Dist object of (1 - overlap) or (1 - jaccard)

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
sigs <- names(testlist)
makeDist(all, "jaccard", sigs)
#>           a         b
#> b 0.6666667          
#> c 0.7500000 0.5000000
makeDist(all, "overlap", sigs)
#>     a   b
#> b 0.0    
#> c 0.5 0.0
```
