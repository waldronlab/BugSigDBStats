# Weighted semantic similarity

Incorporation of weights into the computation of semantic similarity

## Usage

``` r
weightedBMA(
  o,
  ic,
  set1,
  set2,
  weights1 = rep(1, length(set1)),
  weights2 = rep(1, length(set2))
)
```

## Arguments

- o:

  ontology. An object of class `ontology_index`.

- ic:

  information content. Typically obtained via
  `ontologySimilarity::descendants_IC(o)`.

- set1:

  character. First term set.

- set2:

  character. Secound term set.

- weights1:

  numeric. Weights for each term of the first term set, ie. must be
  parallel to `set1`.

- weights2:

  numeric. Weights for each term of the second term set, ie. must be
  parallel to `set2`.

## Value

a numeric value expressing semantic similarity between the two input
term sets, weighted by the given individual term weights.
