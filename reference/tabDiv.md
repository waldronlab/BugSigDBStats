# Table diversity column

Table diversity column

## Usage

``` r
tabDiv(exps, div.col, spl.col, min.exps = 5, perc = FALSE)
```

## Arguments

- exps:

  a `data.frame` storing BugSigDB experiment data

- div.col:

  character. A column of `exps` that should be tabled.

- spl.col:

  character. A column of `exps` that should be used to split `div.col`
  into groups.

- min.exps:

  integer. Minimum number of experiments for a category in `spl.col` to
  be included in the result. Defaults to `5`.

- perc:

  logical. Return absolute frequencies or relative frequencies? Defaults
  to `FALSE` which will return absolute frequencies.

## Value

A sorted table with table stats for the chosen diversity column
