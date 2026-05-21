# Table column

Table column

## Usage

``` r
tabCol(df, col, n = Inf, perc = FALSE)
```

## Arguments

- df:

  a `data.frame` storing BugSigDB data

- col:

  character. A column of `dat` that should be tabled.

- n:

  integer. Number of categories to show. Defaults to `Inf` which will
  then show all categories in the chosen column.

- perc:

  logical. Return absolute frequencies or relative frequencies? Defaults
  to `FALSE` which will return absolute frequencies.

## Value

A sorted table with table stats for the chosen metadata column
