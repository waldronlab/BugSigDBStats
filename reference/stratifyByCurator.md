# Stratify curation output by curator

Counts the number of signatures and papers each person contributed to.
Columns such as `"Revision editor"` can list several people per
signature (e.g. `"Claregrieve1,Chloe,WikiWorks743"`); these entries are
split so that each person is credited individually, regardless of their
position in the list.

## Usage

``` r
stratifyByCurator(
  dat,
  col = "Curator",
  sep = ",",
  exclude = "^WikiWorks",
  short.names = TRUE
)
```

## Arguments

- dat:

  a `data.frame` storing BugSigDB data

- col:

  character. The column of `dat` listing the people to count. Defaults
  to `"Curator"`; use `"Revision editor"` to count revision
  contributions.

- sep:

  character. Separator between multiple people in `col`.

- exclude:

  character. A regular expression; people whose name matches it are
  dropped. Defaults to dropping the `WikiWorks` maintenance accounts.
  Use `NULL` to keep everyone.

- short.names:

  logical. Shorten names to the first word for display (e.g.
  `"Aleru Divine"` becomes `"Aleru"`)? Names whose first word is shared
  by more than one person are kept in full.

## Value

a numeric `matrix` with rows `"signatures"` and `"papers"` and one
column per person, ordered by increasing number of papers. If no people
remain after removing empty and excluded entries, a warning is issued
and a matrix with zero columns is returned.

## Examples

``` r
dat <- data.frame(PMID = c(1, 1, 2, 3),
                  `Revision editor` = c("Ann,WikiWorks743", "Ann,Bob",
                                        "Bob", "WikiWorks"),
                  check.names = FALSE)
stratifyByCurator(dat, col = "Revision editor")
#>            Ann Bob
#> signatures   2   2
#> papers       1   2
```
