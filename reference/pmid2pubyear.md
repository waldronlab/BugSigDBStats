# Get the publication year for a PMID

Get the publication year for a PMID

## Usage

``` r
pmid2pubyear(pmid)
```

## Arguments

- pmid:

  character. A PubMed ID.

## Value

A character containg the publication year for the given `pmid`.

## Examples

``` r
 pmid2pubyear("32026945")
#> Warning: URL 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?tool=bioconductor&rettype=xml&retmode=text&db=PubMed&id=32026945': status was 'Couldn't connect to server'
#> Error in file(file, "r"): cannot open the connection to 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?tool=bioconductor&rettype=xml&retmode=text&db=PubMed&id=32026945'
```
