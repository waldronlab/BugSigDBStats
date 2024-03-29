% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcSemanticSimilarity.R
\name{testReplicability}
\alias{testReplicability}
\title{Test replicability of microbiome changes across studies}
\usage{
testReplicability(
  df,
  onto = NULL,
  multi.contrast = c("all", "first", "largest", "merge"),
  min.studies = 2,
  min.taxa = 5
)
}
\arguments{
\item{df}{\code{data.frame} storing BugSigDB data. Typically obtained via
\code{\link{importBugSigDB}}.}

\item{onto}{ontology. An object of class \code{ontology_index} storing the 
NCBI Taxonomy. Typically obtained via \code{\link{getNcbiTaxonomyObo}}.}

\item{multi.contrast}{character. How to treat studies that report multiple
contrasts for a study to avoid detection of duplication within studies as
replication? Select one of \itemize{
\item \code{"all"} incorporates all contrasts reported by a study. This is 
equivalent to do nothing, ie computes similarity between signatures taking all,
and thus potential similar / duplicated, signatures for each study into account.  
\item \code{"first"} take only the first contrast reported by a study, ie the
first signature with increased abundance and the first signature with decreased
abundance in the study group.
\item \code{"largest"} take only the largest signature with increased and
decreased abundance in the study group, respectively. 
\item \code{"merge"} merge signatures with the same direction of abundance
change (increased or decreased) within studies. 
}}

\item{min.studies}{integer. Minimum number of studies for a condition to be
tested. Defaults to 2, which will then only test a condition investigated by
at least two studies.}

\item{min.taxa}{integer. Minimum size for a signature to be included. Defaults to 5, which
will then only include signatures containing at least 5 taxa.}
}
\value{
A data.frame reporting semantic similarity and re-sampling p-value for
each condition under investigation. Results are stratified by direction of abundance
change in the study group.
}
\description{
Function that implements a statistical test for replicability of microbiome
changes for experimental conditions in BugSigDB
}
\examples{
 dat <- bugsigdbr::importBugSigDB(version = "10.5281/zenodo.5904281")
 dat.feces <- subset(dat, `Body site` == "feces")
 res <- testReplicability(dat.feces)
}
