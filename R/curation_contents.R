############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-01-15 17:59:25
# 
# descr: curation content analysis
# 
############################################################

#' Signature stats
#'
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param col character. A column of \code{dat} that should be tabled. 
#' @return A sorted table with signature stats for the chosen metadata column
#' @export
sigStat <- function(dat, col) sort(table(dat[,col]), decreasing=TRUE)

#' Paper stats
#'
#' @param dat a \code{data.frame} storing BugSigDB data
#' @param col character. A column of \code{dat} that should be tabled. 
#' @return A sorted table with paper stats for the chosen metadata column
#' @export
paperStat <- function(dat, col)
{
    l <- split(dat[,col], dat[,"PMID"])
    ul <- vapply(l, function(x) x[1], dat[1,col])
    sort(table(ul), decreasing=TRUE)
}


