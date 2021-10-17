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

#' Table column
#'
#' @param df a \code{data.frame} storing BugSigDB data
#' @param col character. A column of \code{dat} that should be tabled. 
#' @param n integer. Number of categories to show. Defaults to \code{Inf}
#' which will then show all categories in the chosen column.
#' @param perc logical. Return absolute frequencies or relative frequencies?
#' Defaults to \code{FALSE} which will return absolute frequencies. 
#' @return A sorted table with table stats for the chosen metadata column
#' @importFrom utils head
#' @export
tabCol <- function(df, col, n = Inf, perc = FALSE)
{
    tab <- table(df[[col]])
    tab <- sort(tab, decreasing = TRUE)
    n <- ifelse(n < Inf, min(n, length(tab)), length(tab)) 
    htab <- head(tab, n)
    if(perc) htab <- signif(htab / sum(tab), digits = 3)
    return(htab)
}

#' Table diversity column
#'
#' @param exps a \code{data.frame} storing BugSigDB experiment data
#' @param div.col character. A column of \code{exps} that should be tabled. 
#' @param spl.col character. A column of \code{exps} that should be used to
#' split \code{div.col} into groups.
#' @param min.exps integer. Minimum number of experiments for a category in
#' \code{spl.col} to be included in the result. Defaults to \code{5}.
#' @param perc logical. Return absolute frequencies or relative frequencies?
#' Defaults to \code{FALSE} which will return absolute frequencies. 
#' @return A sorted table with table stats for the chosen diversity column
#' @export
tabDiv <- function(exps, div.col, spl.col, min.exps = 5, perc = FALSE)
{
    spl <- split(exps[[div.col]], exps[[spl.col]])
    spl <- lapply(spl, function(x) x[!is.na(x)])
    spl <- spl[lengths(spl) >= min.exps]
    tab <- lapply(spl, table)
    tab <- vapply(tab, function(x) x[c("increased", "decreased", "unchanged")], integer(3))
    tab[is.na(tab)] <- 0
    rownames(tab) <- c("increased", "decreased", "unchanged")
    if(perc) tab <- apply(tab, 2, function(x) signif(x / sum(x), digits = 2))
    return(t(tab))
}
