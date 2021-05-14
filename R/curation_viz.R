############################################################
# 
# author: Ludwig Geistlinger
# date: 2019-03-13 17:01:25
# 
# descr: visualization of curated signatures
# 
###########################################################

#' Plot curation output as a function of time
#'
#' @param dat a \code{data.frame} storing BugSigDB data.
#' @param col character. A column of \code{dat} that contain the curation date
#' in dmy format. 
#' @param diff logical. Display only the difference between months? Defaults to 
#' \code{FALSE} which will then display total cumulative numbers.
#' @return None. Plots to a graphics device.
#' @importFrom graphics barplot par text
#' @export
plotProgressOverTime <- function(dat, col = "Curated date", diff = FALSE)
{
    dates <- dat[,col]
    # fill empty dates
    for(i in seq_along(dates)) if(dates[i] == "") dates[i] <- dates[i-1]
    # signatures
    dbm <- substring(dates, 1, 7)
    dat[,col] <- dbm
    dbm <- sort(dbm)
    cdbm <- cumsum(table(dbm))
    # papers
    ind <- !duplicated(dat[,"PMID"])
    pbm <- split(dat[ind,"PMID"], dat[ind,col])
    pbm <- lengths(lapply(pbm, unique))
    cpbm <- cumsum(pbm)
    # plot
    dbm <- rbind(cdbm, cpbm)
    main <- ""
    if(diff)
    { 
        dbm <- dbm[,2:ncol(dbm)] - dbm[,1:(ncol(dbm)-1)]
        means <- round(rowMeans(dbm[,2:(ncol(dbm) - 1)]))
        main <- paste("Mean (papers / sigs):", means[2], "/", means[1])
    }

    par(las = 2)
    par(las = 1)
    bp <- barplot(dbm, beside=TRUE, horiz=TRUE, xlim=c(0, max(dbm[1,]) + ifelse(diff, 20, 100)),
            legend.text=c("signatures", "papers"), args.legend=list(x="bottomright"),
            main = main)
    for(i in 1:ncol(dbm)) text(y=bp[,i], x=dbm[,i], labels=dbm[,i], pos=4) 
}

#' Plot curation output for each curator
#'
#' @param dat a \code{data.frame} storing BugSigDB data.
#' @param npc character. A column of \code{dat} that contain the curation date
#' in dmy format. 
#' @return None. Plots to a graphics device.
#' @export
plotCuratorStats <- function(dat, npc)
{
    par(las=1)
    par(mar=c(5, 5, 4, 1))
    title <- paste("#papers:", length(unique(dat[,"PMID"])), 
                    ", #signatures:", nrow(dat))
    bp <- barplot(npc, xlim=c(0, max(npc) + 25), beside=TRUE, horiz=TRUE, 
                    main=title, legend.text=c("signatures", "papers"), 
                    args.legend=list(x="bottomright"))
    par(cex=0.8)
    for(i in 1:ncol(npc)) text(y=bp[,i], x=npc[,i], labels=npc[,i], pos=4) 
}


