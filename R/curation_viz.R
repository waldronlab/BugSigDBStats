############################################################
# 
# author: Ludwig Geistlinger
# date: 2019-03-13 17:01:25
# 
# descr: visualization of curated signatures
# 
###########################################################

# network similarity plot
#' @param smat numeric similarity matrix. Must be symmetric.
#' @param sigs signatures. Named list of character vectors. 
#' @param min.sim numeric. Minimum pairwise similarity in [0,1] for an edge to be included.
#' @param color numeric. Named numeric vector parallel to sigs. 
library(igraph)
library(ggraph)
networkSimPlot <- function(smat, sigs, min.sim = 0.2, color, lwd = 1)
{
    stopifnot(all(rownames(smat) == colnames(smat)))
    stopifnot(all(rownames(smat) == names(sigs)))
    rownames(smat) <- colnames(smat) <- names(sigs) <- vapply(rownames(smat), .getShortName, character(1))

    if (!is.numeric(min.sim) | min.sim < 0 | min.sim > 1) {
        stop("\"min.sim\" should be a number between 0 and 1.")
    }

    wd <- reshape2::melt(smat)
    wd <- wd[wd[, 1] != wd[, 2], ]
    wd <- wd[!is.na(wd[, 3]), ] 
    g <- igraph::graph.data.frame(wd[, -3], directed = FALSE)

    igraph::E(g)$width <- sqrt(wd[, 3] * 5) * lwd
    igraph::E(g)$weight <- wd[, 3]
    g <- igraph::delete.edges(g, igraph::E(g)[wd[, 3] < min.sim])
    idx <- unlist(sapply(igraph::V(g)$name, function(x) which(x == names(sigs))))
    cnt <- lengths(sigs[idx])
    
    igraph::V(g)$size <- cnt
    colVar <- color[idx] 
    igraph::V(g)$color <- colVar

    p <- ggraph(g, layout = "nicely")
    p <- p + geom_edge_link(alpha = 0.8, aes_(width = ~I(width)), 
            colour = "darkgrey")

    # p <- add_category_nodes(p = p, cex_category = cex_category, color = color)
    p <- p + ggnewscale::new_scale_fill() + geom_point(shape = 21, 
        aes_(x = ~x, y = ~y, fill = ~color, size = ~size)) + 
        scale_size_continuous(name = "number of genes", range = c(3, 
            8) * cex_category) + scale_fill_continuous(low = "red", 
        high = "blue", name = color, guide = guide_colorbar(reverse = TRUE)) + 
        theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + 
        theme(panel.background = element_blank())

    p + coord_equal() + guides(size = guide_legend(order = 1), 
        fill = guide_colorbar(order = 2))
    return(g)
}

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

#' Plot curation output for each curator
#'
#' @param dat a \code{data.frame} storing BugSigDB data.
#' @param date.col character. A column of \code{dat} that contain the curation date
#' in dmy format. 
#' @return None. Plots to a graphics device.
#' @export
plotUniqueMicrobesOverTime <- function(dat,
                                       date.col = "Curated date")
{
    dat <- dat[dat[,date.col] != "",]
    msc <- dat[["MetaPhlAn taxon names"]]
    dates <- dat[,date.col]
    dbm <- substring(dates, 1, 7)
    msc.spl <- split(msc, dbm)
    
    msc.spl <- lapply(msc.spl, function(x) unique(unname(unlist(x))))
    for(i in 2:length(msc.spl)) msc.spl[[i]] <- union(msc.spl[[i - 1]], msc.spl[[i]]) 
    nums <- lengths(msc.spl)
    df <- data.frame(date = names(nums), nr.microbes = nums)

    ggpubr::ggscatter(df[-nrow(df),], x = "date", y = "nr.microbes", 
                      ylab = "Number of unique microbes", xlab = "", 
                      ggtheme = ggplot2::theme_bw(), color = "darkblue") + 
                      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


.getNrCommonSignatures <- function(m1, m2, cmat, antagonistic)
{
    by <- ifelse(antagonistic, 2, 1)
    grid1 <- seq(1, ncol(cmat), by = by)
    grid2 <- seq(by, ncol(cmat), by = by)
    sum(cmat[m1,grid1] == 1 & cmat[m2,grid2] == 1)
}

.cooc <- function(umicrobes, cmat, antagonistic)
{
    len <- length(umicrobes)
    grid <- seq_len(len)
    cooc.mat <- matrix(0, nrow = len, ncol = len)    
    
    for(i in grid)
    {
        for(j in grid)
            cooc.mat[i, j] <- .getNrCommonSignatures(umicrobes[i], 
                                                     umicrobes[j], 
                                                     cmat,
                                                     antagonistic)
    }
    
    return(cooc.mat)
}    

#' Plot microbe co-occurrence in a heatmap 
#'
#' @param dat a \code{data.frame} storing BugSigDB data.
#' @param sig.type character. Signature type. Use either \code{"increased"} or
#' \code{"decreased"} to subset to signatures with either increased or decreased abundance
#' in the exposed group, respectively. Default is \code{"both"} which will not 
#' subset by the direction of abundance change.
#' @param tax.level character. Either \code{"mixed"} or any subset of
#' \code{c("kingdom", "phylum", "class", "order", "family", "genus", "species",
#' "strain")}. This full vector is equivalent to \code{"mixed"}.
#' @param exact.tax.level logical. Should only the exact taxonomic level
#' specified by \code{tax.level} be returned? Defaults to \code{TRUE}.
#' If \code{FALSE}, a more general \code{tax.level} is extracted for
#' microbes given at a more specific taxonomic level.
#' @param antagonistic logical. Antagonistic co-occurrence, ie occurring together
#' but one microbe up the other one is down? Defaults to \code{FALSE}.
#' @param anno character. Taxonomic level that should be displayed as an annotation
#' bar. Defaults to \code{"phylum"}.
#' @param top integer. The number of microbes to display. Defaults to \code{100}.
#' @param fontsize integer. Fontsize for text. Defaults to \code{6}.
#' @param ... additional arguments to \code{ComplexHeatmap::Heatmap}.
#' @return Plots to a graphics device. Returns the co-oocurence matrix invisibly.
#' @export
microbeHeatmap <- function(dat, 
                           sig.type = c("both", "increased", "decreased"),
                           tax.level = "mixed",
                           exact.tax.level = FALSE,
                           antagonistic = FALSE,
                           anno = "phylum",
                           top = 100,
                           fontsize = 6,
                           ...)
{
    if(!requireNamespace("safe"))
        stop("Please install the 'safe' package to use 'microbeHeatmap'")

    # restrict by sig type, to paired UP/DOWN signatures, and by tax level
    sig.type <- match.arg(sig.type)
    if(sig.type %in% c("increased", "decreased")) 
        dat <- subset(dat, `Abundance in Group 1` == sig.type)
    dat <- bugsigdbr::restrictTaxLevel(dat, tax.level, exact.tax.level, min.size = 2)
    sigs <- bugsigdbr::getSignatures(dat, tax.id.type = "metaphlan")    

    # get connectivity matrix (signature <-> microbes)
    sink(tempfile())
    cmat <- safe::getCmatrix(sigs, as.matrix = TRUE, 
                             min.size = 0, prune = FALSE)
    sink()
    stopifnot(all(names(sigs) == colnames(cmat)))
    
    # restrict to most frequently co-occurring microbes 
    rs <- rowSums(cmat)
    if(top > length(rs)) top <- length(rs) 
    top <- sort(rs, decreasing = TRUE)[top]
    cmat <- cmat[rs > top,] 
    ind <- colSums(cmat) > 1
    cmat <- cmat[,ind]
    dat <- dat[ind,]    
    
    # co-occurrence also by occurring together inversely 
    # (ie antagonistic = one up, one down)
    if(antagonistic){ 
        stex <- paste(dat$Study, dat$Experiment)
        tab <- table(stex)
        paired <- names(tab)[tab == 2] 
        ind <- stex %in% paired 
        dat <- dat[ind,]
        cmat <- cmat[,ind]
    } 
    
    # calculate co-occurrence matrix
    umicrobes <- rownames(cmat)
    cooc.mat <- .cooc(umicrobes, cmat, antagonistic)
    
    # add annotation 
    uanno <- bugsigdbr::extractTaxLevel(umicrobes,
                                        tax.id.type = "taxname",
                                        tax.level = anno,
                                        exact.tax.level = FALSE)
    anno <- ComplexHeatmap::HeatmapAnnotation(phylum = uanno)
    
    # plot the heatmap
    n <- bugsigdbr::extractTaxLevel(umicrobes,
                                    tax.id.type = "taxname",
                                    tax.level = tax.level)
    rownames(cooc.mat) <- colnames(cooc.mat) <- n 
    print(ComplexHeatmap::Heatmap(log10(cooc.mat + 0.1), 
                            name = "log10 Co-occurence",
                            top_annotation = anno,
                            row_names_gp = gpar(fontsize = fontsize),
                            column_names_gp = gpar(fontsize = fontsize), 
                            ...))
    return(invisible(cooc.mat))
}
