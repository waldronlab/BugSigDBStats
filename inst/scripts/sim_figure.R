############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-12-20 16:33:05
# 
# descr: 
# 
############################################################

library(bugsigdbr)
library(BugSigDBStats)
library(ComplexHeatmap)
library(ontologyIndex)
library(ontologySimilarity)

# (0) Obtain the data
dat <- bugsigdbr::importBugSigDB()
onto <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/ncbitaxon.obo")

# (1) Semantic similarity on mixed taxonomic levels
dat.feces <- subset(dat, `Body site` == "feces")
ind <- lengths(dat.feces[["NCBI Taxonomy IDs"]]) > 4
dat.feces <- dat.feces[ind,]
sigs <- bugsigdbr::getSignatures(dat.feces, tax.id.type = "ncbi")
sigs <- lapply(sigs, function(s) paste0("NCBITaxon:", s))
utax <- unique(unlist(sigs))
nt <- utax[!(utax %in% onto$id)]
sigs <- lapply(sigs, function(s) setdiff(s, nt))
sim.mat <- ontologySimilarity::get_sim_grid(ontology = onto, term_sets = sigs)

# parse out condition from signature name
spl <- strsplit(rownames(sim.mat), "_")
spl <- vapply(spl, `[`, character(1), x = 2)
spl <- strsplit(spl, ":")
spl <- vapply(spl, `[`, character(1), x = 1)
head(spl)

# hierarchical clustering on the semantic similarity matrix
hc <- stats::hclust(stats::as.dist(1 - sim.mat), method = "ward.D")
clus <- stats::cutree(hc, 5)
head(clus)
d <- data.frame(label = names(clus), count = lengths(sigs[names(clus)]))
head(d)

# (2) Jaccard matrix on genus level
sigs.genus <- bugsigdbr::getSignatures(dat.feces,
                                       tax.id.type = "metaphlan",
                                       tax.level = "genus",
                                       exact.tax.level = FALSE)
mydists <- BugSigDBStats::calcPairwiseOverlaps(sigs.genus)
signames <- unique(c(mydists$name1, mydists$name2))
jmat <- matrix(NA, nrow=length(signames), ncol=length(signames), dimnames=list(signames, signames))
diag(jmat) <- 1
for (i in seq_len(nrow(mydists)))
{
    jmat[mydists[i, "name2"], mydists[i, "name1"]] <- mydists[i, "jaccard"]
    jmat[mydists[i, "name1"], mydists[i, "name2"]] <- mydists[i, "jaccard"]
}


# Correlation
.fcor <- function(i) cor(jmat[i,], sim.mat[rownames(jmat)[i], colnames(jmat)])
cors <- vapply(seq_len(nrow(jmat)), .fcor, numeric(1))
summary(cors)

# Heatmap

# color ramp from 0.01 quantile to 0.99 quantile
quantile(as.vector(sim.mat), 0.01)
quantile(as.vector(sim.mat), 0.99)
col2 <- circlize::colorRamp2(c(0,0.8012), c("#EEEEEE", "red"))

h2 <- ComplexHeatmap::Heatmap(sim.mat[rownames(jmat),rownames(jmat)], 
                              name = "semsim",
                              col = col2,
                              row_title = "signatures",
                              column_title = "signatures"
                              show_row_names=FALSE,
                              show_column_names=FALSE)

h1 <- Heatmap(jmat[,column_order(h2)], 
              col = col,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              name = "jaccard")

h2 + h1

# zoom into a cluster
h2 <- draw(h2)
ro <- ComplexHeatmap::row_order(h2)
rd <- ComplexHeatmap::row_dend(h2)
str(rd, max.level = 2)

# cut dendrogram at specified height
hd <- cut(rd, 2.90)

# plot dendrogram
par(mar = c(0,0,4,20))
par(cex = 0.5)
plot(hd$lower[[1]], horiz = TRUE)
