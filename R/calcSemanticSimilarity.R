#' Calculate semantic similarity
#'
#' \code{calcSemanticSimilarity} calculates the pairwise semantic similarity
#' between two sets of signature files exported from BugSigDB.org.
#'
#' @param file1 Signature exported from BugSigDB.org.
#' @param file2 Signature exported from BugSigDB.org.
#' @export
#'
#' @examples 
#'
#' \dontrun{
#' file1 <- "path/to/exported/file1.csv"
#' file2 <- "path/to/exported/file2.csv"
#' output <- calcSemanticSimilarity(file1, file2)
#' head(output)
#' }
#' 
calcSemanticSimilarity <- function(file1, file2) {
    
    on.exit(gc())
    
    sig1 <- .bugSigDBExportToSig(file1)
    sig2 <- .bugSigDBExportToSig(file2)
    ncbi_onto <- .getNcbiTaxonomyObo() 
    
    utax <- unique(unlist(sig1))
    nt <- utax[!(utax %in% ncbi_onto$id)]
    sig1 <- lapply(sig1, function(s) setdiff(s, nt))
    
    utax <- unique(unlist(sig2))
    nt <- utax[!(utax %in% ncbi_onto$id)]
    sig2 <- lapply(sig2, function(s) setdiff(s, nt))
    
    mat <- ontologySimilarity::get_sim_grid(
        ontology = ncbi_onto, term_sets = sig1, term_sets2 = sig2
    )
    
    output <- as.data.frame(as.table(mat))
    names(output) <- c("sig1", "sig2", "semantic_similarity")
    output
}

#' Get NCBI Taxonomy in OBO format
#'
#' \code{.getNcbiTaxonomyObo} gets the NCBI taxonomy in OBO format.
#'
#' @return NCBI taxonomy in OBO format
#' @importFrom utils download.file
#' @export
#' 
.getNcbiTaxonomyObo <- function() {

    onto <- .getResourceFromCache("ncbi.onto")

    if (is.null(onto)) {

        url <- paste0("https://sandbox.zenodo.org/record/954943/files/",
                      "ncbitaxon.rds?download=1")
        temp_file <- tempfile()
        message("Getting NCBI taxonomy ontology.")
        utils::download.file(url = url, destfile = temp_file)
        onto <- readRDS(temp_file)
        .cacheResource(onto, "ncbi.onto")
        message("NCBI taxonomy ontology has been saved in cache.")
        return(onto)

    } else {
        
        message("Retrieveing NCBI taxonomy ontology from cache.")
        return(onto)

    }
}

.bugSigDBExportToSig <- function(fname) {
    df <- utils::read.csv(
        fname, header = TRUE, comment.char = "#", check.names = FALSE
    )
    
    id_cols <- c("Signature page name", "Experiment", "Study")
    ids <- lapply(df[, id_cols], function(x) sub("^.+ ", "", x))
    sig_names <- vector("character", length(ids[[1]]))
    for (i in seq_along(sig_names)) {
        sig_names[i] <- 
            paste0("bsdb:", ids[[1]][i], "/", ids[[2]][i], "/", ids[[3]][i])
    }
    sig <- lapply(df[["NCBI Taxonomy IDs"]], function(x) {
        sig <- strsplit(x, ";")[[1]] |>
            vapply(\(y) sub("^.+\\|", "", y), character(1))
        names(sig) <- NULL
        paste0("NCBITaxon:", sig)
    })
    names(sig) <- sig_names
    sig
}





