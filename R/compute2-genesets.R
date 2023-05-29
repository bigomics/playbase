##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Normalize geneset matrix by row
#'
#' @param g a sparse matrix (eg. genesets x genes)
#'
#' @return
#' @export
#'
#' @examples
normalize_matrix_by_row <- function(G){
        # efficient normalization using linear algebra
        row_sums <- Matrix::rowSums(G)
        D <- Matrix::Diagonal(x = 1/row_sums)
        G_scaled <- D %*% G
        rownames(G_scaled) <- rownames(G)
        return(G_scaled)
    }

#' Title
#'
#' @param pgx value
#' @param max.features value
#' @param test.methods value
#' @param remove.outputs value
#'
#' @return
#' @export
#'
#' @examples
compute_testGenesets <- function(pgx,
                                 max.features=1000,
                                 custom.geneset=NULL,
                                 test.methods = c("gsva","camera","fgsea"),
                                 remove.outputs=TRUE )
{
    ## Rewritten 24.12.2019. Now much faster but needs gset-sparseG-XL
    ## precomputed.
    ##
    ##
    if(!"X" %in% names(pgx)) {
        stop("[compute_testGenesets] FATAL : object must have normalized matrix X")
    }

    # Load custom genesets (if user provided)

    if(!is.null(custom.geneset$gmt)) {
        # convert gmt standard to SPARSE matrix
        custom_gmt <- playbase::createSparseGenesetMatrix(custom.geneset$gmt, min_gene_frequency=1)

    }

    ##-----------------------------------------------------------
    ## Load huge geneset matrix
    ##-----------------------------------------------------------

    G <- playdata::GSET_SPARSEG_XL

    ##-----------------------------------------------------------
    ## Filter genes
    ##-----------------------------------------------------------

    ## filter genes only in dataset

    ##
    ##GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    ##genes = Matrix::head(as.character(unlist(as.list(org.Hs.egSYMBOL))),1000)
    genes = unique(as.character(pgx$genes$gene_name))
    genes <- toupper(genes)  ## handle mouse genes...
    G <- G[,colnames(G) %in% genes]
    
    # Normalize G after removal of genes

    G <- playbase::normalize_matrix_by_row(G)

    if(!is.null(custom.geneset$gmt)) {
        
        custom_gmt <- custom_gmt[,colnames(custom_gmt) %in% genes]
        custom_gmt <- playbase::normalize_matrix_by_row(custom_gmt)
        
        # combine standard genesets with custom genesets
        G <- playbase::merge_sparse_matrix(G,custom_gmt)
    }

    # Transpose G

    G <- Matrix::t(G)

    ##-----------------------------------------------------------
    ## Filter gene sets
    ##-----------------------------------------------------------

    ## filter gene sets on size
    cat("Filtering gene sets on size...\n")
    gmt.size = Matrix::colSums(G!=0)
    size.ok <- (gmt.size >= 15 & gmt.size <= 400 )
    G <- G[, which(size.ok)]
    dim(G)
    table(sub(":.*","",colnames(G)))

    ##-----------------------------------------------------------
    ## create the full GENE matrix (always collapsed by gene)
    ##-----------------------------------------------------------

    single.omics <- !any(grepl("\\[",rownames(pgx$counts)))
    single.omics=TRUE  ## !!! for now...
    if(single.omics) {
        ## normalized matrix
        X <- pgx$X
    } else {
        data.type <- gsub("\\[|\\].*","",rownames(pgx$counts))
        jj <- which(data.type %in% c("gx","mrna"))
        if(length(jj)==0) {
            stop("FATAL. could not find gx/mrna values.")
        }
        X <- pgx$X[jj,]
    }

    ## if reduced samples
    ss <- rownames(pgx$model.parameters$exp.matrix)
    X <- X[,ss,drop=FALSE]
    dim(X)

    ##-----------------------------------------------------------
    ## create the GENESETxGENE matrix
    ##-----------------------------------------------------------
    cat("Matching gene set matrix...\n")
    gg <- toupper(rownames(X)) ## accomodate for mouse...
    ii <- intersect(gg,rownames(G))
    G <- G[ii,,drop=FALSE]
    xx <- setdiff(gg,rownames(G))
    matX <- Matrix::Matrix(0, nrow=length(xx), ncol=ncol(G), sparse=TRUE)
    rownames(matX) <- xx
    colnames(matX) <- colnames(G)
    G <- rbind(G, matX)
    G <- G[match(gg,rownames(G)),,drop=FALSE]
    rownames(G) <- rownames(X) ## original name (e.g. mouse)
    dim(G)
    dim(X)

    ##-----------------------------------------------------------
    ## Prioritize gene sets by fast rank-correlation
    ##-----------------------------------------------------------
    if(is.null(max.features)) max.features <- 20000
    if(max.features < 0) max.features <- 20000
    max.features

    if(max.features > 0) {
        cat("Reducing gene set matrix...\n")

        ## Reduce gene sets by selecting top varying genesets. We use the
        ## very fast sparse rank-correlation for approximate single sample
        ## geneset activation.
        cX <- X - rowMeans(X, na.rm=TRUE)  ## center!
        dim(cX)
        gsetX = qlcMatrix::corSparse( G[,], apply( cX[,],2,rank) )
        ## gsetX = limma::normalizeQuantiles(gsetX) ##???
        ##grp <- pgx$samples$group
        grp <- pgx$model.parameters$group
        gsetX.bygroup <- NULL
        if(!is.null(grp)) {
            gsetX.bygroup <- Matrix::t(apply(gsetX,1,function(x) tapply(x,grp,mean)))
            sdx <- apply(gsetX.bygroup,1,sd)
        } else {
            sdx <- apply(gsetX,1,sd)
        }
        names(sdx) <- colnames(G)
        jj = Matrix::head(order(-sdx), max.features)
        must.include <- "hallmark|kegg|^go|^celltype|^pathway|^custom"
        jj = unique( c(jj, grep(must.include,colnames(G),ignore.case=TRUE)))
        jj = jj[order(colnames(G)[jj])]
        G = G[,jj,drop=FALSE]
    }
    dim(G)

    ##-----------------------------------------------------------
    ## get design and contrast matrix
    ##-----------------------------------------------------------
    design = pgx$model.parameters$design
    contr.matrix = pgx$model.parameters$contr.matrix
    ##contr.matrix
    ##exp.matrix = (design %*% contr.matrix)

    ALL.GSET.METHODS=c("fisher","ssgsea","gsva", "spearman", "camera", "fry",
                       "fgsea","gsea.permPH","gsea.permGS","gseaPR")
    test.methods

    ##-----------------------------------------------------------
    ## Run methods
    ##-----------------------------------------------------------
    cat(">>> Testing gene sets with methods:",test.methods,"\n")
    test.methods

    ## convert to gene list
    gmt <- lapply(apply(G!=0, 2, which),names)
    Y <- pgx$samples
    gc()

    gset.meta = playbase::gset.fitContrastsWithAllMethods(
        gmt = gmt, X = X, Y = Y, G = G,
        design = design, ## genes=GENES,
        contr.matrix = contr.matrix, methods = test.methods,
        mc.threads = 1, mc.cores = NULL, batch.correct = TRUE
    )

    rownames(gset.meta$timings) <- paste("[test.genesets]",rownames(gset.meta$timings))
    print(gset.meta$timings)

    pgx$timings <- rbind(pgx$timings, gset.meta$timings)
    pgx$gset.meta <- gset.meta

    if(0) {
        ## average expression of geneset members
        n1 <- Matrix::rowSums(G!=0)
        n2 <- Matrix::rowMeans(abs(pgx$X) > 1)
        table( n1 > 20 & n2 > 0.05 )
        ##ii  <- which( n1 > 20 & n2 > 0.05 ) ## make faster...
        ii <- Matrix::head(order( -1*n1*n2 ),4000) ## make faster...
        G1 <- Matrix::t(G[ii,]!=0)
        X1 <- pgx$X[ii,]
        ng <- Matrix::colSums(G[ii,]!=0)
        meta.matrix <- as.matrix(G1 %*% X1) / ng
        dim(meta.matrix)
    }

    names(pgx$gset.meta$matrices)
    ##pgx$gsetX = pgx$gset.meta$matrices[["fc"]]  ## META or average FC??!
    pgx$gsetX = pgx$gset.meta$matrices[["meta"]]  ## META or average FC??!
    pgx$GMT <- G[,rownames(pgx$gsetX)]

    # calculate gset info and store as pgx$gset.meta

    gset.size <- Matrix::colSums(pgx$GMT != 0)

    gset.size.raw <- playdata::GSET_SIZE

    # combine standard genesets with custom genesets size vector

    if(!is.null(custom.geneset$gmt)) {
        gset.size.raw <- c(gset.size.raw, custom.geneset$info$GSET_SIZE)
    }

    gset.idx <- match(names(gset.size), names(gset.size.raw))

    gset.fraction <-  gset.size/gset.size.raw[gset.idx]

    pgx$gset.meta$info <- data.frame(
        gset.size.raw = gset.size.raw[gset.idx],
        gset.size = gset.size,
        gset.fraction = gset.fraction
    )
    
    ##-----------------------------------------------------------------------
    ##------------------------ clean up -------------------------------------
    ##-----------------------------------------------------------------------

    ## remove large outputs... (uncomment if needed)
    if(remove.outputs) {
        pgx$gset.meta$outputs <- NULL
    }

    remove(X)
    remove(Y)
    remove(G)
    remove(gmt)

    return(pgx)
}



#' Clean gmt file function
#'
#' @param gmt.all a list of gmt files
#' @param gmt.db a vector of gmt database names
#'
#' @return
#' @export
#'
#' @examples
clean_gmt <- function(gmt.all, gmt.db){
    gmt.db = toupper(gmt.db)
    for(i in 1:length(gmt.all)) {
    names(gmt.all[[i]]) <- sub("\\(GO:","(GO_",names(gmt.all[[i]]))
    names(gmt.all[[i]]) <- gsub("%","_",names(gmt.all[[i]])) # substitute % sign in wikipathways
    names(gmt.all[[i]]) <- sub(":","",names(gmt.all[[i]]))
    ## names(gmt.all[[i]]) <- tolower(names(gmt.all[[i]]))
    names(gmt.all[[i]]) <- paste0(toupper(gmt.db[i]),":",names(gmt.all[[i]]))
    }
    j0 = grep("_up", names(gmt.all))
    j1 = grep("_down", names(gmt.all))
    for(i in j0) {
        names(gmt.all[[i]]) <- paste0(names(gmt.all[[i]])," (up)")
    }
    for(i in j1) {
        names(gmt.all[[i]]) <- paste0(names(gmt.all[[i]])," (down)")
    }
    names(gmt.all) <- NULL
    gmt.all <- unlist(gmt.all,recursive=FALSE, use.names=TRUE)

    ## get rid of trailing numeric values
    gmt.all <-  parallel::mclapply(gmt.all, function(x) gsub("[,].*","",x), mc.cores=1)

    ## order by length and take out duplicated sets (only by name)
    gmt.all <- gmt.all[order(-sapply(gmt.all,length))]
    gmt.all <- gmt.all[!duplicated(names(gmt.all))]

    gmt.all <- gmt.all[order(names(gmt.all))]

    return(gmt.all)
}

#' Build sparse gene set matrix
#'
#' @param gmt.all a list of gmt files
#' @param gmt.db a vector of gmt database names
#'
#' @return
#' @export
#'
#' @examples
createSparseGenesetMatrix <- function(
    gmt.all,
    min.geneset.size = 15,
    max.geneset.size = 500,
    min_gene_frequency= 10) {
    
    ## ----------- Get all official gene symbols
    symbol <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
    known.symbols <- sort(unique(unlist(symbol)))
    
    ##------------- filter by size
    gmt.size <- sapply(gmt.all,length)

    gmt.all <- gmt.all[which(gmt.size >= 15 & gmt.size <= 1000)]
    
    ## ------------- filter genes by minimum frequency and chrom
    genes.table <- table(unlist(gmt.all))
    
    #genes <- names(which( genes.table >= 10 & genes.table <= 1000  ))
    genes <- names(which( genes.table >= min_gene_frequency ))
    genes <- genes[grep("^LOC|RIK$",genes,invert=TRUE)]
    genes <- intersect(genes, known.symbols)
    annot <- playbase::ngs.getGeneAnnotation(genes)
    genes <- genes[ !is.na(annot$chr) ]
    
    ## Filter genesets with permitted genes (official and min.sharing)
    gmt.all <- parallel::mclapply(gmt.all, function(s) intersect(s,genes))
    gmt.size <- sapply(gmt.all,length)
    # gmt.all <- gmt.all[which(gmt.size >= 15 & gmt.size <= 200)]
    gmt.all <- gmt.all[which(gmt.size >= min.geneset.size & gmt.size <= max.geneset.size)] #legacy
    ## build huge sparsematrix gene x genesets
    genes <- sort(genes)
    idx.j <- parallel::mclapply(gmt.all[], function(s) match(s,genes))
    idx.i <- lapply(1:length(gmt.all), function(i) rep(i,length(idx.j[[i]])))
    ii <- unlist(idx.i)
    jj <- unlist(idx.j)
    
    G <- Matrix::sparseMatrix( i=ii, j=jj, x=rep(1,length(ii)),
                      dims = c(length(gmt.all), length(genes)) )
    colnames(G) = genes
    rownames(G) = names(gmt.all)

    return(G)
}

#' Merge two sparse matrices
#'
#' @param m1 matrix 1 in Matrix format
#' @param m2 matrix 2 in Matrix format
#'
#' @return
#' @export
#'
#' @examples
merge_sparse_matrix <- function(m1, m2) {
    num_cols1 <- ncol(m1)
    num_cols2 <- ncol(m2)

    gene_vector <- unique(c(colnames(m1), colnames(m2)))
  
    if (num_cols1 < length(gene_vector)) {
        num_missing_cols <- length(gene_vector) - num_cols1
        zero_cols <- Matrix::Matrix(0, nrow = nrow(m1), ncol = num_missing_cols, sparse = TRUE)
        genes_missing_in_m1 <- setdiff(gene_vector, colnames(m1))
        colnames(zero_cols) <- genes_missing_in_m1

        m1 <- cbind(m1, zero_cols)

        m1 <- m1[, gene_vector]
    }
  
    if (num_cols2 < length(gene_vector)) {
        num_missing_cols <- length(gene_vector) - num_cols2
        zero_cols <- Matrix::Matrix(0, nrow = nrow(m2), ncol = num_missing_cols, sparse = TRUE)
        genes_missing_in_m2 <- setdiff(gene_vector, colnames(m2))
        colnames(zero_cols) <- genes_missing_in_m2
        m2 <- cbind(m2, zero_cols)
        m2 <- m2[, gene_vector]
    }
    
    combined_gmt <- Matrix::rbind2(m1, m2)

  return(combined_gmt)
}
