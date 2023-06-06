##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

##================================================================================
##========================= CONNECTIVITY FUNCTIONS ===============================
##================================================================================


#' @export
pgx.computeConnectivityScores <- function(pgx, sigdb, ntop=1000, contrasts=NULL,
                                          remove.le=FALSE, inmemory=FALSE )
{

    meta = pgx.getMetaFoldChangeMatrix(pgx, what="meta")
    colnames(meta$fc)

    is.h5ref <- grepl("h5$",sigdb)
    if(!is.h5ref) {
        cat("[pgx.computeConnectivityScores] ERROR: must be H5 formatted file\n")
        return(NULL)
    }
    if(inmemory) {
        cat("[pgx.computeConnectivityScores] *** using in-memory ***\n")
    }

    h5.file <- NULL
    if(file.exists(sigdb)) {
        h5.file <- sigdb
    }
    if(is.null(h5.file)) {
        cat("[pgx.computeConnectivityScores] ERROR: could not open H5 file\n")
        return(NULL)
    }

    if(is.null(contrasts)) {
        contrasts <- colnames(meta$fc)
    }
    contrasts <- intersect(contrasts, colnames(meta$fc))

    if(inmemory) {
        F = meta$fc[,contrasts]
        scores <- pgx.correlateSignatureH5.inmemory(
            meta$fc, h5.file = h5.file,
            nsig=100, ntop=ntop, nperm=9999)
    } else {
        scores <- list()
        ct <- contrasts[1]
        for(ct in contrasts) {
            fc <- meta$fc[,ct]
            names(fc) <- rownames(meta$fc)
            names(fc) <- toupper(names(fc)) ## for MOUSE!!
            res <- pgx.correlateSignatureH5(
                fc, h5.file = h5.file,
                nsig=100, ntop=ntop, nperm=9999)
            dim(res)
            scores[[ct]] <- res
        }
    }
    if(is.null(names(scores))) names(scores) <- contrasts

    ## remove leadingEdge (take too much memory!!!)
    if(remove.le) {
        for(j in 1:length(scores)) scores[[j]]$leadingEdge <- NULL
    }

  names(scores)
  return(scores)
}

## ntop=1000;nsig=100;nperm=10000
#' @export
pgx.correlateSignatureH5.inmemory <- function(F, h5.file, nsig=100, ntop=1000, nperm=1000)
{
    ##
    ##
    ##
    ##

    if(NCOL(F)==1 && class(F)=="numeric") {
        rn <- names(F)
        F <- matrix(F,ncol=1)
        rownames(F) <- rn
    }

    if(is.null(rownames(F))) stop("F must have rownames")
    ## mouse... mouse...
    rownames(F) <- toupper(rownames(F))

    ## or instead compute correlation on top100 fc genes (read from file)
    rn <- rhdf5::h5read(h5.file,"data/rownames")
    cn <- rhdf5::h5read(h5.file,"data/colnames")

    ## Entire matrix in memory????
    matG <- rhdf5::h5read(h5.file, "data/matrix")  ### whole matrix!!!!
    matG[which(matG < -999999)] <- NA
    rownames(matG) <- rn
    colnames(matG) <- cn

    mem1 <- round(object.size(matG)/1e9,2)
    cat("[pgx.correlateSignatureH5] object.size(matG)=",mem1,"Gb\n")  ## gigabytes....

    ## ---------------------------------------------------------------
    ## Compute simple correlation between query profile and signatures
    ## ---------------------------------------------------------------
    res <- list()
    i=1
    for(i in 1:ncol(F)) {

        gg <- intersect(rownames(F),rn)
        fc1 <- sort(F[gg,i])
        gg <- unique(names(c(Matrix::head(fc1,nsig), Matrix::tail(fc1,nsig))))
        ## gg <- intersect(gg,rn)
        remove(fc1)
        row.idx <- match(gg,rn)
        length(row.idx)

        ##G <- rhdf5::h5read(h5.file, "data/matrix", index=list(row.idx,1:length(cn)))  ### SLOW!!!
        rG <- matG[row.idx,,drop=FALSE]
        rG <- apply( rG,2,rank, na.last="keep" )
        dim(rG)
        dimnames(rG) <- list(rn[row.idx],cn)

        ## this FC signature
        fc <- F[,i]
        ##rG  <- apply( G[gg,], 2, rank, na.last="keep" )
        rfc <- rank( fc[gg], na.last="keep" ) ## rank correlation??
        ##rho <- stats::cor(rG, rfc, use="pairwise")[,1]
        rG[is.na(rG)] <- 0  ## NEED RETHINK: are missing values to be treated as zero???
        rfc[is.na(rfc)] <- 0
        rho <- stats::cor(rG, rfc, use="pairwise")[,1]

        remove(rG,rfc)

        ##--------------------------------------------------
        ## test tops signatures using fGSEA
        ##--------------------------------------------------

        sel <- Matrix::head(names(sort(-abs(rho))), ntop)
        sel.idx <- match(sel, cn)
        sig100.up <- rhdf5::h5read(h5.file, "signature/sig100.up",
                            index = list(NULL, sel.idx) )
        sig100.dn <- rhdf5::h5read(h5.file, "signature/sig100.dn",
                            index = list(NULL, sel.idx) )
        ##head(sig100.up,2)

        ## combine up/down into one (unsigned GSEA test)
        gmt <- rbind(sig100.up, sig100.dn)
        gmt <- unlist(apply(gmt, 2, list),recursive=FALSE)
        names(gmt) <- cn[sel.idx]
      
        ## use entire fc vector
        system.time(res1 <- fgsea::fgseaSimple(gmt, abs(fc), nperm=nperm))  ## really unsigned???
        dim(res1)

        ## ---------------------------------------------------------------
        ## Combine correlation+GSEA by combined score (NES*rho)
        ## ---------------------------------------------------------------
        jj <- match( res1$pathway, names(rho))
        res1$rho  <- rho[jj]
        res1$R2 <- rho[jj]**2
        res1$score <- (res1$R2 * res1$NES)

        fn <- colnames(F)[i]
        res[[fn]] <- res1[order(res1$score, decreasing=TRUE),,drop=FALSE]

        gc()

    }
    remove(matG)
    gc()

    return(res)
}

##ntop=1000;nsig=100;nperm=10000
#' @export
pgx.correlateSignatureH5 <- function(fc, h5.file, nsig=100, ntop=1000, nperm=10000,
                                     h5.data = "data/matrix", h5.rn="data/rownames",
                                     h5.cn="data/colnames")
{
    ##
    ##
    ##
    ##


    if(is.null(names(fc))) stop("fc must have names")
    ## mouse... mouse...
    names(fc) <- toupper(names(fc))

    ## or instead compute correlation on top100 fc genes (read from file)
    rhdf5::h5closeAll()
    rn <- rhdf5::h5read(h5.file,"data/rownames")
    cn <- rhdf5::h5read(h5.file,"data/colnames")

    ## ---------------------------------------------------------------
    ## Compute simple correlation between query profile and signatures
    ## ---------------------------------------------------------------
    gg <- intersect(names(fc),rn)
    fc1 <- sort(fc[gg])
    gg <- unique(names(c(Matrix::head(fc1,nsig), Matrix::tail(fc1,nsig))))
    ## gg <- intersect(gg,rn)
    length(gg)
    row.idx <- match(gg,rn)
    rhdf5::h5closeAll()
    G <- rhdf5::h5read(h5.file, "data/matrix", index=list(row.idx,1:length(cn)))
    dim(G)
    ##head(G[,1])
    G[which(G < -999999)] <- NA
    ##G[is.na(G)] <- 0  ## NEED RETHINK: are missing values to be treated as zero???
    dim(G)
    dimnames(G) <- list(rn[row.idx],cn)

    ## rank correlation??
    rG  <- apply( G[gg,], 2, rank, na.last="keep" )
    rfc <- rank( fc[gg], na.last="keep" )
    ##rho <- stats::cor(rG, rfc, use="pairwise")[,1]
    rG[is.na(rG)] <- 0  ## NEED RETHINK: are missing values to be treated as zero???
    rfc[is.na(rfc)] <- 0
    suppressWarnings( rho <- stats::cor(rG, rfc, use="pairwise")[,1] )

    remove(G,rG,rfc)

    ## --------------------------------------------------
    ## test tops signatures using fGSEA
    ## --------------------------------------------------

    sel <- Matrix::head(names(sort(-abs(rho))), ntop)
    sel.idx <- match(sel, cn)
    sig100.up <- rhdf5::h5read(h5.file, "signature/sig100.up",
                        index = list(1:100, sel.idx) )
    sig100.dn <- rhdf5::h5read(h5.file, "signature/sig100.dn",
                        index = list(1:100, sel.idx) )
    ##head(sig100.up,2)

    ## combine up/down into one (unsigned GSEA test)
    gmt <- rbind(sig100.up, sig100.dn)
    gmt <- unlist(apply(gmt, 2, list),recursive=FALSE)
    names(gmt) <- cn[sel.idx]
    length(gmt)

    ##system.time( res <- fgsea::fgsea(gmt, fc, nperm=10000))
    system.time( res <- fgsea::fgseaSimple(gmt, abs(fc), nperm=nperm))  ## really unsigned???
    dim(res)

    ## ---------------------------------------------------------------
    ## Combine correlation+GSEA by combined score (NES*rho)
    ## ---------------------------------------------------------------
    jj <- match( res$pathway, names(rho))
    res$rho  <- rho[jj]
    res$R2 <- rho[jj]**2
    res$score <- res$R2 * res$NES
    res <- res[order(res$score, decreasing=TRUE),]

    if(0) {
        res$rho.p <- cor.pvalue(res$rho, n=length(gg))
        res$meta.p  <- apply( res[,c("pval","rho.p")], 1, function(p) metap::sumz(p)$p)
        res <- res[order(res$meta.p),]
    }

    Matrix::head(res)
    return(res)
}

## ntop=1000;nsig=500;nperm=10000
#' @export
pgx.correlateSignature.matrix <- function(fc, refmat, nsig=100, ntop=1000, nperm=10000)
{
    ##
    ##
    ##
    ##

    if(is.null(names(fc))) stop("fc must have names")

    ## mouse... mouse...
    names(fc) <- toupper(names(fc))

    ## or instead compute correlation on top100 fc genes (read from file)
    ##refmat = PROFILES$FC
    rn <- rownames(refmat)
    cn <- colnames(refmat)

    ## ---------------------------------------------------------------
    ## Compute simple correlation between query profile and signatures
    ## ---------------------------------------------------------------
    gg <- intersect(rn, names(fc))
    fc1 <- sort(fc[gg])
    gg <- unique(names(c(Matrix::head(fc1,nsig), Matrix::tail(fc1,nsig))))
    ##gg <- intersect(names(fc),rn)
    ##gg <- intersect(gg,rn)
    G  <- refmat[gg,,drop=FALSE]
    dim(G)

    ## rank correlation??
    rG  <- apply( G[gg,], 2, rank, na.last="keep" )
    rfc <- rank( fc[gg], na.last="keep" )
    ##rho <- stats::cor(rG, rfc, use="pairwise")[,1]

    rG[is.na(rG)] <- 0  ## NEED RETHINK: are missing values to be treated as zero???
    rfc[is.na(rfc)] <- 0
    rho <- stats::cor(rG, rfc, use="pairwise")[,1]

    remove(G,rG,rfc)

    ## --------------------------------------------------
    ## test all signature on query profile using fGSEA
    ## --------------------------------------------------


    sel <- Matrix::head(names(sort(-abs(rho))),ntop)
    notx <- setdiff(sel,colnames(refmat))
    sel <- intersect(sel, colnames(refmat))
    X <- refmat[,sel,drop=FALSE]
    dim(X)
    X[is.na(X)] <- 0
    orderx <- apply(X,2,function(x) {
        idx=order(x);
        list(DN=head(idx,100),UP=rev(Matrix::tail(idx,100)))
    })
    sig100.dn <- sapply(orderx,"[[","DN")
    sig100.dn <- apply(sig100.dn, 2, function(i) rn[i])
    sig100.up <- sapply(orderx,"[[","UP")
    sig100.up <- apply(sig100.up, 2, function(i) rn[i])
    dim(sig100.dn)

    ## ---------------------------------------------------------------
    ## combine up/down into one (unsigned GSEA test)
    ## ---------------------------------------------------------------
    gmt <- rbind(sig100.up, sig100.dn)
    gmt <- unlist(apply(gmt, 2, list),recursive=FALSE)
    names(gmt) <- colnames(X)
    length(gmt)

    ##system.time( res <- fgsea::fgsea(gmt, fc, nperm=10000))
    suppressMessages( suppressWarnings(
        res <- fgsea::fgseaSimple(gmt, abs(fc), nperm=nperm)
    ))
    dim(res)

    ## ---------------------------------------------------------------
    ## Combine correlation+GSEA by combined score (NES*rho)
    ## ---------------------------------------------------------------
    jj <- match(res$pathway, names(rho))
    res$rho <- rho[jj]
    res$R2 <- rho[jj]**2
    res$score <- res$R2*res$NES
    res <- res[order(res$score, decreasing=TRUE),]

    if(0) {
        res$rho.p <- cor.pvalue(res$rho, n=length(gg))
        res$meta.p  <- apply( res[,c("pval","rho.p")], 1, function(p) metap::sumz(p)$p)
        res <- res[order(res$meta.p),]
    }

    Matrix::head(res)
    return(res)
}


chunk=100
#' @export
pgx.createCreedsSigDB <- function(gmt.files, h5.file, update.only=FALSE)
{

    h5exists <- function(h5.file, obj) {
        xobjs <- apply(rhdf5::h5ls(h5.file)[,1:2],1,paste,collapse="/")
        obj %in% gsub("^/|^//","",xobjs)
    }

    if(update.only && h5exists(h5.file, "data/matrix")) {
        X  <- rhdf5::h5read(h5.file, "data/matrix")
        rn <- rhdf5::h5read(h5.file,"data/rownames")
        cn <- rhdf5::h5read(h5.file,"data/colnames")
        rownames(X) <- rn
        colnames(X) <- cn
    } else {
        ##--------------------------------------------------
        ## make big FC signature matrix
        ##--------------------------------------------------
        F <- list()
        sig100.dn <- list()
        sig100.up <- list()
        cat("reading gene lists from",length(gmt.files),"gmt files ")
        i=1
        for(i in 1:length(gmt.files)) {
            if(!file.exists(gmt.files[i])) next()
            cat(".")
            try.error <- try( gmt <- read.gmt(gmt.files[i], add.source=TRUE) )
            if(class(try.error)=="try-error") next()
            ##gmt <- Matrix::head(gmt,30)  ## ONLY FOR TESTING

            j1 <- grep("-up ", names(gmt))
            j2 <- grep("-dn ", names(gmt))
            f1 <- lapply( gmt[j1], function(gg) {x=length(gg):1;names(x)=gg;x})
            f2 <- lapply( gmt[j2], function(gg) {x=-length(gg):-1;names(x)=gg;x})

            s1 <- gmt[j1]
            s2 <- gmt[j2]

            ff <- lapply(1:length(f1),function(i) c(f1[[i]],f2[[i]]))
            sig.names <- sub("-up","",names(f1))
            prefix <- gsub(".*/|single_|_perturbations|.gmt|_signatures","",gmt.files[i])
            sig.names <- paste0("[CREEDS:",prefix,"] ",sig.names)

            names(s1) <- names(s2) <- names(ff) <- sig.names
            sig100.up <- c(sig100.up, s1)
            sig100.dn <- c(sig100.dn, s2)
            F <- c(F, ff)
        }
        cat("\n")

        genes <- as.vector(unlist(sapply(F[],names)))
        genes <- sort(unique(toupper(genes)))
        length(genes)

        ## Filter out genes (not on known chromosomes...)
        gannot <- ngs.getGeneAnnotation(genes)
        table(!is.na(gannot$chr))
        sel <- which(!is.na(gannot$chr))
        genes <- sort(genes[sel])

        X <- lapply(F, function(x) x[match(genes,names(x))])
        X <- do.call(cbind, X)
        dim(X)
        rownames(X) <- genes
        remove(F)

        h5.file
        pgx.saveMatrixH5(X, h5.file, chunk=c(nrow(X),1))

        na100 <- rep(NA,100)
        msig100.up <- sapply(sig100.up, function(g) Matrix::head(c(intersect(g,genes),na100),100))
        msig100.dn <- sapply(sig100.dn, function(g) Matrix::head(c(intersect(g,genes),na100),100))

        if(!h5exists(h5.file, "signature")) rhdf5::h5createGroup(h5.file,"signature")
        rhdf5::h5write( msig100.up, h5.file, "signature/sig100.up")  ## can write list??
        rhdf5::h5write( msig100.dn, h5.file, "signature/sig100.dn")  ## can write list???
        remove(sig100.up,sig100.dn,msig100.up,msig100.dn)

        ## check NA!!! sometimes it is set to large negative
        if(1) {

            rhdf5::h5ls(h5.file)
            X  <- rhdf5::h5read(h5.file, "data/matrix")
            Matrix::head(X[,1])
            X[which(X < -999999)] <- NA
            Matrix::head(X[,1])
            dim(X)
            rhdf5::h5write( X, h5.file, "data/matrix")  ## can write list??
            rhdf5::h5closeAll()
        }

    }
    dim(X)


    ##--------------------------------------------------
    ## Precalculate t-SNE/UMAP
    ##--------------------------------------------------

    if(!update.only || !h5exists(h5.file, "clustering")) {

        X[is.na(X)] <- 0
        pos <- pgx.clusterBigMatrix(
            abs(X),  ## on absolute foldchange!!
            methods=c("pca","tsne","umap"),
            dims=c(2,3),
            reduce.sd = 2000,
            reduce.pca = 200 )
        names(pos)

        if(!h5exists(h5.file, "clustering")) rhdf5::h5createGroup(h5.file,"clustering")
        rhdf5::h5ls(h5.file)
        rhdf5::h5write( pos[["pca2d"]], h5.file, "clustering/pca2d")
        rhdf5::h5write( pos[["pca3d"]], h5.file, "clustering/pca3d")
        rhdf5::h5write( pos[["tsne2d"]], h5.file, "clustering/tsne2d")
        rhdf5::h5write( pos[["tsne3d"]], h5.file, "clustering/tsne3d")
        rhdf5::h5write( pos[["umap2d"]], h5.file, "clustering/umap2d")
        rhdf5::h5write( pos[["umap3d"]], h5.file, "clustering/umap3d")

    }

    rhdf5::h5closeAll()
    ## return(X)
}  ## end of pgx.createCreedsSigDB

#' @export
pgx.createSignatureDatabaseH5 <- function(h5.file, pgx.files, update.only=FALSE)
{

    h5exists <- function(h5.file, obj) {
        xobjs <- apply(rhdf5::h5ls(h5.file)[,1:2],1,paste,collapse="/")
        obj %in% gsub("^/|^//","",xobjs)
    }

    if(update.only && h5exists(h5.file, "data/matrix")) {
        X  <- rhdf5::h5read(h5.file, "data/matrix")
        rn <- rhdf5::h5read(h5.file,"data/rownames")
        cn <- rhdf5::h5read(h5.file,"data/colnames")
        rownames(X) <- rn
        colnames(X) <- cn
    } else {
        ##--------------------------------------------------
        ## make big FC signature matrix
        ##--------------------------------------------------
        F <- list()
        cat("[pgx.createSignatureDatabaseH5] reading FC from",length(pgx.files),"pgx files ")
        i=1
        for(i in 1:length(pgx.files)) {
            if(!file.exists(pgx.files[i])) next()
            cat(".")
            ##try.error <- try( load(pgx.files[i], verbose=0) )
            pgx <- try(local(get(load(pgx.files[i], verbose=0)))) ## override any name
            if("try-error" %in% class(pgx)) next()
            meta <- pgx.getMetaFoldChangeMatrix(pgx, what="meta")
            rownames(meta$fc) <- toupper(rownames(meta$fc))  ## mouse-friendly
            pgxfile <- gsub(".*[/]|[.]pgx$","",pgx.files[i])
            colnames(meta$fc) <- paste0("[",pgxfile,"] ",colnames(meta$fc))
            F[[ pgxfile ]] <- meta$fc
        }
        cat("\n")

        genes <- as.vector(unlist(sapply(F,rownames)))
        genes <- sort(unique(toupper(genes)))
        length(genes)
        F <- lapply(F, function(x) x[match(genes,rownames(x)),,drop=FALSE])
        X <- do.call(cbind, F)
        rownames(X) <- genes

        ## Filter out genes (not on known chromosomes...)
        genes <- rownames(X)
        gannot <- ngs.getGeneAnnotation(genes)
        table(is.na(gannot$chr))
        sel <- which(!is.na(gannot$chr))
        X <- X[sel,,drop=FALSE]
        remove(F)
    }
    dim(X)

    pgx.createSignatureDatabaseH5.fromMatrix(
      X = X,
      h5.file = h5.file,
      update.only = update.only
    )
    
}

#' @export
pgx.createSignatureDatabaseH5.fromMatrix <- function(h5.file, X, update.only=FALSE)
{

  file.exists(h5.file)
  if(file.exists(h5.file)) unlink(h5.file)
  ##chunk=c(nrow(X),1)
  pgx.saveMatrixH5(X, h5.file, chunk=c(nrow(X),1))
  
  ##--------------------------------------------------
  ## Calculate top100 gene signatures
  ##--------------------------------------------------
  cat("[pgx.createSignatureDatabaseH5.fromMatrix] Creating top-100 signatures...")
  
  if(!update.only || !h5exists(h5.file, "signature")) {
    ## X  <- rhdf5::h5read(h5.file, "data/matrix")
    rn <- rhdf5::h5read(h5.file,"data/rownames")
    cn <- rhdf5::h5read(h5.file,"data/colnames")
    rhdf5::h5ls(h5.file)
    
    dim(X)
    ##X <- X[,1:100]
    X[is.na(X)] <- 0
    orderx <- apply(X,2,function(x) {
      idx=order(x);
      list(DN=head(idx,100),UP=rev(Matrix::tail(idx,100)))
    })
    sig100.dn <- sapply(orderx,"[[","DN")
    sig100.dn <- apply(sig100.dn, 2, function(i) rn[i])
    sig100.up <- sapply(orderx,"[[","UP")
    sig100.up <- apply(sig100.up, 2, function(i) rn[i])
    
    if(!h5exists(h5.file, "signature")) rhdf5::h5createGroup(h5.file,"signature")
    rhdf5::h5write( sig100.dn, h5.file, "signature/sig100.dn")  ## can write list???
    rhdf5::h5write( sig100.up, h5.file, "signature/sig100.up")  ## can write list??
    
    remove(orderx)
        remove(sig100.dn)
    remove(sig100.up)
  }

    ##--------------------------------------------------
    ## Precalculate t-SNE/UMAP
    ##--------------------------------------------------
    if(!update.only || !h5exists(h5.file, "clustering")) {
      
      if(!h5exists(h5.file, "clustering")) rhdf5::h5createGroup(h5.file,"clustering")
      rhdf5::h5ls(h5.file)
      
      pos <- pgx.clusterBigMatrix(
        abs(X),  ## on absolute foldchange!!
        methods=c("pca","tsne","umap"),
        dims=c(2,3),
        reduce.sd = 2000,
        reduce.pca = 200 )
      names(pos)
      
      rhdf5::h5write( pos[["pca2d"]], h5.file, "clustering/pca2d")  ## can write list??
      rhdf5::h5write( pos[["pca3d"]], h5.file, "clustering/pca3d")  ## can write list??
      rhdf5::h5write( pos[["tsne2d"]], h5.file, "clustering/tsne2d")  ## can write list??
      rhdf5::h5write( pos[["tsne3d"]], h5.file, "clustering/tsne3d")  ## can write list??
      rhdf5::h5write( pos[["umap2d"]], h5.file, "clustering/umap2d")  ## can write list??
      rhdf5::h5write( pos[["umap3d"]], h5.file, "clustering/umap3d")  ## can write list??
    }
  
    rhdf5::h5closeAll()
    ## return(X)
}

#' @export
pgx.addEnrichmentSignaturesH5 <- function(h5.file, X=NULL, mc.cores=0,
                                          methods = c("gsea","rankcor")[2] )
{
    h5exists <- function(h5.file, obj) {
        xobjs <- apply(rhdf5::h5ls(h5.file)[,1:2],1,paste,collapse="/")
        obj %in% gsub("^/|^//","",xobjs)
    }

    if(is.null(X)) {
        X  <- rhdf5::h5read(h5.file, "data/matrix")
        rn <- rhdf5::h5read(h5.file,"data/rownames")
        cn <- rhdf5::h5read(h5.file,"data/colnames")
        rownames(X) <- rn
        colnames(X) <- cn
        X[which(X < -999999)] <- NA
    }

    ##---------------------- ONLY HALLMARK FOR NOW -----------------------
    ##sig100.dn <- rhdf5::h5read(h5.file, "signature/sig100.dn")
    ##sig100.up <- rhdf5::h5read(h5.file, "signature/sig100.up")
    G <- playdata::GSET_SPARSEG_XL
    dim(G)
    sel <- grep("HALLMARK|C[1-9]|^GO", rownames(G))
    sel <- grep("HALLMARK", rownames(G))
    length(sel)
    genes <- intersect(colnames(G),rownames(X))
    G <- G[sel,genes,drop=FALSE]
    X <- X[genes,]
    gmt <- apply(G, 1, function(x) colnames(G)[which(x!=0)])

    if(!h5exists(h5.file, "enrichment")) {
        rhdf5::h5createGroup(h5.file,"enrichment")
    }
    if(h5exists(h5.file, "enrichment/genesets")) {
        rhdf5::h5delete(h5.file, "enrichment/genesets")
    }
    rhdf5::h5write(names(gmt), h5.file, "enrichment/genesets")

    if("gsea" %in% methods) {
        cat("[pgx.addEnrichmentSignaturesH5] starting fGSEA for",length(gmt),"gene sets...")
        i=1
        F1 <- parallel::mclapply(colnames(X), function(i) {
            xi <- X[,i]
            xi <- xi[!is.na(xi)] 
            xi <- xi + 1e-3*rnorm(length(xi))
            res1 <- fgsea::fgseaSimple(gmt, xi, nperm=10000, nproc=mc.cores)
            r <- res1$NES
            names(r) <- res1$pathway
            r
        })
        F1 <- do.call(cbind, F1)
        F1 <- F1[match(names(gmt),rownames(F1)),,drop=FALSE]
        F1[is.na(F1)] <- 0
        colnames(F1) <- colnames(X)
        if(h5exists(h5.file, "enrichment/GSEA")) rhdf5::h5delete(h5.file, "enrichment/GSEA")
        rhdf5::h5write(F1, h5.file, "enrichment/GSEA")
    }

    if("rankcor" %in% methods) {
        F3 <- gset.rankcor(X, Matrix::t(G))$rho
        F3 <- F3[match(names(gmt),rownames(F3)),,drop=FALSE]
        F3[is.na(F3)] <- 0
        if(h5exists(h5.file, "enrichment/rankcor")) rhdf5::h5delete(h5.file, "enrichment/rankcor")
        rhdf5::h5write(F3, h5.file, "enrichment/rankcor")
    }

    rhdf5::h5ls(h5.file)
    rhdf5::h5closeAll()

    cat("[pgx.addEnrichmentSignaturesH5] done!\n")

}

#' @export
pgx.ReclusterSignatureDatabase <- function(h5.file, reduce.sd=1000, reduce.pca=100)
{


    h5exists <- function(h5.file, obj) {
        xobjs <- apply(rhdf5::h5ls(h5.file)[,1:2],1,paste,collapse="/")
        obj %in% gsub("^/|^//","",xobjs)
    }

    X  <- rhdf5::h5read(h5.file, "data/matrix")
    rn <- rhdf5::h5read(h5.file,"data/rownames")
    cn <- rhdf5::h5read(h5.file,"data/colnames")
    rownames(X) <- rn
    colnames(X) <- cn
    X[which(X < -999999)] <- NA

    ##--------------------------------------------------
    ## Precalculate t-SNE/UMAP
    ##--------------------------------------------------
    dim(X)

    if(!h5exists(h5.file, "clustering")) rhdf5::h5createGroup(h5.file,"clustering")

    pos <- pgx.clusterBigMatrix(
        abs(X),  ## on absolute foldchange!!
        methods = c("pca","tsne","umap"),
        dims = c(2,3),
        reduce.sd = reduce.sd,
        reduce.pca = reduce.pca )
    names(pos)

    rhdf5::h5write( pos[["pca2d"]], h5.file, "clustering/pca2d")  ## can write list??
    rhdf5::h5write( pos[["pca3d"]], h5.file, "clustering/pca3d")  ## can write list??
    rhdf5::h5write( pos[["tsne2d"]], h5.file, "clustering/tsne2d")  ## can write list??
    rhdf5::h5write( pos[["tsne3d"]], h5.file, "clustering/tsne3d")  ## can write list??
    rhdf5::h5write( pos[["umap2d"]], h5.file, "clustering/umap2d")  ## can write list??
    rhdf5::h5write( pos[["umap3d"]], h5.file, "clustering/umap3d")  ## can write list??
    rhdf5::h5closeAll()
}


##-------------------------------------------------------------------
## Pre-calculate geneset expression with different methods
##-------------------------------------------------------------------

#' @export
pgx.computeMultiOmicsGSE <- function(X, gmt, omx.type,
                                     method=NULL, center=TRUE)
{
    if(0) {
        omx.type <- c("MRNA","MIR")[1+grepl("^MIR",rownames(X))]
        table(omx.type)
        omx.type <- sample(c("MRNA","CNV"),nrow(X),replace=TRUE)
    }
    if(is.null(omx.type))
        omx.type <- gsub("[:=].*","",rownames(X))
    omx.types <- setdiff(unique(omx.type),c("MIR",""))
    omx.types

    sx <- list()
    for(tp in omx.types) {
        x1 <- X[which(omx.type==tp),]
        rownames(x1) <- sub(":.*","",rownames(x1))
        sx[[tp]] <- pgx.computeGeneSetExpression(x1, gmt, method=method, center=center)
        sx[[tp]] <- lapply(sx[[tp]],function(x) {
            rownames(x)=paste0(tp,"=",rownames(x))
            x
        })
    }

    ## concatenate all omx-types
    cx <- sx[[1]]
    for(j in 1:length(sx[[1]])) {
        cx[[j]] <- do.call(rbind, lapply(sx,"[[",j))
    }

    return(cx)
}

#' @export
pgx.computeGeneSetExpression <- function(X, gmt, method=NULL,
                                         min.size=10, center=TRUE)
{

    ALL.METHODS <- c("gsva","spearman","average")
    ALL.METHODS <- c("gsva","ssgsea","spearman","average")
    if(is.null(method))
        method <- ALL.METHODS
    ## this is important!!! centering on genes (GSVA does)
    if(center) {
        X <- X - rowMeans(X,na.rm=TRUE)
    }
    dim(X)

    gmt.size <- sapply(gmt, function(x) sum(x %in% rownames(X)))
    gmt <- gmt[ gmt.size >= min.size ]
    length(gmt)

    S <- list()
    if("gsva" %in% method) {
        S[["gsva"]] <- GSVA::gsva(X, gmt, method="gsva")
    }
    if("ssgsea" %in% method) {
        S[["ssgsea"]] <- GSVA::gsva(X, gmt, method="ssgsea", min.sz=1)
    }
    if(any(method %in% c("spearman","average"))) {
        gg <- rownames(X)
        G <- gmt2mat(gmt, bg=gg)
        if("spearman" %in% method) {
            ##rho <- stats::cor(as.matrix(G[gg,]), apply(X[gg,],2,rank))
            rho <- t(G[gg,]) %*% scale(apply(X[gg,],2,rank)) / sqrt(nrow(X)-1)
            rho[is.na(rho)] <- 0
            S[["spearman"]] <- rho
        }
        if("average" %in% method) {
            ##rho <- stats::cor(as.matrix(G[gg,]), apply(G[gg,],2,rank))
            avg.X <- t(G[gg,]) %*% X[gg,] / Matrix::colSums(G[gg,])
            avg.X[is.na(avg.X)] <- 0
            S[["average"]] <- avg.X
        }
    }

    ## compute meta score
    S1 <- lapply(S,function(x) apply(x,2,rank)) ## rank by sample
    S[["meta"]] <- scale(Reduce('+',S1)/length(S1))
    gs <- Reduce(intersect, lapply(S,rownames))
    S <- lapply(S, function(x) x[gs,])

    return(S)
}

##================================================================================
##========================= SIGDB H5 FUNCTIONS ===================================
##================================================================================

## IK 2.6.23: going to rewrite following functions without global
## SIGDB.DIR. For the moment seems these functions are duplicated in
## omicsplayground/components/board.connectivity

getSIGDB.DIR <- function() {
  if(!exists("SIGDB.DIR") && exists("FILESX")) {
    SIGDB.DIR <- c(FILESX,file.path(FILESX,"sigdb"))
  } else {
    SIGDB.DIR
  }
}


#' @export
sigdb.getConnectivityFullPath.DEPRECATED <- function(sigdb) {
    db.exists <- sapply(getSIGDB.DIR(), function(d) file.exists(file.path(d,sigdb)))
    db.exists
    db.dir <- names(which(db.exists))[1]
    db.dir
    file.path(db.dir, sigdb)
}

#' @export
sigdb.getConnectivityContrasts <- function(sigdb, path=NULL) {
    if(!is.null(path)) {
      sigdb <- file.path(path, sigdb)
    }
    if(!file.exists(sigdb)) return(NULL)
    rhdf5::h5read(sigdb, "data/colnames")
}

#' @export
sigdb.getConnectivityMatrix <- function(sigdb, select=NULL, genes=NULL, path=NULL)
{
    if(sigdb=="" || is.null(sigdb)) {
        warning("[getConnectivityMatrix] ***WARNING*** sigdb=",sigdb)
        return(NULL)
    }
    if(!is.null(path)) {
      sigdb <- file.path(path, sigdb)
    }
    if(!file.exists(sigdb)) {
        warning("[getConnectivityMatrix] ***WARNING*** file ",sigdb,"not found")
        return(NULL)
    }

    ##db.exists <- sapply( getSIGDB.DIR(), function(d) file.exists(file.path(d,sigdb)))
    if(grepl("csv$",sigdb)) {
      X <- read.csv(sigdb, row.names=1, check.names=FALSE)
      X <- as.matrix(X)
      X <- X[,colMeans(is.na(X)) < 0.99,drop=FALSE]  ## omit empty columns
      if(!is.null(genes)) X <- X[intersect(genes,rownames(X)),,drop=FALSE]
      if(!is.null(select)) X <- X[, intersect(select,colnames(X))]
    } else if(grepl("h5$",sigdb)) {
      cn <- rhdf5::h5read(sigdb, "data/colnames")
      rn <- rhdf5::h5read(sigdb, "data/rownames")
      rowidx <- 1:length(rn)
      colidx <- 1:length(cn)
      if(!is.null(genes)) rowidx <- match(intersect(genes,rn),rn)
      if(!is.null(select)) colidx <- match(intersect(select,cn),cn)      
      nr <- length(rowidx)
      nc <- length(colidx)
      dbg("[sigdb.getConnectivityMatrix] *WARNING* reading large H5 file:",nr,"x",nc,"")      
      X  <- rhdf5::h5read(sigdb, "data/matrix", index = list(rowidx,colidx) )
      rownames(X) <- rn[rowidx]
      colnames(X) <- cn[colidx]
    } else {
        cat("[getConnectivityMatrix] WARNING: could not retrieve matrix\n")
    }
    return(X)
}

#' @export
sigdb.getEnrichmentMatrix <- function(sigdb, select=NULL, path=NULL,
                                      which=c("gsea","rankcor")[1])
{
    if(sigdb=="" || is.null(sigdb)) {
        dbg("[getEnrichmentMatrix] ***WARNING*** sigdb=",sigdb)
        return(NULL)
    }
    if(!is.null(select)) {
        dbg("[getEnrichmentMatrix] length(select)=",length(select))
        dbg("[getEnrichmentMatrix] Matrix::head(select)=",head(select))
    }
    if(!grepl("h5$",sigdb)) {
        stop("getEnrichmentMatrix:: only for H5 database files")
        return(NULL)
    }
    if(!is.null(path)) {
      sigdb <- file.path(path, sigdb)
    }
    if(!file.exists(sigdb)) {
        warning("[sigdb.getEnrichmentMatrix] *WARNING* file ",sigdb,"not found")
        return(NULL)
    }

    h5exists <- function(h5.file, obj) {
        xobjs <- apply(rhdf5::h5ls(h5.file)[,1:2],1,paste,collapse="/")
        obj %in% gsub("^/|^//","",xobjs)
    }

    cn <- rhdf5::h5read(sigdb, "data/colnames")
    has.gs   <- h5exists(sigdb, "enrichment/genesets")
    if(!has.gs) {
        dbg("[getEnrichmentMatrix] WARNING: PGX object has no enrichment results")
        return(NULL)
    }

    rn <- rhdf5::h5read(sigdb, "enrichment/genesets")
    rowidx <- 1:length(rn)
    colidx <- 1:length(cn)
    if(!is.null(select)) colidx <- match(intersect(select,cn),cn)

    has.gsea <- h5exists(sigdb, "enrichment/GSEA")
    has.rankcor <- h5exists(sigdb, "enrichment/rankcor")
    Y <- NULL
    if(which == "gsea" && has.gsea) {
      Y  <- rhdf5::h5read(sigdb, "enrichment/GSEA", index = list(rowidx,colidx) )
    }
    if(which == "rankcor" && has.rankcor) {
      Y  <- rhdf5::h5read(sigdb, "enrichment/rankcor", index = list(rowidx,colidx) )
    }
    if(is.null(Y) || nrow(Y)==0) {
        return(NULL)
    }
    rownames(Y) <- rn[rowidx]
    colnames(Y) <- cn[colidx]

    return(Y)
}

#' @export
sigdb.getSignatureMatrix <- function(sigdb, path=NULL) {

    if(sigdb=="" || is.null(sigdb)) {
        dbg("[getEnrichmentMatrix] ***WARNING*** sigdb=",sigdb)
        return(NULL)
    }
    if(!grepl("h5$",sigdb)) {
        stop("getEnrichmentMatrix:: only for H5 database files")
    }
    if(!is.null(path)) {
      sigdb <- file.path(path, sigdb)
    }
    if(!file.exists(sigdb)) {
        warning("[sigdb.getEnrichmentMatrix] *WARNING* file ",sigdb,"not found")
        return(NULL)
    }

    rhdf5::h5ls(sigdb)
    cn <- rhdf5::h5read(sigdb, "data/colnames")
    dn <- rhdf5::h5read(sigdb, "signature/sig100.dn")
    up <- rhdf5::h5read(sigdb, "signature/sig100.up")
    colnames(dn) <- cn
    colnames(up) <- cn

    list(up=up, dn=dn)
}

##================================================================================
##=============================== END OF FILE ====================================
##================================================================================

