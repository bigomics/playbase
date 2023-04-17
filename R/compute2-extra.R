##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Title
#'
#' @param ngs
#' @param extra
#' @param lib.dir
#' @param sigdb
#'
#' @return
#' @export
#'
#' @examples
compute_extra <- function(ngs, extra=c("meta.go","deconv","infer","drugs", ## "graph",
                                       "connectivity","wordcloud"), lib.dir, sigdb=NULL) {

    timings <- c()
    libx.dir <- paste0(sub("/$","",lib.dir),'x')  ## ../libx yikes....
    message("[compute_extra] setting libx.dir =",libx.dir)

    extra <- intersect(extra, EXTRA.MODULES)
    if(length(extra)==0) {
        return(ngs)
    }

    ## detect if it is single or multi-omics
    single.omics <- !any(grepl("\\[",rownames(ngs$counts)))
    single.omics
    if(single.omics) {
        message(">>> computing extra for SINGLE-OMICS")
        rna.counts <- ngs$counts
    } else {
        message(">>> computing extra for MULTI-OMICS")
        data.type <- gsub("\\[|\\].*","",rownames(ngs$counts))
        jj <- which(data.type %in% c("gx","mrna"))
        length(jj)
        if(length(jj)==0) {
            stop("FATAL. could not find gx/mrna values.")
        }
        rna.counts <- ngs$counts[jj,]
        ##rownames(rna.counts) <- gsub(".*:|.*\\]","",rownames(rna.counts))
        is.logged <- ( min(rna.counts, na.rm=TRUE) < 0 ||
                       max(rna.counts, na.rm=TRUE) < 50 )
        if(is.logged) {
            message("expression data seems log. undoing logarithm")
            rna.counts <- 2**rna.counts
        }
    }

    if("meta.go" %in% extra) {
        message(">>> Computing GO core graph...")
        tt <- system.time({
            ngs$meta.go <- pgx.computeCoreGOgraph(ngs, fdr=0.20)
        })
        timings <- rbind(timings, c("meta.go", tt))
        message("<<< done!")
    }

    if("deconv" %in% extra) {
        message(">>> computing deconvolution")
        tt <- system.time({
            ngs <- compute_deconvolution(
                ngs, rna.counts=rna.counts,
                full=FALSE)
        })
        timings <- rbind(timings, c("deconv", tt))
        message("<<< done!")
    }

    if("infer" %in% extra) {
        message(">>> inferring extra phenotypes...")
        tt <- system.time({
            ngs <- compute_cellcycle_gender(ngs, rna.counts=rna.counts)
        })
        timings <- rbind(timings, c("infer", tt))
        message("<<< done!")
    }

    if("drugs" %in% extra) {
        ngs$drugs <- NULL  ## reset??
        cmap.dir <- file.path(libx.dir,"cmap")
        if(!dir.exists(cmap.dir)) {
            cmap.dir <- file.path(lib.dir,"cmap")  ## look for default lib
        }
        if(!dir.exists(cmap.dir)) {
            message("Warning:: missing CMAP files. Skipping drug connectivity analysis!")
        }
        dbg("[compute_extra] cmap.dir = ",cmap.dir)

        if(dir.exists(cmap.dir)) {

            message(">>> Computing drug activity enrichment...")
            tt <- system.time({
                ngs <- compute_drugActivityEnrichment(ngs)
            })
            timings <- rbind(timings, c("drugs", tt))

            message(">>> Computing drug sensitivity enrichment...")
            tt <- system.time({
                ngs <- compute_drugSensitivityEnrichment(ngs, cmap.dir)
            })
            timings <- rbind(timings, c("drugs-sx", tt))

            ## message(">>> Computing gene perturbation enrichment...")
            ## tt <- system.time({
            ##     ngs <- compute.genePerturbationEnrichment(ngs, lib.dir=cmap.dir)
            ## })
            ##timings <- rbind(timings, c("drugs-gene", tt))
        }
        message("<<< done!")
    }

    if("graph" %in% extra) {
        message(">>> computing OmicsGraphs...")
        tt <- system.time({
            ngs <- compute_omicsGraphs(ngs)
        })
        timings <- rbind(timings, c("graph", tt))
        message("<<< done!")
    }

    if("wordcloud" %in% extra) {
        message(">>> computing WordCloud statistics...")
        tt <- system.time({
            res <- pgx.calculateWordCloud(ngs, progress=NULL, pg.unit=1)
        })
        timings <- rbind(timings, c("wordcloud", tt))
        ngs$wordcloud <- res
        remove(res)
        message("<<< done!")
    }

    if("connectivity" %in% extra) {
        message(">>> computing connectivity scores...")

        ## ngs$connectivity <- NULL  ## clean up
        if(is.null(sigdb)) {
            lib.dir2 <- unique(c(lib.dir,libx.dir))  ### NEED BETTER SOLUTION!!!
            sig.dir <- c(SIGDB.DIR,lib.dir2)
            sigdb <- dir(sig.dir, pattern="^sigdb-.*h5$", full.names=TRUE)
            sigdb
        }

        db <- sigdb[1]
        for(db in sigdb) {
            if(file.exists(db)) {
                ntop = 10000
                ntop = 1000
                message("computing connectivity scores for ",db)
                ## in memory for many comparisons
                meta = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
                inmemory <- ifelse(ncol(meta$fc)>50,TRUE,FALSE)
                inmemory
                tt <- system.time({
                    scores <- pgx.computeConnectivityScores(
                        ngs, db, ntop=ntop, contrasts=NULL,
                        remove.le=TRUE, inmemory=inmemory )
                })
                timings <- rbind(timings, c("connectivity", tt))

                db0 <- sub(".*/","",db)
                ngs$connectivity[[db0]] <- scores
                remove(scores)
            }
        }
        names(ngs$connectivity)
    }

    ##------------------------------------------------------
    ## pretty collapse all timings
    ##------------------------------------------------------
    ##timings0 <- do.call(rbind, timings)
    timings <- as.matrix(timings)
    rownames(timings) <- timings[,1]
    timings0 <- apply(as.matrix(timings[,-1,drop=FALSE]),2,as.numeric)
    if(nrow(timings)==1) {
        timings0 <- matrix(timings0,nrow=1)
        colnames(timings0) <- colnames(timings)[-1]
        rownames(timings0) <- rownames(timings)
    }
    rownames(timings0) <- rownames(timings)
    timings0 <- apply( timings0, 2, function(x) tapply(x,rownames(timings0),sum))
    if(is.null(nrow(timings0))) {
        cn <- names(timings0)
        rn <- unique(rownames(timings))
        timings0 <- matrix(timings0, nrow=1)
        colnames(timings0) <- cn
        rownames(timings0) <- rn[1]
    }
    rownames(timings0) <- paste("[extra]",rownames(timings0))

    ngs$timings <- rbind(ngs$timings, timings0)
    message("<<< done!")

    return(ngs)
}

## -------------- deconvolution analysis --------------------------------

#' Title
#'
#' @param ngs
#' @param rna.counts
#' @param full
#'
#' @return
#' @export
#'
#' @examples
compute_deconvolution <- function(ngs, rna.counts=ngs$counts, full=FALSE) {

    ## list of reference matrices
    refmat <- list()
    #readSIG <- function(f) read.csv(file.path(lib.dir,"sig",f), row.names=1, check.names=FALSE)
    refmat[["Immune cell (LM22)"]] <- playbase::LM22 # read.csv(file.path(lib.dir,"sig/LM22.txt"),sep="\t",row.names=1)
    refmat[["Immune cell (ImmProt)"]] <- playbase::IMMPROT_SIGNATURE1000 #readSIG("immprot-signature1000.csv")
    refmat[["Immune cell (DICE)"]] <- playbase::DICE_SIGNATURE1000 #readSIG("DICE-signature1000.csv")
    refmat[["Immune cell (ImmunoStates)"]] <- playbase::IMMUNOSTATES_MATRIX #readSIG("ImmunoStates_matrix.csv")
    refmat[["Tissue (HPA)"]]       <- playbase::RNA_TISSUE_MATRIX #readSIG("rna_tissue_matrix.csv")
    refmat[["Tissue (GTEx)"]]      <- playbase::GTEX_RNA_TISSUE_TPM #readSIG("GTEx_rna_tissue_tpm.csv")
    refmat[["Cell line (HPA)"]]    <- playbase::HPA_RNA_CELLINE #readSIG("HPA_rna_celline.csv")
    refmat[["Cell line (CCLE)"]]   <- playbase::CCLE_RNA_CELLINE #readSIG("CCLE_rna_celline.csv")
    refmat[["Cancer type (CCLE)"]] <- playbase::CCLE_RNA_CANCERTYPE #readSIG("CCLE_rna_cancertype.csv")

    ## list of methods to compute
    ##methods = DECONV.METHODS
    methods = c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor","CIBERSORT","EPIC")
    ## methods <- c("NNLM","cor")

    if(full==FALSE) {
        ## Fast methods, subset of references
        sel = c("Immune cell (LM22)","Immune cell (ImmunoStates)",
                "Immune cell (DICE)","Immune cell (ImmProt)",
                "Tissue (GTEx)","Cell line (HPA)","Cancer type (CCLE)")
        refmat <- refmat[intersect(sel,names(refmat))]
        methods <- c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor")
    }

    ##counts <- ngs$counts
    counts <- rna.counts
    rownames(counts) <- toupper(ngs$genes[rownames(counts),"gene_name"])
    res <- pgx.multipleDeconvolution(counts, refmat=refmat, method=methods)

    ngs$deconv <- res$results
    if(!is.null(res$timings)) {
        rownames(res$timings) <- paste0("[deconvolution]",rownames(res$timings))
        res$timings
        ngs$timings <- rbind(ngs$timings, res$timings)
    }
    remove(refmat)
    remove(res)

    return(ngs)
}

## -------------- infer sample characteristics --------------------------------

#' Title
#'
#' @param ngs
#' @param rna.counts
#'
#' @return
#' @export
#'
#' @examples
compute_cellcycle_gender <- function(ngs, rna.counts=ngs$counts)
{
    pp <- rownames(rna.counts)
    is.mouse = (mean(grepl("[a-z]",gsub(".*:|.*\\]","",pp))) > 0.8)
    is.mouse
    if(!is.mouse) {
        if(1) {
            message("estimating cell cycle (using Seurat)...")
            ngs$samples$cell.cycle <- NULL
            ngs$samples$.cell.cycle <- NULL
            ##counts <- ngs$counts
            counts <- rna.counts
            rownames(counts) <- toupper(ngs$genes[rownames(counts),"gene_name"])
            res <- try(pgx.inferCellCyclePhase(counts) )  ## can give bins error
            if(class(res)!="try-error") {
                ngs$samples$.cell_cycle <- res
                table(ngs$samples$.cell_cycle)
            }
        }
        if(!(".gender" %in% colnames(ngs$samples) )) {
            message("estimating gender...")
            ngs$samples$.gender <- NULL
            X <- log2(1+rna.counts)
            gene_name <- ngs$genes[rownames(X),"gene_name"]
            ngs$samples$.gender <- pgx.inferGender( X, gene_name )
            table(ngs$samples$.gender)
        } else {
            message("gender already estimated. skipping...")
        }
        Matrix::head(ngs$samples)
    }
    return(ngs)
}


#' Title
#'
#' @param ngs
#'
#' @return
#' @export
#'
#' @examples
compute_drugActivityEnrichment <- function(ngs) {

    ## -------------- drug enrichment
    # get drug activity databases
    ref.db <- grep('L1000_*',data(package='playbase')$results[,'Item'],value=T)
    ref.db <- ref.db[ref.db != 'L1000_REPRURPOSING_DRUGS']

    if(length(ref.db)==0) {
        message("[compute_drugActivityEnrichment] Warning:: missing drug activity database")
        return(ngs)
    }
    names(ref.db) <- sub("-","/",gsub("_.*","",ref.db))

    for(i in 1:length(ref.db)) {

        f <- ref.db[i]
        message("[compute_drugActivityEnrichment] reading L1000 reference: ",f)
        X <- get(f)
        xdrugs <- gsub("[_@].*$","",colnames(X))
        ndrugs <- length(table(xdrugs))
        message("number of profiles: ",ncol(X))
        message("number of drugs: ",ndrugs)
        is.drug <- grepl("activity|drug|ChemPert",f,ignore.case=TRUE)

        NPRUNE=250
        fname <- names(ref.db)[i]
        out1 <- pgx.computeDrugEnrichment(
            ngs, X, xdrugs, methods=c("GSEA","cor"),
            nmin=3, nprune=NPRUNE, contrast=NULL )

        if(is.null(out1)) {
            message("[compute_drugActivityEnrichment] WARNING:: pgx.computeDrugEnrichment failed!")
            next()
        }

        ## --------------- attach annotation
        annot0 <- NULL
        if(is.drug) {
            annot0 <- playbase::L1000_REPRURPOSING_DRUGS
            annot0$drug <- annot0$pert_iname
            rownames(annot0) <- annot0$pert_iname
        } else {
            ## gene perturbation OE/LIG/SH
            dd <- rownames(out1[["GSEA"]]$X)
            d1 <- dd
            d2 <- sub("-.*","",dd)
            annot0 <- data.frame(drug=dd, moa=d1, target=d2)
            rownames(annot0) <- dd
        }

        if(1) {
            annot0 <- annot0[match(rownames(out1[["GSEA"]]$X),rownames(annot0)),]
            rownames(annot0) <- rownames(out1[["GSEA"]]$X)
            Matrix::head(annot0)
        }

        ## --------------- attach results to object
        db <- names(ref.db)[i]
        ## ngs$drugs[["activity/L1000"]]  <- res.mono[["GSEA"]]
        ## ngs$drugs[["activity/L1000"]][["annot"]] <- annot0[,c("drug","moa","target")]
        ngs$drugs[[db]] <- out1[["GSEA"]]
        ngs$drugs[[db]][["annot"]] <- annot0[,c("drug","moa","target")]
        ngs$drugs[[db]][["clust"]] <- out1[["clust"]]
        ngs$drugs[[db]][["stats"]] <- out1[["stats"]]
    }

    remove(X)
    remove(xdrugs)
    return(ngs)
}

#' Title
#'
#' @param ngs
#' @param cmap.dir
#'
#' @return
#' @export
#'
#' @examples
compute_drugSensitivityEnrichment <- function(ngs, cmap.dir)
{

    ref.db <- dir(cmap.dir, pattern='sensitivity.*rds$')
    ref.db
    if(length(ref.db)==0) {
        message("[compute_drugSensitivityEnrichment] Warning:: missing drug sensitivity database")
        return(ngs)
    }
    names(ref.db) <- sub("-","/",gsub("_.*","",ref.db))
    ref.db
    ref <- ref.db[1]
    for(i in 1:length(ref.db)) {
        ref <- ref.db[i]
        ##X <- readRDS(file=file.path(cmap.dir,"drugSX-GDSC-t25-g1000.rds"))
        ##X <- readRDS(file=file.path(cmap.dir,"drugSX-CTRPv2-t25-g1000.rds"))
        X <- readRDS(file=file.path(cmap.dir,ref))
        xdrugs <- gsub("[@_].*$","",colnames(X))
        length(table(xdrugs))
        dim(X)

        NPRUNE=-1
        NPRUNE=250
        out1 <- pgx.computeDrugEnrichment(
            ngs, X, xdrugs, methods=c("GSEA","cor"),
            nmin=10, nprune=NPRUNE, contrast=NULL )

        if(!is.null(out1)) {
            ## attach annotation
            db <- sub("-.*","",ref)
            annot0 <- read.csv(file.path(cmap.dir,paste0(db,"-drugs.csv")))
            Matrix::head(annot0)
            rownames(annot0) <- annot0$drug
            annot0 <- annot0[match(rownames(out1[["GSEA"]]$X),rownames(annot0)),]
            rownames(annot0) <- rownames(out1[["GSEA"]]$X)
            dim(annot0)

            s1 <- names(ref.db)[i]
            ngs$drugs[[s1]] <- out1[["GSEA"]]
            ngs$drugs[[s1]][["annot"]] <- annot0[,c("moa","target")]
            ngs$drugs[[s1]][["clust"]] <- out1[["clust"]]
            ngs$drugs[[s1]][["stats"]] <- out1[["stats"]]
        }

    } ## end of for rr

    names(ngs$drugs)

    remove(X)
    remove(xdrugs)
    return(ngs)
}

## ------------------ Omics graphs --------------------------------

#' Title
#'
#' @param ngs
#'
#' @return
#' @export
#'
#' @examples
compute_omicsGraphs <- function(ngs) {
    ## gr1$layout <- gr1$layout[igraph::V(gr1)$name,]  ## uncomment to keep entire layout
    ngs$omicsnet <- pgx.createOmicsGraph(ngs)
    ngs$pathscores <- pgx.computePathscores(ngs$omicsnet, strict.pos=FALSE)

    ## compute reduced graph
    ngs$omicsnet.reduced <- pgx.reduceOmicsGraph(ngs)
    ngs$pathscores.reduced <- pgx.computePathscores(ngs$omicsnet.reduced, strict.pos=FALSE)
    ##save(ngs, file=rda.file)
    return(ngs)
}

