##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' @export
pgx.load <- function(file, verbose=0) {
  local(get(load(file, verbose=verbose)))
}

#' @export
ngs.save <- function(ngs, file, update.date=TRUE, light=TRUE, system=FALSE) {
  message("warning: ngs.save() is deprecated. please use pgx.save()")
  pgx.save(ngs, file=file, update.date=update.date, light=light, system=system)
}

#' @export
pgx.save <- function(pgx, file, update.date=TRUE, light=TRUE, system=FALSE) {

    if(update.date||is.null(pgx$date)) pgx$date <- Sys.Date()

    if(light) {
        ## ------- make a light version
        pgx$gx.meta$outputs <- NULL
        pgx$gset.meta$outputs <- NULL
        pgx$model.parameters$efit <- NULL
        pgx$gmt.all <- NULL
        ##pgx$families <- NULL
        pgx$collections <- NULL
        ## pgx$counts <- NULL
        pgx$gset.meta$matrices <- NULL
    }
    if(system==FALSE) {
        ## remove system (is big...)
        pgx$omicsnet <- NULL
        pgx$omicsnet.reduced <- NULL
    }
    sort(sapply(pgx, object.size)) / 1e9
    sum(sapply(pgx, object.size)) / 1e9

    cat(">>> saving PGX file to",file,"\n")
    file <- iconv(file, from = '', to = 'ASCII//TRANSLIT')
    save(pgx, file=file)
}


##filter.get="omicsplayground";from=NULL;to=NULL;unique=TRUE
#' @export
pgx.parseAccessLogs <- function(logs.dir, from=NULL, to=NULL,
                                filter.get=NULL, unique = TRUE)
{

    logs.dir <- logs.dir[dir.exists(logs.dir)]
    if(length(logs.dir)==0) return(NULL)
    logs.dir
    access.files <- lapply(logs.dir, dir, pattern="access.log",
                           full.names=TRUE, recursive=TRUE )
    access.files
    access.files <- unlist(access.files)
    ##access.files <- grep("shinyproxy",access.files,value=TRUE,invert=TRUE)
    access.files

    access.logs <- lapply(access.files, function(f)
        suppressMessages(suppressWarnings(try(read.table(f)))))
    ##access.logs <- lapply(access.files, function(f)
    ##    suppressMessages(suppressWarnings(try(data.table::fread(f,sep=" ")))))
    access.logs <- access.logs[sapply(access.logs,class)!="try-error"]
    access.logs <- access.logs[sapply(access.logs,nrow)>0]
    length(access.logs)
    if(length(access.logs)==0) return(NULL)

    i=3
    for(i in 1:length(access.logs)) {
        df <- data.frame(access.logs[[i]])
        cols <- c("V1","V4","V6")
        if(mean(grepl(":80$",df[,"V1"])) > 0.9) {
            cols <- c("V2","V5","V7")
        }
        access.logs[[i]] <- df[,cols]
        colnames(access.logs[[i]]) <- c("ip","date","get")
    }
    remove(df)

    ## Filter access log
    acc <- do.call(rbind, access.logs)
    dim(acc)

    if(!is.null(filter.get)) {
        sel <- grep(filter.get,acc[,"get"])
        acc <- acc[sel,]
    }
    dim(acc)
    Matrix::head(acc)

    ## Extract visiting period
    ## if the operating system is not windows set the timezone to LC_TIME
    if(Sys.info()["sysname"] != "Windows") {
        Sys.setlocale("LC_TIME","en_US.UTF-8")
    }
    ##Sys.setlocale("LC_TIME","C") ## just to make sure
    acc$date <- gsub("[:].*|\\[","",as.character(acc[,"date"]))
    acc$date <- as.Date(acc$date, format = "%d/%b/%Y")
    acc <- acc[order(acc$date),]

    from.date <- Matrix::head(acc$date,1)
    to.date <- Matrix::tail(acc$date,1)
    from.to <- paste(from.date,"-",to.date)
    from.to

    ## Extract IP
    acc.ip <- as.character(acc[,"ip"])
    ##loc <- ip_api(unique(acc.ip))
    unique.ip <- unique(acc.ip)
    ## unique.hostname <- ip_to_hostname(unique.ip)  ## takes loooonnnggg... time!!
    ## names(unique.hostname)  <- unique.ip

    ## create lookup-table for IP to country

    file <- system.file("extdata","GeoLite2-Country.mmdb", package = "rgeolocate")
    loc  <- rgeolocate::maxmind(unique.ip, file, c("country_code","country_name"))
    loc$ip <- unique.ip
    ##loc <- rgeolocate::maxmind(ip, file, c("country_code", "country_name", "city_name"))
    loc$country_name[which(loc$ip %in% c("127.0.0.1"))] <- "<local.ip>"
    ##loc$country_code[which(loc$ip %in% c("127.0.0.1"))] <- "<local.ip>"
    loc$country_name[is.na(loc$country_name)] <- "(unknown)"

    country_codes <- unique(loc$country_code)
    names(country_codes) <- loc[match(country_codes,loc$country_code),"country_name"]
    country_codes["(unknown)"] = "(unknown)"

    ## now map IP to country_code
    acc$country_code <- loc$country_code[match(acc.ip,loc$ip)]
    acc$country_name <- loc$country_name[match(acc.ip,loc$ip)]
    acc$country_code[is.na(acc$country_code)] <- "(unknown)"
    acc$country_name[is.na(acc$country_name)] <- "(unknown)"
    Matrix::tail(sort(table(acc$country_code)),40)

    if(0) {
        getDodgy <- function(acc0,n=100) {
            ii <- grep("omicsplayground",acc0[,"get"],invert=TRUE)
            ii <- ii[which(nchar(as.character(acc0[ii,"get"])) > 20)]
            Matrix::head(acc0[ii,],n=n)
        }
        Matrix::tail(sort(table(acc$country_code)),40)
        ii <- which(acc.cc=="US")
        getDodgy(acc[ii,],100)
        ii <- which(acc.cc=="CN")
        getDodgy(acc[ii,],100)
    }

    ## cumulative table: counting total number of hits since start
    acc$days <- acc$date - from.date + 1
    ndays <- max(acc$days)
    ncountries <- length(unique(acc$country_code))
    M <- matrix(0, nrow=ncountries, ncol=ndays)
    rownames(M) <- sort(unique(acc$country_code))
    d = 1
    for(d in 1:ndays) {
        jj <- which(acc$days <= d)
        if(unique) jj <- jj[which(!duplicated(acc$ip[jj]))]  ## unique per day
        tt <- table(as.character(acc$country_code[jj]))
        tt <- tt[match(rownames(M),names(tt))]
        tt[is.na(tt)] <- 0
        M[,d] <- tt
    }
    M[is.na(M)] <- 0
    colnames(M) <- as.character(from.date + 1:ncol(M) -1 )

    ## weekly table: counting hits per week since start
    week.since <- seq(from.date, max(acc$date), 7)
    W <- matrix(0, nrow=ncountries, ncol=length(week.since)-1)
    rownames(W) <- sort(unique(acc$country_code))
    colnames(W) <- paste0("week.",1:ncol(W))
    i = 1
    acc0 <- acc
    ## acc0 <- acc[!duplicated(acc$ip),] ## unique, new visitors
    for(i in 1:(ncol(W)-1)) {
        jj <- which(acc0$date >= week.since[i] & acc0$date < week.since[i+1] )
        if(unique) jj <- jj[which(!duplicated(acc0$ip[jj]))]
        tt <- table(as.character(acc0$country_code[jj]))
        tt <- tt[match(rownames(W),names(tt))]
        tt[is.na(tt)] <- 0
        W[,i] <- tt
    }

    if(0) {
        barplot(Matrix::colSums(W,na.rm=TRUE),las=3, main="unique.IP")
    }

    ## final table
    if(unique) {
        tt <- table(as.character(loc$country_name))  ## unique IP
    } else {
        tt <- table(as.character(acc$country_name))  ## full access log
    }
    sum(tt)
    cc <- country_codes[names(tt)]
    names(cc) <- names(tt)
    df <- data.frame( country_name = names(tt),
                     country_code = cc,
                     count = (as.integer(tt)))
    df <- df[order(-df$count),]
    sum(df$count,na.rm=TRUE)
    M <- M[match(df$country_code, rownames(M)),]
    W <- W[match(df$country_code, rownames(W)),]

    res <- list(visitors=df, period=from.to, from=from.date, to=to.date,
                table=M, weekly=W, access=acc)
    return(res)
}

#' @export
h5exists <- function(h5.file, obj) {
    xobjs <- apply(rhdf5::h5ls(h5.file)[,1:2],1,paste,collapse="/")
    obj %in% gsub("^/|^//","",xobjs)
}

##h5.file="test.h5";chunk=100
#' @export
pgx.saveMatrixH5 <- function(X, h5.file, chunk=NULL )
{
    if(file.exists(h5.file)) unlink(h5.file)

    rhdf5::h5createFile(h5.file)
    rhdf5::h5createGroup(h5.file,"data")

    if(is.null(chunk)) {
        rhdf5::h5write( X, h5.file, "data/matrix")
    } else {
        rhdf5::h5createDataset(
            h5.file, "data/matrix",
            c(nrow(X),ncol(X)),
            ##storage.mode = "integer",
            chunk = chunk,
            level = 7
        )
        rhdf5::h5write(
            X,
            file = h5.file,
            name = "data/matrix",
            index = list(1:nrow(X),1:ncol(X))
        )
    }
    rhdf5::h5write( rownames(X), h5.file, "data/rownames")
    rhdf5::h5write( colnames(X), h5.file, "data/colnames")

    rhdf5::h5closeAll()
}

#' @export
pgx.readMatrixH5 <- function(h5.file, select=NULL, rows=NULL) {
    if(is.null(select) && is.null(rows)) {
        X  <- rhdf5::h5read(h5.file, "data/matrix")
        rn <- rhdf5::h5read(h5.file,"data/rownames")
        cn <- rhdf5::h5read(h5.file,"data/colnames")
    }
    if(!is.null(select) || !is.null(rows)) {
        X  <- rhdf5::h5read(h5.file, "data/matrix", index = list(rows, select))
        rn <- rhdf5::h5read(h5.file,"data/rownames", index = list(rows))
        cn <- rhdf5::h5read(h5.file,"data/colnames", index = list(select))
    }
    rownames(X) <- rn
    colnames(X) <- cn
    X[which(X < -999999)] <- NA
    as.matrix(X)
}

#' @export
pgx.readOptions <- function(file = "./OPTIONS") {
    if(!file.exists(file)) return(NULL)
    opt <- read.table(file, sep="=", row.names=1)
    opt <- gsub("^[ ]*|[ ]*$","",apply(opt,1,c)) ## strip leading/post spaces
    opt <- sapply(opt,list)
    opt <- sapply(opt,strsplit,split="[;,]")
    is.bool <- sapply(opt, function(x) all(tolower(x) %in% c("true","false")))
    is.bool
    opt[is.bool] <- sapply(opt[is.bool], function(x) tolower(x) %in% c("true"))
    names(opt) <- trimws(names(opt))
    opt
}

#' @export
pgx.readDatasetProfiles <- function(pgx.dir, file="datasets-allFC.csv",
                                     verbose=TRUE)
{
    if(!dir.exists(pgx.dir)) {
        stop(paste("[readDatasetProfiles1] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }
    fn <- file.path(pgx.dir,file)
    fn
    if(!file.exists(fn)) {
        stop("FATAL : could not find profiles matrix. please create first with initDatasetFolder().\n")
        ## pgx.updateDatasetProfiles(pgx.dir, file=file)
        return()
    } else {
        if(verbose) message("[readDatasetProfiles1] Found existing dataset profiles matrix")
    }

    allFC <- fread.csv(file=file.path(pgx.dir, file), row.names=1, check.names=FALSE)
    allFC <- as.matrix(allFC)
    if(verbose) message("[readDatasetProfiles1] dataset profiles matrix : dim=",dim(allFC))
    dim(allFC)
    return(allFC)
}

#' @export
pgx.scanInfoFile <- function(pgx.dir, file="datasets-info.csv", force=FALSE, verbose=TRUE)
{
  INITDATASETFOLDER = TRUE
  
  pgx.files <- dir(pgx.dir, pattern="[.]pgx$")
  if(length(pgx.files)==0) return(NULL)  ## no files!

  ## only run pgx.initDatasetFolder if pgx are changed

  if(!dir.exists(pgx.dir)) {
        stop(paste("[initDatasetFolder] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }
    dbg("[pgx.initDatasetFolder] *** called ***")
  
    ## all public datasets
    pgx.dir <- pgx.dir[1]  ## only one folder!!!
    pgx.files <- dir(pgx.dir, pattern="[.]pgx$")
    pgx.files <- sub("[.]pgx$","",pgx.files) ## strip pgx

    dbg("[pgx.initDatasetFolder] length(pgx.files) = ",length(pgx.files))
    if(length(pgx.files)==0) {
        allfc.file1 <- file.path(pgx.dir, allfc.file)
        info.file1  <- file.path(pgx.dir, info.file)
        if(delete.old) {
          unlink(info.file1)
          unlink(allfc.file1)
        }
        return(NULL)
    }
  
    ##----------------------------------------------------------------------
    ## If an allFC file exists
    ##----------------------------------------------------------------------

    allfc.file1 <- file.path(pgx.dir, allfc.file)
    has.fc <- file.exists(allfc.file1)

    info.file1 <- file.path(pgx.dir, info.file)
    has.info <- file.exists(info.file1)

    ##----------------------------------------------------------------------
    ## If an allFC file exits, check if it is done for all PGX files
    ##----------------------------------------------------------------------

    pgxinfo <- NULL
    pgx.missing0 <- pgx.files
    pgx.missing1 <- pgx.files
    pgx.delete0 <- c()
    pgx.delete1 <- c()

    allFC <-NULL
    if(!force && has.fc) {
        if(verbose) message("[initDatasetFolder] checking which pgx files already done in allFC...")
        allFC <- data.table::fread(allfc.file1,check.names=FALSE,nrows=1) ## HEADER!!!        
        pgx.done <- gsub("^\\[|\\].*","",colnames(allFC)[-1])
        pgx.missing0 <- setdiff(pgx.missing0, pgx.done)
        pgx.delete0  <- setdiff(pgx.done, pgx.files)
        allFC <-NULL        
    }

    if(!force && has.info) {
        if(verbose) message("[initDatasetFolder] checking which pgx files already in PGX info...")
        ## do not use fread! quoting bug
        pgxinfo = read.csv(info.file1, stringsAsFactors=FALSE, row.names=1, sep=',')        
        pgxinfo.files <- unique(sub(".pgx$","",pgxinfo$dataset))
        pgx.missing1 <- setdiff(pgx.missing1, pgxinfo.files)
        pgx.delete1  <- setdiff(pgxinfo.files, pgx.files)
    }

    ##----------------------------------------------------------------------
    ## Check if it is done for all PGX files
    ##----------------------------------------------------------------------

    ## files to be done either for allFC or missing in INFO
    pgx.missing <- unique(c(pgx.missing0, pgx.missing1))
    pgx.delete  <- unique(c(pgx.delete0, pgx.delete1))

    pgxinfo.changed = FALSE
    pgxfc.changed = FALSE

    ## these pgx need forced update
    if(!is.null(new.pgx)) {
       new.pgx <- sub(".pgx$","",new.pgx)
       new.pgx <- intersect(new.pgx, pgx.files)  ## only existing pgx
       pgx.delete  <- union(pgx.delete, new.pgx)
       pgx.missing <- union(pgx.missing, new.pgx)
       pgx.missing0 <- union(pgx.missing0, new.pgx)
       pgx.missing1 <- union(pgx.missing1, new.pgx)
       sel1 <- which(pgxinfo$dataset %in% new.pgx)
       if(length(sel1)) {
           pgxinfo <- pgxinfo[-sel1,]
           pgxinfo.changed <- TRUE       
       }
    }
  
    if(length(pgx.missing)==0 && length(pgx.delete)==0) {
        if(verbose) message("[initDatasetFolder] no update required. use FORCE=1 for forced update.")
        return(INITDATASETFOLDER = FALSE)
  }

  return(INITDATASETFOLDER)

  ## before reading the info file, we need to update for new files
  pgx.initDatasetFolder(pgx.dir, force=force, verbose=TRUE)  

  
}

#' @export
pgx.initDatasetFolder <- function( pgx.dir,
                                  allfc.file = "datasets-allFC.csv",
                                  info.file = "datasets-info.csv",
                                  force = FALSE, delete.old = FALSE,
                                  new.pgx = NULL,
                                  verbose = TRUE)
{
    ##
    ## Initialize file information file for SINGLE folder
    ##
    ##
    
    if(verbose) message("[initDatasetFolder] folder has ",length(pgx.missing)," new PGX files")
    if(verbose) message("[initDatasetFolder] info-file has ",length(pgx.delete)," old items")

    ##----------------------------------------------------------------------
    ## Reread allFC file. Before we only read the header.
    ##----------------------------------------------------------------------
    allFC <-NULL
    if(!force && file.exists(allfc.file1) && length(pgx.missing)>0) {
        ##allFC <- read.csv(allfc.file1,row.names=1,check.names=FALSE)
        allFC <- fread.csv(allfc.file1,row.names=1,check.names=FALSE)
    }
    dim(allFC)

    ## these pgx need forced update, so remove
    if(!is.null(allFC) && !is.null(new.pgx)) {
       new.pgx <- sub(".pgx$","",new.pgx)
       allfc.pgx <- gsub("^\\[|\\].*","",colnames(allFC))
       sel2 <- which(allfc.pgx %in% new.pgx)
       if(length(sel2)) {
          allFC <- allFC[,-sel2,drop=FALSE]
          pgxfc.changed <- TRUE         
       }
    }  
  
    ##----------------------------------------------------------------------
    ## For all new PGX files, load the PGX file and get the meta FC
    ## matrix.
    ##----------------------------------------------------------------------

    info.cols <- NULL
    missing.FC <- list()
    message("[initDatasetFolder] missing pgx = ",paste(pgx.missing,collapse=" "))
    pgxfile = pgx.missing[1]
  
    ngs <- NULL
    for(pgxfile in pgx.missing) {

        cat(".")
        pgxfile1 <- file.path(pgx.dir,pgxfile)
        pgxfile1 <- paste0(sub("[.]pgx$","",pgxfile1),".pgx")
        ##try.error <- try( load(file.path(pgx.dir,pgxfile),verbose=0) )
        pgx <- try(local(get(load(pgxfile1, verbose = 0)))) ## override any name

        if("try-error" %in% class(pgx)) {
            message(paste("[initDatasetFolder] ERROR in loading PGX file:",pgxfile1,". skipping\n"))
            next()
        }

        if(!pgx.checkObject(pgx)) {
            message(paste("[initDatasetFolder] INVALID PGX object",pgxfile,". Skipping"))
            next()
        }

        ##---------------------------------------------
        ## extract the meta FC matrix
        ##---------------------------------------------
        ## rownames(pgx$X) <- toupper(sub(".*:","",rownames(pgx$X)))
        if(pgxfile %in% pgx.missing0) {
          meta <- pgx.getMetaFoldChangeMatrix(pgx, what="meta")
          rownames(meta$fc) <- toupper(rownames(meta$fc))
          missing.FC[[pgxfile]] <- meta$fc
          pgxfc.changed <- TRUE
        }

        ##---------------------------------------------
        ## compile the info for update
        ##---------------------------------------------
        if(pgxfile %in% pgx.missing1) {
          pgx$name <- sub(".pgx$","",pgxfile)  ## force filename as name
          pgxinfo <- pgx.updateInfoPGX(pgxinfo, pgx)
          pgxinfo.changed <- TRUE
        }
    }

    pgx <- NULL

    ##----------------------------------------------------------------------
    ## Update the INFO meta file
    ##----------------------------------------------------------------------
    rownames(pgxinfo) <- NULL
    pgxinfo <- data.frame(pgxinfo, check.names=FALSE)

    ## remove unneccessary entries if forced.
    sel.delete <- which(!sub(".pgx$","",pgxinfo$dataset) %in% pgx.files)
    if(length(sel.delete) && delete.old) {
      pgxinfo <- pgxinfo[-sel.delete,,drop=FALSE]
      pgxinfo.changed <- TRUE
    }

    if(pgxinfo.changed) {
      if(verbose) message("[initDatasetFolder] writing updated PGX.INFO file to ",info.file1,"...")
      write.csv(pgxinfo, file = info.file1)
      Sys.chmod(info.file1, "0666")
    }

    ##----------------------------------------------------------------------
    ## Update the ALL.FC meta file
    ##----------------------------------------------------------------------

    ## remove unneccessary entries if forced.
    if(!is.null(allFC) && delete.old) {
      fc.done <- gsub("^\\[|\\].*","",colnames(allFC))
      sel.deleteFC <- which(!fc.done %in% pgx.files)
      if(length(sel.deleteFC)) {
        allFC <- allFC[, -sel.deleteFC, drop=FALSE]
        pgxfc.changed <- TRUE
      }
    }

    if(length(missing.FC)==0 && !pgxfc.changed) {
        ## no change in info
        dbg("[initDatasetFolder] allFC complete. no change needed.")
        return(NULL)
    }

    if(length(missing.FC)>0) {
      ## find most common genes
      all.gg <- toupper(as.character(unlist(sapply(missing.FC, rownames))))
      gg.tbl <- table(all.gg)
      table(gg.tbl)

      ## Conform the multiple metaFC matrices
      gg <- names(gg.tbl)
      length(gg)
      missing.FC <- lapply(missing.FC, function(x) {
        x <- x[match(gg,toupper(rownames(x))),,drop=FALSE]
        rownames(x) <- gg
        return(x)
      })

      ## append file name in front of contrast names
      id <- paste0("[",sub("[.]pgx","",names(missing.FC)),"]")
      id
      for(i in 1:length(missing.FC)) {
        colnames(missing.FC[[i]]) <- paste0(id[i]," ",colnames(missing.FC[[i]]))
      }
      allFC.new <- do.call(cbind, missing.FC)
      allFC.new <- as.matrix(allFC.new)

      if(is.null(allFC)) {
        allFC <- allFC.new
      } else {
        ## Add any new FC profiles to the existing allFC
        gg <- sort(unique(c(rownames(allFC), rownames(allFC.new))))
        j1 <- match(gg, rownames(allFC))
        j2 <- match(gg, rownames(allFC.new))
        allFC <- allFC[j1,,drop=FALSE]
        allFC.new <- allFC.new[j2,,drop=FALSE]
        allFC <- cbind(allFC, allFC.new)
        rownames(allFC) <- gg
      }

      ## restrict to 20000 genes
      allfc.sd <- apply(allFC, 1, sd, na.rm=TRUE)
      allfc.nna <- rowMeans(!is.na(allFC))
      jj <- Matrix::head( order(-allfc.sd * allfc.nna), 20000)
      allFC <- allFC[jj,,drop=FALSE]
      dim(allFC)
      pgxfc.changed <- TRUE
    }
  
    ## save modified allFC
    if(pgxfc.changed) {
      ## check for duplicates
      if(verbose) message("[initDatasetFolder] allFC changed. updating file to",allfc.file1,"...")
      allFC <- allFC[,!duplicated(colnames(allFC)),drop=FALSE]
      allFC <- allFC[,order(colnames(allFC)),drop=FALSE]
      AA <- data.frame(rownames=rownames(allFC), allFC, check.names=FALSE)
      ##write.csv(allFC, file=allfc.file1)
      data.table::fwrite(AA, file=allfc.file1)
      Sys.chmod(allfc.file1, "0666")
      remove(AA)
    }

    ## update sigdb or create if not exists
    sigdb <- file.path(pgx.dir, "datasets-sigdb.h5")  
    if(!file.exists(sigdb) || pgxfc.changed) {
       ## NEED RETHINK!!!! HERE???
       if(file.exists(sigdb)) unlink(sigdb)
       if(verbose) message("[initDatasetFolder] creating signature DB to",sigdb,"...")       
       pgx.createSignatureDatabaseH5.fromMatrix(sigdb, X=allFC)
       if(verbose) message("[initDatasetFolder] add enrichment signature to",sigdb,"...")       
       pgx.addEnrichmentSignaturesH5(sigdb, X=allFC, methods = "rankcor")
    }

    pgxinfo.file <- file.path(pgx.dir, file)
    if(!file.exists(pgxinfo.file)) return(NULL)  ## no info??
    ## do not use fread.csv or fread here!! see issue #441
    pgxinfo = read.csv(pgxinfo.file, stringsAsFactors=FALSE, row.names=1, sep=',')
    pgxinfo$path <- pgx.dir
    return(pgxinfo) 
}


#' Update PGX-table with new pgx object. 
#'
#' @export
pgx.updateInfoPGX <- function(pgxinfo, pgx, remove.old=TRUE)
{

    cond = grep("title|source|group|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",
                colnames(pgx$samples),invert=TRUE,value=TRUE)
    cond = grep("title|source|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",
                colnames(pgx$samples),invert=TRUE,value=TRUE)
    cond

    is.mouse = (mean(grepl("[a-z]",pgx$genes$gene_name))>0.8)
    organism = c("human","mouse")[1 + is.mouse]
    if("organism" %in% names(pgx)) organism <- pgx$organism

    this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    date = ifelse(is.null(pgx$date), this.date, as.character(pgx$date))
    dataset.name <- pgx$name
    ## dataset.name <- ifelse(is.null(pgx$name), pgxfile, pgx$name)
    creator <- ifelse("creator" %in% names(pgx), pgx$creator, "")

    this.info <- c(
        dataset = dataset.name,
        ## author = "", ## add author? maintainer? owner??
        creator = creator,
        ## collection = pgx$collection,
        datatype = ifelse(is.null(pgx$datatype),"", pgx$datatype),
        description = ifelse(is.null(pgx$description),"", pgx$description),
        organism = organism,
        nsamples = nrow(pgx$samples),
        ngenes = nrow(pgx$X),
        nsets = nrow(pgx$gsetX),
        conditions = paste(cond,collapse=" "),
        date = as.character(date),
        path = NULL
    )

    ## force to be character...
    if(!is.null(pgxinfo) && NCOL(pgxinfo)>0 && nrow(pgxinfo)>0 )
    {
        if("date" %in% colnames(pgxinfo)) {
            pgxinfo$date <- as.character(pgxinfo$date)
        }
        if(!"creator" %in% colnames(pgxinfo)) {
            pgxinfo$creator <- ""
        }
        which.factor <- which(sapply(pgxinfo,is.factor))
        which.factor
        for(i in which.factor) {
            pgxinfo[,i] <- as.character(pgxinfo[,i])
        }

        ## remove existing entries??
        if(remove.old && nrow(pgxinfo)>0 ) {
            d1 <- sub("[.]pgx$","",pgxinfo$dataset)
            d2 <- sub("[.]pgx$","",this.info["dataset"])
            if(!is.null(d2) && !is.na(d2) && d2 %in% d1 && d2!="") {
              sel <- which(d1!=d2)
              pgxinfo <- pgxinfo[sel,,drop=FALSE]
            }
        }

        ## merge with same columns
        info.cols <- colnames(pgxinfo)
        info.cols <- unique(c(info.cols, names(this.info)))
        this.info = this.info[match(info.cols,names(this.info))]
        names(this.info) = info.cols
        pgxinfo1 <- pgxinfo
        for(f in setdiff(info.cols,colnames(pgxinfo1))) {
            pgxinfo1[[f]] <- NA
        }
        match(info.cols,colnames(pgxinfo1))
        pgxinfo1 = pgxinfo1[,match(info.cols,colnames(pgxinfo1)),drop=FALSE]
        colnames(pgxinfo1) = info.cols
        pgxinfo <- rbind( pgxinfo1, this.info)
    } else {
        pgxinfo <- data.frame(rbind(this.info))
    }

    pgxinfo
}

##================================================================================
## PGXINFO methods (preferred API)
##================================================================================

if(0) {
  pgx.dir="~/Playground/omicsplayground/data"
  file="datasets-info.csv"
  pgxname="example-data"
}

#' Delete pgx entry in datasets-info table in files
#'
#' @export
pgxinfo.deletePgx <- function(pgx.dir, pgxname, file="datasets-info.csv",
                              delete.fc = FALSE)
{
    info.file <- file.path(pgx.dir, file)
    pgxname <- sub("[.]pgx$","",pgxname)
    pgxinfo <- read.csv(info.file, row.names=1)
    pgxinfo
    
    info_datasets <- sub("[.]pgx$","",pgxinfo$dataset)
    sel <- which(info_datasets == pgxname)
    sel
    if(!is.null(sel) && length(sel)>0) {
      pgxinfo <- pgxinfo[-sel,,drop=FALSE]
      write.csv(pgxinfo, file=info.file)
    }

    ## Should we also delete the entry in allFC and sigdb? This will
    ## take some overhead, so if not needed better skip.
    allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
    tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
    h5.file <- file.path(pgx.dir, "datasets-sigdb.h5")
    
    ## delete columns from allFC file
    if(delete.fc && file.exists(allfc.file)) {
      allFC <- data.table::fread(allfc.file,check.names=FALSE) ## HEADER!!!        
      allFC <- as.matrix(allFC, rownames=1)
      del <- grep(paste0("\\[",pgxname,"\\]"), colnames(allFC))
      del
      if(length(del)>0) {
        allFC <- allFC[,-del,drop=FALSE]
        allFC1 <- data.frame(gene=rownames(allFC),allFC,check.names=FALSE)
        data.table::fwrite(allFC1, allfc.file) ## HEADER!!!        
      }
    }

    ## delete rows from t-SNE file
    if(delete.fc && file.exists(tsne.file)) {
      tsne <- data.table::fread(tsne.file,check.names=FALSE) ## HEADER!!!        
      tsne <- as.matrix(tsne, rownames=1)
      del <- grep(paste0("\\[",pgxname,"\\]"), rownames(tsne))
      if(length(del)) {
        tsne <- tsne[-del,,drop=FALSE]
        write.csv(tsne, file=tsne.file)
      }
    }
    
    ## delete dataset from H5 file
    if(delete.fc && file.exists(h5.file)) {
      ##rhdf5::h5ls(h5.file)
      sigdb.removeDataset(h5.file, pgxname)
      ##rhdf5::h5ls(h5.file)
    }
   
    ## return(pgxinfo)
}


#' Update PGX-info file with new pgx object
#'
#' @export
pgxinfo.addPgx <- function(pgx.dir, pgx, file="datasets-info.csv",
                           update.fc=TRUE )
{
    info.file <- file.path(pgx.dir, file)
    pgxinfo <- read.csv(info.file, row.names=1)
    pgxinfo <- pgx.updateInfoPGX(pgxinfo, pgx, remove.old=TRUE)
    write.csv(pgxinfo, file=info.file)

    allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
    tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
    h5.file <- file.path(pgx.dir, "datasets-sigdb.h5")
    pgxname <- sub("[.]pgx$","",pgx$name)
    pgxname

    ## add FC columns to allFC file
    if(update.fc && !file.exists(allfc.file)) {
        ## complete update of all files
        pgx.initDatasetFolder(pgx.dir)
        return()
    }
    
    ## add FC columns to allFC file
    allfC <- NULL
    if(update.fc && file.exists(allfc.file)) {
      allFC <- data.table::fread(allfc.file,check.names=FALSE) ## HEADER!!!        
      allFC <- as.matrix(allFC, rownames=1)
      del <- grep(paste0("\\[",pgxname,"\\]"), colnames(allFC))
      if(length(del)>0) {
        allFC <- allFC[,-del,drop=FALSE]
      }
      F1 <- playbase::pgx.getMetaMatrix(pgx)$fc
      F1 <- F1[match(rownames(allFC),rownames(F1)),]
      colnames(F1) <- paste0("[",pgxname,"] ",colnames(F1))
      allFC <-cbind(allFC, F1)
      allFC1 <- data.frame(gene=rownames(allFC),allFC, check.names=FALSE)
      data.table::fwrite(allFC1, allfc.file) ## HEADER!!!        

      ## NEED RETHINK!!!! This could be done perhaps more efficient
      ## when updating with one extra dataset.
      pgx.createSignatureDatabaseH5.fromMatrix(h5.file, X=allFC,
        update.only=FALSE)
      
      ## update t-SNE file (in the future we do not need this file)
      ## Just get the tSNE from sigdb
      cn   <- rhdf5::h5read(h5.file, "data/colnames")
      tsne <- rhdf5::h5read(h5.file, "clustering/tsne2d")
      rownames(tsne) <- cn
      colnames(tsne) <- paste0("tsne.",1:ncol(tsne))
      write.csv(tsne, file=tsne.file)
    } 
}

#' @export
pgxinfo.read <- function(pgx.dir, file="datasets-info.csv", match=TRUE)
{
  pgx.files <- dir(pgx.dir, pattern="[.]pgx$")
  if(length(pgx.files)==0) return(NULL)  ## no files!

  pgxinfo.file <- file.path(pgx.dir, file)
  if(!file.exists(pgxinfo.file)) return(NULL)  ## no info??
  ## do not use fread.csv or fread here!! see issue #441
  pgxinfo = read.csv(pgxinfo.file, stringsAsFactors=FALSE, row.names=1, sep=',')
  if(match) {
    pgx.files1 <- sub("[.]pgx$","",pgx.files)
    pgxinfo.datasets <- sub("[.]pgx$","",pgxinfo$dataset)
    sel <- pgxinfo.datasets %in% pgx.files1
    pgxinfo <- pgxinfo[sel,,drop=FALSE]
  }
  return(pgxinfo)
}
