##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##
##


##======================================================================
##==================== TORCH ===========================================
##======================================================================

#'
#'
#' @export
MultiOmicsSAE_torch <- R6::R6Class(
  "MultiOmicsSAE.torch",
  private = list(),
  public = list(
    X = NULL,
    y = NULL,
    nviews = NULL,
    labels = NULL,        
    model = NULL,
    x_train = NULL,
    y_train = NULL,        
    sdx = NULL,
    loss_weights = NULL,
    
    initialize = function(X, y, ae_dims=c(100,20), dropout=0.0,
                          loss_weights=c(y=1, ae=10, l1=0.1, l2=0.01)
                          ) {
      self$X <- X
      self$nviews <- length(X)
      self$sdx <- lapply(X, matrixStats::rowSds)
      self$x_train <- lapply( lapply(X,t), scale)
      self$x_train <- lapply(self$x_train, torch::torch_tensor)

      self$y <- factor(y)
      self$labels <- levels(self$y)
      self$y_train <- as.integer(factor(self$y))
      self$y_train <- torch::torch_tensor(self$y_train)
      ydim <- length(unique(self$labels))

      self$model <- torch.multiOmicsSAE_module(
        self$x_train, ydim, ae_dims=ae_dims, dropout = dropout)

      if(length(loss_weights)!=4) stop("loss_weights wrong length (must be 4)")
      self$loss_weights <- loss_weights
    },
    fit = function(niter=200) {      
      ##niter=200
      model <- self$model
      optimizer <- torch::optim_adam( model$parameters, lr = 0.01 )
      xx <- self$x_train
      y <- self$y_train
      
      for (t in 1:niter) {
        output <- model(xx)
        enc_weights <-  lapply(model$encoder, function(m) m$parameters[[1]])
        y_loss  <- torch::nnf_cross_entropy( output$y_pred, y)
        ae_loss <- Reduce('+',mapply( torch::nnf_mse_loss, output$ae, xx))
        l1_loss <- Reduce('+',sapply(enc_weights, function(w) torch::torch_norm(w,1)))
        l2_loss <- Reduce('+',sapply(enc_weights, function(w) torch::torch_norm(w,2)))
        ##loss <- 1*y_loss + 1e1*ae_loss + 0.1*l1_loss + 0.01*l2_loss
        w <- self$loss_weights
        loss <- w['y']*y_loss + w['ae']*ae_loss + w['l1']*l1_loss + w['l2']*l2_loss         
        if (t %% 50 == 0) {
          cat("Y loss: ", as.numeric(y_loss ), "\n")
          cat("AE loss: ", as.numeric(ae_loss), "\n")
        }
        optimizer$zero_grad()
        loss$backward()
        optimizer$step()
      }

    },
    predict = function(X = NULL) {
      if(!is.null(X)) {
        xx <- lapply(X, function(x) scale(t(x)))
        xx <- lapply(xx, torch::torch_tensor)        
      } else {
        X <- self$X
        xx <- self$x_train
      }
      output <- self$model(xx)
      y_pred <- as.matrix(output$y_pred)
      rownames(y_pred) <- colnames(X[[1]])
      colnames(y_pred) <- self$labels
      y_pred
    },
    
    get_redux = function(xx = NULL, tsne=FALSE) {

      if(is.null(xx)) {
        xx <- self$x_train
      } else {
        xx <- lapply(xx, function(x) scale(t(x)))
        xx <- lapply(xx, torch::torch_tensor)        
      }
      output <- self$model(xx)
      redux <- c( output$enc, "multi-omics" = output$combined )
      redux <- lapply( redux, as.matrix )
      if(tsne) {
        px <- min(30, nrow(redux[[1]])/4)
        redux <- lapply( redux, function(r)
          Rtsne::Rtsne(r, perplexity=px, check_duplicates=FALSE)$Y )
      }
      redux
    },
    get_gradients = function(sd.weight=TRUE) {

      sdx <- self$sdx
      x_train <- self$x_train
      
      ## calculate feature perturbation response 
      grad <- list()
      k=1
      for(k in names(x_train)) {
        x <- x_train[[k]]
        dim(x)  
        delta <- rbind(0, diag(ncol(x)))
        xx <- lapply(x_train, function(x) matrix(0, nrow(delta), ncol(x)))
        xx[[k]] <- delta
        xx <- lapply(xx, torch::torch_tensor)        
        output <- self$model(xx)
        grad_k <- as.matrix(output$y_pred)
        colnames(grad_k) <- self$labels
        rownames(grad_k) <- c("ref", rownames(self$X[[k]]))
        
        ## convert perturbation response to gradient
        ##z <- log( grad_k / (1 - grad_k) )
        z <- grad_k
        dz <- t(t(z) - z["ref",])
        dz <- dz[-1,,drop=FALSE] ## remove ref
        if(sd.weight) {
          ##dz <- dz[names(sdx[[k]]),] * sdx[[k]]  ## SD is really important!!
          dz <- dz * sdx[[k]][rownames(dz)]  ## SD is really important!!
        }
        grad[[k]] <- dz
      }
      grad      
    }
  )
)

torch.multiOmicsSAE_module <- torch::nn_module(
  "multiomicsSAE",
  initialize = function(xx, ydim, ae_dims=c(100,20), dropout=0.0) {
    self$encoder <- torch::nn_module_list()
    self$decoder <- torch::nn_module_list()
    for(k in names(xx)) {
      self$encoder[[k]] <- torch::nn_sequential(
        torch::nn_linear(ncol(xx[[k]]), ae_dims[1]),
        torch::nn_relu(),
        torch::nn_dropout(dropout),        
        torch::nn_linear(ae_dims[1], ae_dims[2]),        
        torch::nn_relu(),
        torch::nn_dropout(dropout)        
      )
    }
    for(k in names(xx)) {
      self$decoder[[k]] <- torch::nn_sequential(
        torch::nn_linear(ae_dims[2], ae_dims[1]),
        torch::nn_relu(),
        torch::nn_dropout(dropout)                ,
        torch::nn_linear(ae_dims[1], ncol(xx[[k]]))
      )
    }
    nviews <- length(xx)    
    nmerge <- nviews*ae_dims[2]
    mlp_dim <- round(nmerge/2)
    self$integrator <- torch::nn_sequential(
      torch::nn_linear(nmerge, mlp_dim),
      torch::nn_relu(),
      torch::nn_dropout(dropout)              
    )
    self$predictor <- torch::nn_sequential(
      torch::nn_linear(mlp_dim, ydim)
    )
  },
  forward = function(xx) {
    enc_output <- list()
    ae_output <- list()
    for(k in names(xx)) {
      enc_output[[k]] <- self$encoder[[k]](xx[[k]])
      ae_output[[k]]  <- self$decoder[[k]](enc_output[[k]])
    }
    combined_output <- self$integrator( torch::torch_cat(enc_output,dim=2))
    pred_output <- self$predictor(combined_output)
    list( y_pred = pred_output, enc = enc_output,
         ae = ae_output, combined = combined_output)
  }
)


#'
#'
#' @export
deep.getConfusionMatrix <- function(net) {
  output <- net$model( net$x_train )
  y_pred <- net$labels[ max.col(output$y_pred) ]
  mat <- table( actual=net$y, y_pred, useNA="always")
  mat <- mat[match(net$labels,rownames(mat)),
             match(net$labels,colnames(mat))]
  rownames(mat) <- colnames(mat) <- net$labels
  mat[is.na(mat)] <- 0
  mat
}

#'
#'
#' @export
deep.plotAutoEncoderReconstruction <- function(net, dtypes=NULL,
                                               main=NULL, par=TRUE) {

  output <- net$model( net$x_train )
  str(output)
  if(is.null(dtypes)) {
    dtypes <- names(net$X)
  } else {
    intersect(dtypes, c("mixed",names(net$X)))
  }
  if(length(dtypes)==0) return(NULL)
  if(par) {
    nc <- ceiling(sqrt(length(dtypes)))
    nr <- ceiling(length(dtypes)/nc)
    par(mfrow=c(nr,nc))
  }
  if("mixed" %in% dtypes) {
    if(is.null(main)) main <- "mixed features"
    dtypes <- names(net$X)
    sX <- t(do.call(rbind, net$X[dtypes]))
    pX <- lapply(output$ae[dtypes], as.matrix)
    pX <- do.call(cbind, pX)
    sX <- scale(sX)
    jj <- sample(length(sX),1000,replace=TRUE)
    plot( sX[jj], pX[jj], cex=0.5, main=main,
         xlab="actual", ylab="reconstructed")
  } else {
    for(k in dtypes) {
      sX <- scale(t(net$X[[k]]))
      pX <- as.matrix(output$ae[[k]])
      jj <- sample(length(sX),1000,replace=TRUE)
      tt <- k
      if(!is.null(main)) tt <- main
      plot( sX[jj], pX[jj], cex=0.5, main=tt,
           xlab="actual", ylab="reconstructed")
    }
  }

}


#'
#'
#' @export
deep.plotRedux <- function(net, views=NULL, par=TRUE, cex=1) {
  redux <- net$get_redux(xx=NULL, tsne=TRUE)
  y <- net$y
  if(!is.null(views)) redux <- redux[names(redux) %in% views]
  if(par) {
    nq <- ceiling(sqrt(length(redux)))
    par(mfrow=c(nq,nq), mar=c(4,4,4,4))
  }
  for(k in names(redux)) {
    plot(redux[[k]], col=factor(y), main=k,
         cex=cex, pch=19, xlab="tSNE-x", ylab="tSNE-y")
  }
}

#'
#'
#' @export
deep.plotMultiOmicsGradients <- function(grad, n=20, cex.names=1, par=TRUE) {
  if(par) {
    ngroup <- ncol(grad[[1]])
    nview <- length(grad)  
    par(mfrow=c(ngroup,nview), mar=c(10,4,0.5,2))
  }
  for(k in 1:ncol(grad[[1]])) {
    for(i in 1:length(grad)) {
      gr <- grad[[i]][,k]
      gr[is.na(gr) | is.infinite(gr)] <- 0
      gr <- head(gr[order(-abs(gr))], n)
      barplot(sort(gr), ylab="gradient", las=3, cex.names=cex.names)
      title( paste0("phenotype=", colnames(grad[[1]])[k],
                    "; datatype=", names(grad)[i]),
            line=-0.5, cex.main=1)
    }
  }
}

#'
#'
#' @export
deep.plotGradientVSFoldchange <- function(grad, X, y, fc=NULL, p.weight=0,
                                          data=FALSE, par=TRUE) {

  ## Calculate correspoding T-test statistics
  if(is.null(fc)) {
    fc <- list()
    pv <- list()
    pheno <- colnames(grad[[1]])
    ct <- sapply( pheno, function(p) ifelse( y==p, as.character(y), "others"))
    for(k in names(grad)) {
      res <- lapply( colnames(ct), function(j) {
        gx.limma(X[[k]], ct[,j], ref='others', lfc=0, fdr=1, sort.by='none')
      })
      fc[[k]] <- sapply(res, function(r) r$logFC)
      ##pv[[k]] <- sapply(res, function(r) r$P.Value)
      pv[[k]] <- sapply(res, function(r) r$adj.P.Val)            
      rownames(fc[[k]]) <- rownames(pv[[k]]) <- rownames(res[[1]])
      colnames(fc[[k]]) <- colnames(pv[[k]]) <- colnames(ct)
    }
  }

  ## align datatypes
  dt <- intersect(names(grad),names(fc))
  grad <- grad[dt]
  fc <- fc[dt]
  pv <- pv[dt]  
  
  ## align features
  for(k in names(grad)) {
    pp <- intersect( rownames(grad[[k]]), rownames(fc[[k]]) )
    ph <- intersect( colnames(grad[[k]]), colnames(fc[[k]]) )
    grad[[k]] <- grad[[k]][pp,ph,drop=FALSE]
    fc[[k]] <- fc[[k]][pp,ph,drop=FALSE]
    pv[[k]] <- pv[[k]][pp,ph,drop=FALSE]        
  }
  
  if(data) {
    G <- do.call( rbind, mofa.prefix(grad))
    F <- do.call( rbind, mofa.prefix(fc))
    P <- do.call( rbind, mofa.prefix(pv))        
    colnames(G) <- paste0("grad.",colnames(G))
    colnames(F) <- paste0("logFC.",colnames(F))
    colnames(P) <- paste0("pval.",colnames(P))        

    df <- data.frame(G, F)
    avg.rnk <- rank(rowMeans(G**2)) + rank(rowMeans(F**2))
    df <- df[order(-avg.rnk),,drop=FALSE]
    return(df)
  }
    
  ## plot gradient vs foldchagne 
  if(par) {
    ngrad <- length(grad)*ncol(grad[[1]])
    nc <- ceiling(sqrt(ngrad))
    nr <- ceiling(ngrad/nc)
    par(mfrow=c(nr,nc), mar=c(4.5,4,3,2))
  }
  i=2
  k=3
  for(i in 1:ncol(grad[[1]])) {
    for(k in names(grad)) {
      f <- fc[[k]][,i]
      p <- pv[[k]][,i]      
      g <- grad[[k]][,i]
      xlab = "logFC"
      if(p.weight==1) {
        f <- f * -log10(p)
        xlab = "-log10(p) * logFC"
      }
      if(p.weight==2) {
        f <- f * (1-p)**4
        xlab = "(1-p) * logFC"
      }
      plot( f, g, xlab = xlab, ylab = "network gradient")
      abline(h=0, v=0, lty=2)
      title( paste0("phenotype=", colnames(grad[[k]])[i],
                    "; datatype=", k))
      sel <- unique(c( head(order(-f),5), head(order(-g),5)))
      text( f[sel], g[sel],  labels = names(f)[sel],
           col="red", pos=2) 
    }
  }
}



#'
#'
#' @export
deep.plotNeuralNet <- function(net, ae_dims=c(100,20), outfile=NULL) {

  if(!dir.exists("/opt/PlotNeuralNet")) {
    message("ERROR: please install PlotNeuralNet in /opt")
    return(NULL)
  }
  
  if(!is.null(outfile)) {
    if(!grepl("svg$",outfile)) {
      message("ERROR: output file must be SVG")
      return(NULL)
    }
  }
  
  header_src = glue::glue("
import sys
sys.path.append('../')
from pycore.tikzeng import *

# defined your arch
arch = [
    to_head( '..' ),
    to_cor(),
    to_begin(),
", .trim=FALSE)
  
  str_at <- function(at) paste0('(',at[1],',',at[2],',',at[3],')')
  input_src <- function(name, at, image, caption) {
    to1 <- str_at(c(at[1]-2.5,at[2],at[3]))
    to2 <- str_at(at)
    glue::glue("
    # input {name}
    to_input( '../images/heatmap.png', to='{to1}', width=9.5, height=2.2 ),
    to_Conv('{name}', '', '', offset='(-4,0,6.5)', to='{to2}', height=0, depth=0, width=0, caption='{caption}'),
", .trim=FALSE)
  }

  autoencoder_src <- function(name, caption, at, dims) {

    to <- str_at(at)
    enc1 <- paste0("enc",name,"1")
    enc2 <- paste0("enc",name,"2")
    btn  <- paste0("btn",name)
    dec1 <- paste0("dec",name,"1")
    dec2 <- paste0("dec",name,"2")    
      
    glue::glue("
    # Auto-encoder {name}
    to_Conv('{enc1}', '{dims[1]}', '', offset='(0,0,0)', to='{to}', height=10, depth=45, width=2),
    to_Conv('{enc2}', '{dims[2]}', '', offset='(2,0,0)', to='({enc1}-east)', height=10, depth=20, width=2, caption='{caption}'),
    to_Pool('{btn}', offset='(2,0,0)', to='({enc2}-east)', height=10, depth=10, width=2),        
    to_Conv('{dec1}', '{dims[2]}', '', offset='(2,0,0)', to='({btn}-east)', height=10, depth=20, width=2),
    to_Conv('{dec2}', '{dims[1]}', '', offset='(2.5,0,0)', to='({dec1}-east)', height=10, depth=45, width=2),
    to_connection('{enc1}', '{enc2}'),
    to_connection('{enc2}', '{btn}'),
    to_connection('{btn}', '{dec1}'),
    to_connection('{dec1}', '{dec2}'),
    to_skip( of='{btn}', to='sum1', pos=1.4),
", .trim = FALSE)

  }
  sum_src = function() {
    nview <- length(views)
    to <- str_at( c(0,0, (nview-1)*22/2))
    glue::glue("    
    # Combine layer
    to_Sum('sum1', offset='(17,0,0)', to='{to}', radius=2.5, opacity=0.6),
    to_Conv('label2', '', '', offset='(16.2,0,0)', to='{to}', height=0, depth=0, width=0, caption='MERGE'),
", .trim = FALSE)
  }
  dense_src = function(ydim, dims=c(100,20)) {
    glue::glue("        
    # Dense prediction layer (2 layers)
    to_Conv('dense1', '{dims[1]}', '', offset='(4,0,0)', to='(sum1-east)', height=10, depth=35, width=2, caption='classifier'),
    to_Conv('dense2', '{dims[2]}', '', offset='(2,0,0)', to='(dense1-east)', height=10, depth=20, width=2),
    to_ConvSoftMax( name='soft1', s_filer = '{ydim}', offset='(2,0,0)', to='(dense2-east)', width=2, height=10, depth=10, caption='PREDICTION'),
    to_connection( 'sum1',  'dense1'),
    to_connection( 'dense1', 'dense2'),
    to_connection( 'dense2', 'soft1'),    
", .trim = FALSE)
  }
  end_src = function(nview) {
    zoff <- (nview-1)*22
    glue::glue("
    ## margins
    to_Conv('empty1', '', '', offset='(0,0,0)', to='(-7,0,{zoff})', height=0, depth=0, width=0),
    to_Conv('empty2', '', '', offset='(5,0,0)', to='(dense2-east)', height=0, depth=0, width=0),
    to_end()
    ]

def main():
    namefile = str(sys.argv[0]).split('.')[0]
    to_generate(arch, namefile + '.tex' )

if __name__ == '__main__':
    main()
", .trim=FALSE)
  }

  
  ## build code text
  ltx <- header_src
  views <- rev(names(net$X))
  nview <- length(views)

  redux <- net$get_redux()
  rdim <- ncol(redux[[1]]) ## bottleneck dimension

  for(i in 1:length(views)) {
    at <- c(0, 0, (i-1)*22 )
    caption <- toupper(views[i])
    if(max(nchar(views))<10) caption <- paste("datatype:~",caption)
    ltx1 <- input_src( paste0("input",views[i]), at, "", caption=caption )
    ltx <- paste(ltx, ltx1)
  }
  
  ltx <- paste(ltx, sum_src())
  for(i in 1:length(views)) {
    at <- c(0, 0, (i-1)*22 )
    xdim <- nrow(net$X[[views[i]]])
    dims <- c(xdim, ae_dims)
    caption <- paste0("autoencoder~",toupper(views[i]))
    ltx1 <- autoencoder_src(views[i], at=at, caption=caption, dims=dims )
    ltx <- paste(ltx, ltx1)
  }

  mdim <- nview*rdim
  ltx <- paste(ltx, dense_src(ydim=10, dims=c(mdim,mdim/2)))  
  ltx <- paste(ltx, end_src(nview=nview))  
  
  ## create PDF using PlotNeuralNet
  pyfile <- tempfile(fileext=".py", tmpdir="/opt/PlotNeuralNet/pyexamples")
  ##pyfile="/opt/PlotNeuralNet/pyexamples/model2.py"
  write(ltx, file=pyfile)
  texfile <- sub("[.]py$",".tex",pyfile)
  auxfile <- sub("[.]py$",".aux",pyfile)
  logfile <- sub("[.]py$",".log",pyfile)  
  pdffile <- sub("[.]py$",".pdf",pyfile)
  svgfile <- sub("[.]py$",".svg",pyfile)    
  cmd <- glue::glue("cd /opt/PlotNeuralNet/pyexamples/ && python {pyfile} && pdflatex {texfile} && pdf2svg {pdffile} {svgfile}")
  suppressMessages( system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE, intern=FALSE) )
  
  if(file.exists(auxfile)) unlink(auxfile)
  if(file.exists(logfile)) unlink(logfile)
  if(file.exists(texfile)) unlink(texfile)
  if(file.exists(pyfile))  unlink(pyfile)
  if(file.exists(pdffile)) unlink(pdffile)    
  
  if(file.exists(svgfile)) {
    if(is.null(outfile)) {
      outfile <- tempfile(fileext=".svg")
    }
    file.rename(svgfile, outfile)
    return(outfile)
  } else {
    return(NULL)
  }
}

##======================================================================
##==================== KERAS ===========================================
##======================================================================

#'
#'
#' @export
MultiOmicsSAE_keras <- R6::R6Class(
  "MultiOmicsSAE.keras",
  private = list(),
  public = list(
    X = NULL,
    y = NULL,
    nviews = NULL,
    output_labels = NULL,        
    pred_model = NULL,
    full_model = NULL,    
    ae_model = NULL,    
    dr_model = NULL,
    bottleneck_model = NULL,    
    history = NULL,
    x_train = NULL,
    y_train = NULL,        
    sdx = NULL,
    
    initialize = function(X, y, l1=0, l2=0) {
      self$X <- X
      self$y <- factor(y)
      self$nviews <- length(X)
      self$output_labels <- levels(y)
      self$create_model(l1=l1, l2=l2)
    },

    create_model = function(l1=0.01, l2=0.01) {

      X <- self$X
      y <- self$y      
      nviews <- self$nviews

      ## prescale variable, save original SD
      self$sdx <- lapply(X, function(x) matrixStats::rowSds(x))
      x_train <- lapply(X, function(x) scale(t(x)))
      ##y_train <- keras::to_categorical(factor(y))
      y_train <- matrix(0, nrow(x_train[[1]]),length(unique(y)))
      y_train[cbind(1:nrow(y_train),as.integer(factor(y)))] <- 1
      
      ## Create inputs
      create_inputs <- function(x) {
        layer_input(shape = ncol(x))
      }
      omics_inputs <- lapply(x_train, create_inputs)
      
      ## Create dense layers for each view
      create_dense <- function(id, x) {
        layer_dense(
          x, units=64, activation="relu",
          kernel_regularizer = regularizer_l1_l2(l1, l2)
        ) %>%
          layer_dropout(0.1) %>%
          layer_dense(units=16, activation="relu", name=id)        
      }      
      omics_dense <- list()
      i=1
      bottlenecks <- c()
      for(i in 1:length(omics_inputs)) {
        bottlenecks[i] <- paste0("bottleneck_",names(omics_inputs)[i])
        omics_dense[[i]] <- create_dense(
          id = bottlenecks[i], omics_inputs[[i]]
        )
      }
      bottlenecks
      
      ## Merge latent vectors
      combined <- layer_concatenate( unname(omics_dense) )
      
      ## integrated classifier
      output <- combined %>%
        layer_dense(24, name="DR1") %>%
        layer_dropout(0.1) %>%        
        layer_dense(ncol(y_train), activation = 'softmax')
      
      ## AE outputs
      create_ae <- function(x, dim) {
        layer_dense(x, units=64, activation="relu") %>%
          layer_dropout(0.1) %>%        
          layer_dense(units = dim)  
      }
      input.sizes <- sapply(omics_inputs,function(x) x$shape[[2]])
      omics_ae <- lapply(1:nviews, function(i)
        create_ae(omics_dense[[i]], input.sizes[i]))
      
      ae_model <- keras_model(
        inputs = unname(omics_inputs),
        outputs = c( omics_ae )
      )      

      pred_model <- keras_model(
        inputs = unname(omics_inputs),
        outputs = output
      )      

      full_model <- keras_model(
        inputs = unname(omics_inputs),
        outputs = c( output, omics_ae )
      )      

      dr_model <- keras_model(
        inputs = unname(omics_inputs),
        outputs = get_layer(pred_model, "DR1")$output
      )

      ## create bottleneck model
      bneck_layers <- lapply(bottlenecks, function(b)
        get_layer(pred_model, b)$output)
      bottleneck_model <- keras_model(
        inputs = unname(omics_inputs),
        outputs = bneck_layers 
      )
      
      self$x_train <- x_train
      self$y_train <- y_train      
      self$pred_model <- pred_model
      self$full_model <- full_model      
      self$ae_model <- ae_model      
      self$bottleneck_model <- bottleneck_model      
      self$dr_model <- dr_model    
    },
      
    fit_ae = function(epochs=20) {      

      x_train <- self$x_train
      y_train <- self$y_train
      nviews <- self$nviews
      
      # compile model
      self$ae_model %>% compile(
        loss = rep('mean_squared_error', nviews),    
        optimizer = "adam"
      )

      # fit model
      self$ae_model %>% fit(
        x = unname(x_train),
        y = unname(x_train),
        epochs = epochs,
        batch_size = 40,
        verbose = 1
      )
      ##self$history <- history
      
    },

    fit_pred = function( epochs=20 ) {      

      x_train <- self$x_train
      y_train <- self$y_train
      
      # compile model
      self$pred_model %>% compile(
        loss = 'categorical_crossentropy',
        optimizer = "adam"
      )

      # fit model
      history <- self$pred_model %>% fit(
        x = unname(x_train),
        y = y_train,
        epochs = epochs,
        batch_size = 40,
        verbose = 1
      )
      
      self$history <- history
    },
    
    fit_full = function( weights=c(1,1), epochs=20 ) {      
      
      x_train <- self$x_train
      y_train <- self$y_train
      nviews <- self$nviews
      
      # compile model
      self$full_model %>% compile(
        loss = c( 'categorical_crossentropy',
                 rep('mean_squared_error', nviews) ),
        loss_weights = c( weights[1], rep(weights[2], nviews)),  
        optimizer = "adam"
      )
      
      # fit model      
      history <- self$full_model %>% fit(
        x = unname(x_train),
        y = list( y_train, unname(x_train) ),
        epochs = epochs,
        batch_size = 40,
        verbose = 1
      )
      
      self$history <- history
    },
    
    get_redux = function(xx = NULL, tsne=FALSE) {
      if(is.null(xx)) xx <- self$X
      xx <- lapply(xx, function(x) scale(t(x)))      
      bn <- predict(self$bottleneck_model, unname(xx))
      if(!is.list(bn)) bn <- list(bn)
      names(bn) <- names(xx)
      dr <- list(predict(self$dr_model, unname(xx)))
      names(bn) <- paste0("bottleneck_",names(xx))
      names(dr) <- paste0("combined_dr")
      redux = c(bn, dr)
      if(tsne) {
        px <- min(30, nrow(redux[[1]])/4)
        redux <- lapply( redux, function(r)
          Rtsne::Rtsne(r, perplexity=px, check_duplicates=FALSE)$Y )
      }
      redux
    },
        
    predict = function(xx = NULL) {
      if(is.null(xx)) xx <- self$X
      xx <- lapply(xx, function(x) scale(t(x)))      
      pred <- predict(self$pred_model, unname(xx))
      rownames(pred) <- colnames(self$X[[1]])
      colnames(pred) <- self$output_labels
      pred
    },

    get_gradients = function(sd.weight=TRUE) {

      sdx <- self$sdx
      x_train <- self$x_train
      
      ## calculate feature perturbation response 
      grad <- list()
      k=1
      for(k in names(x_train)) {
        x <- x_train[[k]]
        dim(x)  
        delta <- rbind(0, diag(ncol(x)))
        xx <- lapply(x_train, function(x) matrix(0, nrow(delta), ncol(x)))
        xx[[k]] <- delta
        grad_k <- self$pred_model %>% predict(unname(xx))
        colnames(grad_k) <- self$output_labels
        rownames(grad_k) <- c("ref", colnames(x))
        
        ## convert perturbation response to gradient
        z <- log( grad_k / (1 - grad_k) )
        dz <- t(t(z) - z["ref",])
        dz <- dz[-1,,drop=FALSE] ## remove ref
        if(sd.weight) {
          ##dz <- dz[names(sdx[[k]]),] * sdx[[k]]  ## SD is really important!!
          dz <- dz * sdx[[k]][rownames(dz)]  ## SD is really important!!
        }
        grad[[k]] <- dz
      }
      grad
    }

  )
)
