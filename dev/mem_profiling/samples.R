library(playbase)

options(app.profile = TRUE)

print(getwd())

output_dir <- getwd() # root playbase
output_dir <- paste0(output_dir, "/dev/mem_profiling")

# read params of a FULL settings computation
params <- readRDS(file.path(output_dir,"params.RData"))

duplicate_samples <- function(n) {
    samples <- playbase::SAMPLES
    counts <- playbase::COUNTS
    contrasts <- playbase::CONTRASTS

    # check that n is a multiple of samples
    if (n %% dim(samples)[1] != 0) {
        stop("n must be a multiple of the number of samples")
    }

    # create a hash with same length as n
    hash <- paste0("sample", 1:n)

    # duplicate samples
    duplications <- rep(1:dim(samples)[1], n/dim(samples)[1])
    samples <- samples[duplications, ]

    # add hash to rownames samples
    rownames(samples) <- paste(rownames(samples), hash, sep="_")

    # duplicate counts
    counts <- counts[,duplications]

    # add hash to rownames counts
    colnames(counts) <- paste(colnames(counts), hash, sep="_")

    # add random noise in counts
    noise_factor <- 0.05
    counts <- counts + abs(rnorm(n = length(counts))* noise_factor)

    INPUTS_CHECKED <- pgx.crosscheckINPUT(samples, counts, contrasts)

    return(list(samples=INPUTS_CHECKED$SAMPLES, counts=INPUTS_CHECKED$COUNTS, contrasts=INPUTS_CHECKED$CONTRASTS))
}


iterations = seq(1,101,10) * 18

for(i in iterations) {
    input <- duplicate_samples(n=i)
    
    pgx <- playbase::pgx.createPGX(
        organism = params$organism,
        counts = input$counts,
        samples = input$samples,
        contrasts = input$contrasts,
        name = params$name,
        datatype = params$datatype,
        description = params$description,
        creator = params$creator,
        X = NULL,
        batch.correct = params$batch.correct,
        normalize = params$normalize,
        prune.samples = params$prune.samples,
        filter.genes = params$filter.genes,
        only.known = params$only.known,
        only.proteincoding = params$only.proteincoding,
        only.hugo = params$only.hugo,
        convert.hugo = params$convert.hugo,
        custom.geneset = params$custom.geneset,
        max.genesets = params$max.genesets,
        annot_table = params$annot_table
    )
    
    write.csv(summaryRprof(memory = "tseries"), paste(output_dir,"/samples_output/","sample_mem_create_tseries_", i, ".csv", sep=""))
    
    pgx <- playbase::pgx.computePGX(
        pgx = pgx,
        max.genes = params$max.genes,
        gx.methods = params$gx.methods,
        gset.methods = params$gset.methods,
        extra.methods = params$extra.methods,
        use.design = params$use.design,        ## no.design+prune are combined
        prune.samples = params$prune.samples,  ##
        do.clustergenes = params$do.cluster,
        do.clustergenesets = params$do.cluster,
        cluster.contrasts = params$cluster.contrasts,
        pgx.dir = params$pgx.save.folder,
        libx.dir = "./libx",
        user_input_dir = output_dir
    )

    write.csv(summaryRprof(memory = "tseries"), paste(output_dir,"/samples_output/","sample_mem_compute_tseries_", i, ".csv", sep=""))
}
