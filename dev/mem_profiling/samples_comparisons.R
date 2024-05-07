library(playbase)

options(app.profile = TRUE)

print(getwd())

output_dir <- getwd() # root playbase
output_dir <- paste0(output_dir, "/dev/mem_profiling")

# read params of a FULL settings computation
params <- readRDS(file.path(output_dir,"params.RData"))

iterations <- as.integer(c(5,10,25,seq(1,10) *50))
iterations_samples <- seq(1,101,10) * 18
for(index in seq_along(iterations)) {
    
    i_contrasts=iterations[index]
    i_samples
    input <- playbase::duplicate_samples_contrasts(n=i_samples,c=i_contrasts)
    
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
    
    dir.create(file.path(output_dir,"contrasts_samples_output"), showWarnings = FALSE)
    
    write.csv(summaryRprof(memory = "tseries"), paste(output_dir,"/contrasts_samples_output/","sample_mem_create_tseries_", i_samples,"_",i_contrasts, ".csv", sep=""))
    
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
    
    write.csv(summaryRprof(memory = "tseries"), paste(output_dir,"/contrasts_samples_output/","sample_mem_compute_tseries_", i_samples,"_",i_contrasts,".csv", sep=""))
}
