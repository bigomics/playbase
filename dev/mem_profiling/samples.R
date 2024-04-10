library(playbase)

options(app.profile = TRUE)

output_dir <- getwd()

# save in dev/mem_profiling
output_dir <- paste0(output_dir, "/dev/mem_profiling/")

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
        samples = input$samples,
        counts = input$counts,
        contrasts = input$contrasts
    )
    write.csv(summaryRprof(memory = "tseries"), paste(output_dir,"sample_mem_create_tseries_", i, ".csv", sep=""))
    pgx <- playbase::pgx.computePGX(
        pgx = pgx
    )
    write.csv(summaryRprof(memory = "tseries"), paste(output_dir,"sample_mem_compute_tseries_", i, ".csv", sep=""))
}
