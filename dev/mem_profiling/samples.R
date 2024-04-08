library(playbase)
library(profmem)


options(profmem.threshold = 20000)
n = dim(samples)[1]


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

    INPUTS_CHECKED <- pgx.crosscheckINPUT(samples, counts, contrasts)

    return(list(samples=INPUTS_CHECKED$SAMPLES, counts=INPUTS_CHECKED$COUNTS, contrasts=INPUTS_CHECKED$CONTRASTS))
}

input <- duplicate_samples(n=18)

pgx <- playbase::pgx.createPGX(
 samples = input$samples,
 counts = input$counts,
 contrasts = input$contrasts
)

# Step 2. Populate pgx object with results


p <- profmem({
    pgx <- playbase::pgx.computePGX(
        pgx = pgx
        )
    })

