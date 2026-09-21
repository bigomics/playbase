
devtools::load_all()
## example files can be accessed via the playbase.ingest::example_file() function
counts <- playbase.ingest::read_counts(playbase.ingest::example_file("counts.csv"))
samples <- playbase.ingest::read_samples(playbase.ingest::example_file("samples.csv"))
contrasts <- playbase.ingest::read_contrasts(playbase.ingest::example_file("contrasts.csv"))

## create a pgx object
pgx <- playbase::pgx.createPGX(counts, samples, contrasts)

## compute a pgx object
pgx <- playbase::pgx.computePGX(pgx)
