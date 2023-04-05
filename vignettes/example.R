
devtools::load_all()
## example files can be accessed via the playbase::example_file() function
counts <- playbase::read_counts(playbase::example_file("counts.csv"))
samples <- playbase::read_samples(playbase::example_file("samples.csv"))
contrasts <- playbase::read_contrasts(playbase::example_file("contrasts.csv"))

## create a pgx object (new style)
pgx <- playbase::pgx.createPGX(counts, samples, contrasts)

## old style still works...
pgx <- playbase::pgx.createPGX(counts, samples, contrasts)
