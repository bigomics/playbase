library(playbase)
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
# real all samples files under samples_output
type <- c("samples", "contrasts","contrasts_gc","samples_contrasts")[3]

samples_output <- list.files(sprintf("dev/mem_profiling/%s_output_libx",type), full.names = TRUE)

samples_output <- samples_output

# read csv
samples <- lapply(samples_output, read.csv)

# add the name of the file (containing number of samples) as a column
samples <- lapply(seq_along(samples), function(i) {
    #i=1
    samples[[i]]$file <- as.numeric(gsub(".*_(\\d+).csv", "\\1", samples_output[i]))
    # create a regex to get the number 1818 from "dev/mem_profiling/samples_output/sample_mem_create_tseries_1818.csv" 
    return(samples[[i]])
})

# merge all samples
samples <- do.call(rbind, samples)

# reomve all special characters from samples$stack.2
samples$stack.2 <- gsub("[^[:alnum:]]", "", samples$stack.2)

# plot stack.2 by vsize.large max
head(samples)

samples_grouped <- 
    samples %>% 
    group_by(stack.2, file) %>% 
    dplyr::summarize(vsize.large = max(vsize.large))

View(samples_grouped)
# convert vsize from bytes to MB
samples_grouped$vsize.large <- samples_grouped$vsize.large / 1024 / 1024

head(samples_grouped)

ggplot(samples_grouped, aes(x = stack.2, y = vsize.large, color = file)) +
    geom_point() +
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red") +
    ylab("vsize.large (MB)") +
    theme(axis.text.x = element_text(angle = 90))


subset(samples_grouped, stack.2 == "pgx.createPGX")