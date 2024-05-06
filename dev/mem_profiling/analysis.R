library(playbase)
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
# real all samples files under samples_output
samples_output <- list.files("dev/mem_profiling/samples_output_libx", full.names = TRUE)

# read csv
samples <- lapply(samples_output, read.csv)

# add the name of the file (containing number of samples) as a column
samples <- lapply(seq_along(samples), function(i) {
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
    summarise(vsize.large = max(vsize.large))


# convert vsize from bytes to MB
samples_grouped$vsize.large <- samples_grouped$vsize.large / 1024 / 1024

head(samples_grouped)

ggplot(samples_grouped, aes(x = stack.2, y = vsize.large, color = file)) +
    geom_point() +
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red") +
    ylab("vsize.large (MB)") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))  # Rotate x-axis labels for better visibility


subset(samples_grouped, stack.2 == "pgx.createPGX")