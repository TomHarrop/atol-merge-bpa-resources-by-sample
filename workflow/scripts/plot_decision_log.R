#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

decision_log_file <- "results/decision_log.csv"

decision_log <- fread(decision_log_file, na.strings = "None")
kept_pd <- melt(
  decision_log[dataset_accepted == TRUE],
  id.vars = "id",
  measure.vars = c("bioplatforms_project_id", "data_context", "platform"),
  variable.name = "atol_key"
)

decision_log[, length(unique(id))]
kept_pd[, length(unique(id))]

gp <- ggplot(kept_pd, aes(x = value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(NULL) +
  facet_wrap(vars(atol_key), scales = "free") +
  geom_histogram(stat = "count")

ggsave("test/kept_datasets.pdf", width = 16, height=9, units = "in" )
