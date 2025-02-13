#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

decision_log_file <- "results/decision_log.csv"

decision_log <- fread(decision_log_file)
decision_log[, kept := all(keep_project, keep_context, keep_platform), by = id]
kept_pd <- melt(
  decision_log[keep_dataset == TRUE],
  id.vars = "id",
  measure.vars = c("bioplatforms_project", "data_context", "platform"),
  variable.name = "atol_key"
)


ggplot(kept_pd, aes(x = value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~atol_key, scales = "free") +
  geom_histogram(stat = "count")
