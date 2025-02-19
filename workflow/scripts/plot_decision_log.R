#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

decision_log_file <- "results/decision_log.csv"

decision_log <- fread(decision_log_file, na.strings = c("None", "", "NA"))

all_pd <- melt(
  decision_log,
  id.vars = "id",
  measure.vars = c("bioplatforms_project_id", "data_context"),
  variable.name = "atol_key"
)

wholeorder <- all_pd[, .N, by = .(value)][order(N, decreasing = TRUE), unique(value)]

all_pd[
  ,
  ordered_value := factor(value, levels = wholeorder)
]

gp1 <- ggplot(all_pd, aes(x = ordered_value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_log10() +
  xlab(NULL) +
  facet_wrap(vars(atol_key), scales = "free") +
  geom_histogram(stat = "count")

ggsave(
  "test/all_datasets.pdf",
  gp1,
  width = 16,
  height = 9,
  units = "in"
)


kept_pd <- melt(
  decision_log[package_accepted == TRUE],
  id.vars = "id",
  measure.vars = c("bioplatforms_project_id", "data_context", "platform"),
  variable.name = "atol_key"
)

keptorder <- kept_pd[, .N, by = .(value)][order(N, decreasing = TRUE), unique(value)]


kept_pd[
  ,
  ordered_value := factor(value, levels = keptorder)
]


decision_log[
  data_context_accepted == TRUE & package_accepted == FALSE & bioplatforms_project_id == "bpa-omg"
]

decision_log[, length(unique(id))]
kept_pd[, length(unique(id))]

gp <- ggplot(kept_pd, aes(x = ordered_value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(NULL) +
  facet_wrap(vars(atol_key), scales = "free") +
  geom_histogram(stat = "count")

ggsave(
  "test/kept_datasets.pdf",
  gp,
  width = 16,
  height = 9,
  units = "in"
)
