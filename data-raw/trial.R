library(dplyr)
library(gtsummary)

set.seed(1123)
trial <-
  gtsummary::trial %>%
  mutate(death = as.numeric(death))

usethis::use_data(trial, overwrite = TRUE)