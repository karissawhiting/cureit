
# adding a competing risks outcome to the gtsummary::trial dataset
set.seed(1123)
trial <-
  gtsummary::trial %>%
  mutate(death = as.numeric(death))

usethis::use_data(trial, overwrite = TRUE)