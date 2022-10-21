library(dplyr)
library(forcats)
library(ISwR)
library(labelled)

melanoma <- ISwR::melanom %>%
  transmute(
    id = no,
    status = case_when(
      status %in% c(1) ~ 1,
      status %in% c(2, 3) ~ 0
    ),
    days,
    years = days/365,
    tumor_thickness =
      case_when(
        thick <= 129 ~ "≤ 129",
        thick <= 322 ~ "> 129 & ≤ 322",
        thick > 322 ~ "> 322"
      ),
    tumor_thickness = fct_relevel(tumor_thickness, "≤ 129", "> 129 & ≤ 322"),
    sex = case_when(
      sex == 2 ~ "Male",
      sex == 1 ~ "Female"
    ),
    ulceration = case_when(
      ulc == 1 ~ "Yes",
      ulc == 2 ~ "No"
    )
  )

melanoma <- melanoma %>%
  set_variable_labels(.labels = snakecase::to_title_case(names(.)) )


usethis::use_data(melanoma, overwrite = TRUE)
