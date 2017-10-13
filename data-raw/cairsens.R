cairsens46 <- data.table::fread('data-raw/cairsens46-22-09-2017.csv')
cairsens47 <- data.table::fread('data-raw/cairsens47-22-09-2017.csv')
devtools::use_data(cairsens46, overwrite = TRUE)
devtools::use_data(cairsens47, overwrite = TRUE)
