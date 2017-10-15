BC <- data.table::fread('data-raw/BC.csv')
devtools::use_data(BC, overwrite = TRUE)
