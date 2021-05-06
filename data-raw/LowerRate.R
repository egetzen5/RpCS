## code to prepare `LowerRate` dataset goes here

LowerRate <- c(rep(2.5,30), 2,2,1.5,1.5, 1.2,1.2,rep(0.8,4),rep(0.5,10),rep(0.1,50))

usethis::use_data(LowerRate, overwrite = TRUE)
