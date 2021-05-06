## code to prepare `EarlyIntervention` dataset goes here

EarlyIntervention <- c(rep(3,25), 2,2,1.5,1.5, 1.2,1.2,rep(0.8,4),rep(0.5,10),rep(0.1,55))

usethis::use_data(EarlyIntervention, overwrite = TRUE)
