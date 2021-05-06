## code to prepare `StandardData` dataset goes here

StandardData <- c(rep(3,30), 2,2,1.5,1.5, 1.2,1.2,rep(0.8,4),rep(0.5,10),rep(0.1,50))

usethis::use_data(StandardData, overwrite = TRUE)

