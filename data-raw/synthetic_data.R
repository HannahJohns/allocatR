rm(list=ls())

sample_size <- 50

set.seed(361458)

synthetic_covariate_data <- data.frame(
  cov1 = runif(sample_size),
  cov2 = rexp(sample_size),
  cov3 = sample(c("A", "B", "C"), size = sample_size, replace = T),
  cov4 = sample(c("Low", "Medium", "High"), size = sample_size, replace = T)
)

synthetic_covariate_data$cov3 <- factor(synthetic_covariate_data$cov3)
synthetic_covariate_data$cov4 <- factor(synthetic_covariate_data$cov4, levels =c("Low","Medium", "High"), ordered = T)

usethis::use_data(synthetic_covariate_data, overwrite = TRUE)
