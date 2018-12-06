# zou strategy vs. sparse weights strategy #

library(sparks)

# zou generation setup #
v1v1 <- matrix(data = 290^2, nrow = 4, ncol = 4)
v2v2 <- matrix(data = 300^2, nrow = 4, ncol = 4)
v1v2 <- matrix(data = 0, nrow = 4, ncol = 4)
v1v3 <- matrix(data = -0.3*290^2, nrow = 4, ncol = 2)
v2v3 <- matrix(data = 0.925*300^2, nrow = 4, ncol = 2)
v3v3 <- matrix(data = 0.3^2*290^2 + 0.925^2*300^2, nrow = 2, ncol = 2)

hi <- cbind(v1v1, v1v2)
hi <- cbind(hi, v1v3)
hi2 <- cbind(v1v2, v2v2, v2v3)
hi3 <- cbind(t(v1v3), t(v2v3), v3v3)

covmat <- rbind(hi, hi2, hi3)

# covmat[,9:10] <- 0
# covmat[9:10,] <- 0

# probable eigenvectors of the covariance matrix #
v1 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
v2 <- c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0)

zou_eigen <- cbind(v1,v2)


# vectors saving information about the wxp data generation
  w_tuckers <- c()
  Tw_tuckers <- c()
  w_hits <- c()
  w_vars <- c()
  w_loss <- c()

# vectors saving information about zou generation
  zou_tuckers <- c()
  Tzou_tuckers <- c()
  zou_hits <- c()
  zou_vars <- c()
  zou_loss <- c()

# vectors to save information about the multistart approach #
  Tmult_tuckers <- c()
  mult_tuckers <- c()
  mult_hits <- c()
  mult_vars <- c()
  mult_loss <- c()

# setting the seed
  set.seed(22)
  xwprandom <- sample(100000, 1000)


for(i in 1:100){
  
  # xwp generation #
  xwp <- Wsparsedata(n = 300, p = 10, k = 2, n_zeros = 6, Winit = NULL, empirical = FALSE, seed = xwprandom[i], VAFx = 0.7)
  
  dat <- xwp$X
  # dat <- scaleData(xwp$X)
  
  # at i = 40, this takes a very long time
  w_lasso <- findLasso(dat = dat, zeros = c(6,6), R = 2, whichfunction = "spca_adj", init = 1000, maxiterOut = 1000, maxiterIn = 1000)$lasso
  
  result <- spca_adj(x = dat, K = 2, para = w_lasso, type = "predictor", sparse = "penalty", inits = "SVD")
  
  wtucker <- RegularizedSCA::TuckerCoef(xwp$W, result$loadings)
  
  w_tuckers[i] <- as.numeric(wtucker$tucker_value)
  
  Tw_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(xwp$Xtrue %*% xwp$W, dat %*% result$Wraw)$tucker_value)
  
  w_hits[i] <- sum((abs(result$loadings[,wtucker$perm]) < 1e-7) + (abs(xwp$W) < 1e-7) == 2) / 12
  
  w_vars[i] <- sum((dat - dat %*% result$Wraw %*% t(result$Pmat))^2) / sum(dat^2)
  
  w_loss[i] <- result$loss
  
  
  
  # zou generation #
  set.seed(xwprandom[i])
  
  zou_dat_TRUE <- matrix(MASS::mvrnorm(n = 300, mu = rep(0,10), Sigma = covmat, empirical = FALSE), nrow = 300, ncol = 10)
  
  zou_dat <- VAFcontrol(data = zou_dat_TRUE, seed = xwprandom[i], VAFx = 0.7)
  
  # zou_dat <- scaleData(zou_dat)
  
  zou_lasso <- findLasso(dat = zou_dat, zeros = c(6,6), R = 2, whichfunction = "spca_adj", init = 1e+8, maxiterOut = 5, maxiterIn = 1000)$lasso
  
  zou_result <- spca_adj(x = zou_dat, K = 2, para = zou_lasso, type = "predictor", sparse = "penalty", inits = "SVD")
  
  zoutucker <- RegularizedSCA::TuckerCoef(zou_eigen, zou_result$loadings)
  
  zou_tuckers[i] <- as.numeric(zoutucker$tucker_value)
  
  Tzou_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(zou_dat %*% zou_eigen, zou_dat %*% zou_result$Wraw)$tucker_value)
  
  zou_hits[i] <- sum((abs(zou_result$loadings[,zoutucker$perm]) < 1e-7) + (abs(zou_eigen) < 1e-7) == 2) / 12
  
  zou_vars[i] <- sum((zou_dat - zou_dat %*% zou_result$Wraw %*% t(zou_result$Pmat))^2) / sum(zou_dat^2)
  
  
  # multistart approach on the sparse weights model #
  
  # mult_result <- spca_adj(x = dat, K = 2, para = rep(4,2), type = "predictor", sparse = "varnum", init = "oracle", oracle = xwp$W + runif(n = 20, min = -0.1, max = 0.1))
  
  mult_result_all <- spca_adj(x = dat, K = 2, para = w_lasso, type = "predictor", sparse = "penalty", init = "multistart", nrstart = 50)
  
  mult_result <- mult_result_all$bunch[[which.min(mult_result_all$LOSS)]]
  
  multtucker <- RegularizedSCA::TuckerCoef(xwp$W, mult_result$loadings)
  
  mult_tuckers[i] <- as.numeric(multtucker$tucker_value)
  
  Tmult_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(dat %*% xwp$W, dat %*% mult_result$Wraw)$tucker_value)
  
  mult_hits[i] <- sum((abs(mult_result$loadings[,multtucker$perm]) < 1e-7) + (abs(xwp$W) < 1e-7) == 2) / 12
  
  mult_vars[i] <- sum((dat - dat %*% mult_result$Wraw %*% t(mult_result$Pmat))^2) / sum(dat^2)
  
  mult_loss[i] <- mult_result$loss
  
}