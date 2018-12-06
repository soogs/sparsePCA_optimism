# sparse eigenvectors vs sparse weights generation #

library(sparks)
library(mvtnorm)
library(MASS)

  # a. low dimension, high sample size (ldhss) // sparse eigenvectors
  # b. hdlss // sparse eigenvectors
  # c. ldhss // sparse weights
  # d. hdlss // sparse weights

# shen sparse eigenvectors low dimensional #
  # p = 10
  # n = 300
  # 2 eigenvectors

  v1 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0.9, 0.9)
  v2 <- c(0, 0, 0, 0, 1, 1, 1, 1, -0.3, 0.3)

  vstar <- cbind(v1, v2)
  vrest <- matrix(runif(n = 80, min = -1, max = 1), nrow = 10, ncol = 8)
  vstar <- cbind(vstar, vrest)

  ldshen_eigenortho <- qr.Q(qr(vstar))
  
  # ldshen_eigenvalues <- c(200, 100, 50, 50, 6, 5, 4, 3, 2, 1)

  ldshen_eigenvalues <- c(200, 100, rep(0.0000001,8))



# shen sparse eigenvectors high dimensional #
  # p = 500
  # n = 50
  # 2 eigenvectors
  # 1st to 10th elements of the first eigenvector is 1
  # 11th to 20th element of the second eigenvector is 1
  # the rest is 0

  hdshen_eigen <- matrix(0, nrow = 500, ncol =500)
  
  hdshen_eigen[1:10,1] <- 1
  hdshen_eigen[11:20,2] <- 1
  
  hdshen_eigenrest <- matrix(runif(n = (500*(500-2)), min = -1, max = 1), nrow = 500, ncol = 498)
  
  hdshen_eigenfull <- cbind(hdshen_eigen, hdshen_eigenrest)
  
  hdshen_eigenortho <- qr.Q(qr(hdshen_eigenfull))
  
  hdshen_eigenvalues <- c(400, 300, rep(0.0000001, 498))
  
  # hdshen_cov <- hdshen_eigenortho %*% diag(hdshen_eigenvalues) %*% t(hdshen_eigenortho)

# vectors to save the results
  ldshen_tuckers <- c()
  Tldshen_tuckers <- c()
  ldshen_nonhits <- c()
  ldshen_vars <- c()
  
  hdshen_tuckers <- c()
  Thdshen_tuckers <- c()
  hdshen_nonhits <- c()
  hdshen_vars <- c()
  
  ldxwp_tuckers <- c()
  Tldxwp_tuckers <- c()
  ldxwp_nonhits <- c()
  ldxwp_vars <- c()
  ldxwp_loss <-c()
  
  hdxwp_tuckers <- c()
  Thdxwp_tuckers <- c()
  hdxwp_nonhits <- c()
  hdxwp_vars <- c()
  hdxwp_loss <- c()
  
  mult_ldxwp_tuckers <- c()
  mult_Tldxwp_tuckers <- c()
  mult_ldxwp_nonhits <- c()
  mult_ldxwp_vars <- c()
  mult_ldxwp_loss <- c()
  
  mult_hdxwp_tuckers <- c()
  mult_Thdxwp_tuckers <- c()
  mult_hdxwp_nonhits <- c()
  mult_hdxwp_vars <- c()
  mult_hdxwp_loss <- c()

# setting the seed
  set.seed(22)
  xwprandom <- sample(100000, 1000)

for(i in 1:100){
  
  # ldhss // sparse eigenvectors #
  set.seed(xwprandom[i])
  ld_z <- mvrnorm(n = 100, mu = rep(0,10), Sigma = diag(10), empirical = TRUE)
  
  ld_z_ortho <- qr.Q(qr(ld_z))
  
  ld_z <- ld_z_ortho*sqrt(100)
  
  ldshen_dat <- t(ldshen_eigenortho %*% diag(sqrt(ldshen_eigenvalues)) %*% t(ld_z))
  
  ldshen_dat <- VAFcontrol(data = ldshen_dat, seed = xwprandom[i], VAFx = 0.7)
  
  # STANDARDIZE OR NOT #
  # ldshen_dat <- scaleData(ldshen_dat)
  
  ldshen_result <- spca_adj(x = ldshen_dat, K = 2, para = c(6,6), type = "predictor", sparse = "varnum", inits = "SVD")
  
  ldshen_tucker <- RegularizedSCA::TuckerCoef(ldshen_eigenortho[,1:2], ldshen_result$loadings)
  
  ldshen_tuckers[i] <- as.numeric(ldshen_tucker$tucker_value)
  
  Tldshen_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(ldshen_dat %*% ldshen_eigenortho[,1:2], ldshen_dat %*% ldshen_result$Wraw)$tucker_value)
  
  ldshen_nonhits[i] <- sum((abs(ldshen_result$loadings[,ldshen_tucker$perm]) > 1e-7) + (abs(ldshen_eigenortho[,1:2]) > 1e-7) == 2) / (12)
  
  ldshen_vars[i] <- sum((ldshen_dat - ldshen_dat %*% ldshen_result$Wraw %*% t(ldshen_result$Pmat))^2) / sum(ldshen_dat^2)
  
  
  # hdlss // sparse eigenvectors #
  set.seed(xwprandom[i])
  hd_z <- rmvnorm(n = 500, mean = rep(0,500), sigma = diag(500))
  
  hd_z_ortho <- qr.Q(qr(hd_z))
  
  hd_z <- hd_z_ortho[1:50,]
  
  hd_z <- normalize(hd_z)
  
  hd_z <- hd_z*sqrt(50)
  
  hdshen_dat <- t(hdshen_eigenortho %*% diag(sqrt(hdshen_eigenvalues)) %*% t(hd_z))
  
  hdshen_dat <- VAFcontrol(data = hdshen_dat, seed = xwprandom[i], VAFx = 0.7)
  
  # STANDARDIZE OR NOT #
  # hdshen_dat <- scaleData(hdshen_dat)
  
  hdshen_result <- spca_adj(x = hdshen_dat, K = 2, para = c(10,10), type = "predictor", sparse = "varnum", inits = "SVD")
  
  hdshen_tucker <- RegularizedSCA::TuckerCoef(hdshen_eigenortho[,1:2], hdshen_result$loadings)
  
  hdshen_tuckers[i] <- as.numeric(hdshen_tucker$tucker_value)
  
  Thdshen_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(hdshen_dat %*% hdshen_eigenortho[,1:2], hdshen_dat %*% hdshen_result$Wraw)$tucker_value)
  
  hdshen_nonhits[i] <- sum((abs(hdshen_result$loadings[,hdshen_tucker$perm]) > 1e-7) + (abs(hdshen_eigenortho[,1:2]) > 1e-7) == 2) / (20)
  
  hdshen_vars[i] <- sum((hdshen_dat - hdshen_dat %*% hdshen_result$Wraw %*% t(hdshen_result$Pmat))^2) / sum(hdshen_dat^2)
  
  
  
  # ldhss // sparse weights #
  ldxwp <- Wsparsedata(n = 300, p = 10, k = 2, n_zeros = 4, Winit = NULL, empirical = FALSE, seed = xwprandom[i], VAFx = 0.7)
  
  ldxwp_dat <- ldxwp$X
  # STANDARDIZE OR NOT #
  # ldxwp_dat <- scaleData(ldxwp$X)
  
  ldxwp_lasso <- findLasso(dat = ldxwp_dat, zeros = c(6,6), R = 2, whichfunction = "spca_adj", init = 1000, maxiterOut = 5, maxiterIn = 1000)
  
  ldxwp_result <- spca_adj(x = ldxwp_dat, K = 2, para = ldxwp_lasso$lasso, type = "predictor", sparse = "penalty", inits = "SVD")
  
  ldxwp_tucker <- RegularizedSCA::TuckerCoef(ldxwp$W, ldxwp_result$loadings)
  
  ldxwp_tuckers[i] <- as.numeric(ldxwp_tucker$tucker_value)
  
  Tldxwp_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(ldxwp_dat %*% ldxwp$W, ldxwp_dat %*% ldxwp_result$Wraw)$tucker_value)
  
  ldxwp_nonhits[i] <- sum((abs(ldxwp_result$loadings[,ldxwp_tucker$perm]) > 1e-7) + (abs(ldxwp$W) > 1e-7) == 2) / (12)
  
  ldxwp_vars[i] <- sum((ldxwp_dat - ldxwp_dat %*% ldxwp_result$Wraw %*% t(ldxwp_result$Pmat))^2) / sum(ldxwp_dat^2)
  
  ldxwp_loss[i] <- ldxwp_result$loss
  
  # hdlss // sparse weights #
  hdxwp <- Wsparsedata(n = 10, p = 500, k = 2, n_zeros = 490, Winit = NULL, empirical = FALSE, seed = xwprandom[i], VAFx = 0.7)
  
  hdxwp_dat <- hdxwp$X
  # STANDARDIZE OR NOT #
  # hdxwp_dat <- scaleData(hdxwp$X)
  
  hdxwp_lasso <- findLasso(dat = hdxwp_dat, zeros = c(490, 490), R = 2, whichfunction = "spca_adj", init = 1000, maxiterOut = 5, maxiterIn = 1000)
  
  hdxwp_result <- spca_adj(x = hdxwp_dat, K = 2, para = hdxwp_lasso$lasso, type = "predictor", sparse = "penalty", inits = "SVD")
  
  hdxwp_tucker <- RegularizedSCA::TuckerCoef(hdxwp$W, hdxwp_result$loadings)
  
  hdxwp_tuckers[i] <- as.numeric(hdxwp_tucker$tucker_value)
  
  Thdxwp_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(hdxwp_dat %*% hdxwp$W, hdxwp_dat %*% hdxwp_result$Wraw)$tucker_value)
  
  hdxwp_nonhits[i] <- sum((abs(hdxwp_result$loadings[,hdxwp_tucker$perm]) > 1e-7) + (abs(hdxwp$W) > 1e-7) == 2) / (20)
  
  hdxwp_vars[i] <- sum((hdxwp_dat - hdxwp_dat %*% hdxwp_result$Wraw %*% t(hdxwp_result$Pmat))^2) / sum(hdxwp_dat^2)
  
  hdxwp_loss[i] <- hdxwp_result$loss
  
  # MULTISTART // ldhss // sparse weights #
  mult_ldxwp_result_all <- spca_adj(x = ldxwp_dat, K = 2, para = ldxwp_lasso$lasso, type = "predictor", sparse = "penalty", inits = "multistart", nrstart = 50)
  
  mult_ldxwp_result <- mult_ldxwp_result_all$bunch[[which.min(mult_ldxwp_result_all$LOSS)]]
  
  mult_ldxwp_tucker <- RegularizedSCA::TuckerCoef(ldxwp$W, mult_ldxwp_result$loadings)
  
  mult_ldxwp_tuckers[i] <- as.numeric(mult_ldxwp_tucker$tucker_value)
  
  mult_Tldxwp_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(ldxwp_dat %*% ldxwp$W, ldxwp_dat %*% mult_ldxwp_result$Wraw)$tucker_value)
  
  mult_ldxwp_nonhits[i] <- sum((abs(mult_ldxwp_result$loadings[,mult_ldxwp_tucker$perm]) > 1e-7) + (abs(ldxwp$W) > 1e-7) == 2) / (12)
  
  mult_ldxwp_vars[i] <- sum((ldxwp_dat - ldxwp_dat %*% mult_ldxwp_result$Wraw %*% t(mult_ldxwp_result$Pmat))^2) / sum(ldxwp_dat^2)
  
  mult_ldxwp_loss[i] <- mult_ldxwp_result$loss
  
  
  # MULTISTART // hdlss // sparse weights #
  mult_hdxwp_result_all <- spca_adj(x = hdxwp_dat, K = 2, para = hdxwp_lasso$lasso, type = "predictor", sparse = "penalty", init = "multistart", nrstart = 50)
  
  mult_hdxwp_result <- mult_hdxwp_result_all$bunch[[which.min(mult_hdxwp_result_all$LOSS)]]
  
  mult_hdxwp_tucker <- RegularizedSCA::TuckerCoef(hdxwp$W, mult_hdxwp_result$loadings)
  
  mult_hdxwp_tuckers[i] <- as.numeric(mult_hdxwp_tucker$tucker_value)
  
  mult_Thdxwp_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(hdxwp_dat %*% hdxwp$W, hdxwp_dat %*% mult_hdxwp_result$Wraw)$tucker_value)
  
  mult_hdxwp_nonhits[i] <- sum((abs(mult_hdxwp_result$loadings[,mult_hdxwp_tucker$perm]) > 1e-7) + (abs(hdxwp$W) > 1e-7) == 2) / (20)
  
  mult_hdxwp_vars[i] <- sum((hdxwp_dat - hdxwp_dat %*% mult_hdxwp_result$Wraw %*% t(mult_hdxwp_result$Pmat))^2) / sum(hdxwp_dat^2)
  
  mult_hdxwp_loss[i] <- mult_hdxwp_result$loss
  
}
