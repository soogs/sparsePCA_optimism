# sparse eigenvectors vs correlated components #

library(sparks)
library(mvtnorm)
library(MASS)

# 2-by-2 conditions #
  # a. low dimension, high sample size (ldhss) // sparse eigenvectors
  # b. hdlss // sparse eigenvectors
  # c. ldhss // correlated components
  # d. hdlss // correlated components
  
# ldhss setup #
  # p = 10
  # n = 300
  # 2 eigenvectors
    
  v1 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0.9, 0.9)
  v2 <- c(0, 0, 0, 0, 1, 1, 1, 1, -0.3, 0.3)
  
  vstar <- cbind(v1, v2)
  vrest <- matrix(runif(n = 80, min = -1, max = 1), nrow = 10, ncol = 8)
  vstar <- cbind(vstar, vrest)
  
  ld_eigenortho <- qr.Q(qr(vstar))
  
  ld_eigenvalues <- c(200, 100, 50, 50, 6, 5, 4, 3, 2, 1)
  # this is the eigenvalues specified in the paper
  
  # ld_eigenvalues <- c(200, 100, rep(0.000001,8))

# hdlss setup
  # p = 500
  # n = 50
  # 2 eigenvectors
  # 1st to 10th elements of the first eigenvector is 1
  # 11th to 20th element of the second eigenvector is 1
  # the rest is 0
  
  hd_eigen <- matrix(0, nrow = 500, ncol =500)
  
  hd_eigen[1:10,1] <- 1
  hd_eigen[11:20,2] <- 1
  
  hd_eigenrest <- matrix(runif(n = (500*(500-2)), min = -1, max = 1), nrow = 500, ncol = 498)
  
  hd_eigenfull <- cbind(hd_eigen, hd_eigenrest)
  
  hd_eigenortho <- qr.Q(qr(hd_eigenfull))
  
  hd_eigenvalues <- c(400, 300, rep(1, 498))
  
  # hdshen_cov <- hd_eigenortho %*% diag(hd_eigenvalues) %*% t(hd_eigenortho)

# covariance matrices for the correlated components #
  ldcov <- diag(1,10)
  ldcov[1,2] <- 0.5
  ldcov[2,1] <- 0.5
  ldcov[ldcov == 0] <- 0.3
  
  hdcov <- diag(1,500)
  hdcov[1,2] <- 0.5
  hdcov[2,1] <- 0.5
  hdcov[hdcov == 0] <- 0.3
  
  
  
  
# vectors to save the results of simulation
  ldshen_tuckers <- c()
  Tldshen_tuckers <- c() # check the tucker congruence of the T matrix
  ldshen_nonhits <- c()
  ldshen_vars <- c()
  ldshen_loss <- c()
  
  hdshen_tuckers <- c()
  Thdshen_tuckers <- c()
  hdshen_nonhits <- c()
  hdshen_vars <- c()
  hdshen_loss <- c()
  
  ldcor_tuckers <- c()
  Tldcor_tuckers <- c()
  ldcor_nonhits <- c()
  ldcor_vars <- c()
  ldcor_loss <- c()
  
  hdcor_tuckers <- c()
  Thdcor_tuckers <- c()
  hdcor_nonhits <- c()
  hdcor_vars <- c()
  hdcor_loss <- c()
  
# results for the multistart #
  mult_ldcor_tuckers <- c()
  mult_Tldcor_tuckers <- c()
  mult_ldcor_nonhits <- c()
  mult_ldcor_vars <- c()
  mult_ldcor_loss <- c()
  
  mult_hdcor_tuckers <- c()
  mult_Thdcor_tuckers <- c()
  mult_hdcor_nonhits <- c()
  mult_hdcor_vars <- c()
  mult_hdcor_loss <- c()



# setting the seed
set.seed(11)
shenrandom <- sample(100000, 1000)

for(i in 1:100){
  
  # ldhss // sparse eigenvectors #
  set.seed(shenrandom[i])
  ld_z <- mvrnorm(n = 100, mu = rep(0,10), Sigma = diag(10), empirical = TRUE)
  
  ld_z_ortho <- qr.Q(qr(ld_z))
  
  ld_z <- ld_z_ortho*sqrt(100)
  
  ldshen_dat <- t(ld_eigenortho %*% diag(sqrt(ld_eigenvalues)) %*% t(ld_z))
  
  ldshen_dat <- VAFcontrol(data = ldshen_dat, seed = shenrandom[i], VAFx = 0.7)
  
  # STANDARDIZE OR NOT #
  # ldshen_dat <- scaleData(ldshen_dat)
  
  ldshen_lasso <- findLasso(dat = ldshen_dat, zeros = c(4,4), R = 2, whichfunction = "rsvd_spca", init = 100, maxiterOut = 1000, maxiterIn = 1000)
  
  ldshen_result <- rsvd_spca(dat = ldshen_dat, R = 2, lambda = ldshen_lasso$lasso, penalty = "soft", maxiter = 1000, inits = "SVD")
  
  ldshen_tucker <- RegularizedSCA::TuckerCoef(ld_eigenortho[,1:2], ldshen_result$V)
  
  ldshen_tuckers[i] <- as.numeric(ldshen_tucker$tucker_value)
  
  Tldshen_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(ld_z[,1:2], ldshen_result$U)$tucker_value)
  
  ldshen_nonhits[i] <- sum((abs(ldshen_result$V[,ldshen_tucker$perm]) > 1e-7) + (abs(ld_eigenortho[,1:2]) > 1e-7) == 2) / (12)
  
  ldshen_vars[i] <- sum((ldshen_dat - ldshen_result$U %*% t(ldshen_result$V))^2) / sum(ldshen_dat^2)
  
  ldshen_loss[i] <- ldshen_result$loss
  
  # hdlss // sparse eigenvectors #
  set.seed(shenrandom[i])
  hd_z <- rmvnorm(n = 500, mean = rep(0,500), sigma = diag(500))
  
  hd_z_ortho <- qr.Q(qr(hd_z))
  
  hd_z <- hd_z_ortho[1:50,]
  
  hd_z <- normalize(hd_z)
  
  hd_z <- hd_z*sqrt(50)
  
  # i think it may be impossible to column-orthogonalize
  # a matrix with size p > n ?
  
  hdshen_dat <- t(hd_eigenortho %*% diag(sqrt(hd_eigenvalues)) %*% t(hd_z))
  
  hdshen_dat <- VAFcontrol(data = hdshen_dat, seed = shenrandom[i], VAFx = 0.7)
  
  # STANDARDIZE OR NOT #
  # hdshen_dat <- scaleData(hdshen_dat)
  
  hdshen_lasso <- findLasso(dat = hdshen_dat, zeros = c(490, 490), R = 2, whichfunction = "rsvd_spca", init = 100, maxiterOut = 1000, maxiterIn = 1000)
  
  hdshen_result <- rsvd_spca(dat = hdshen_dat, R = 2, lambda = hdshen_lasso$lasso, penalty = "soft", maxiter = 1000, inits = "SVD")
  
  hdshen_tucker <- RegularizedSCA::TuckerCoef(hd_eigenortho[,1:2], hdshen_result$V)
  
  hdshen_tuckers[i] <- as.numeric(hdshen_tucker$tucker_value)
  
  Thdshen_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(hd_z[,1:2], hdshen_result$U)$tucker_value)
  
  hdshen_nonhits[i] <- sum((abs(hdshen_result$V[,hdshen_tucker$perm]) > 1e-7) + (abs(hd_eigenortho[,1:2]) > 1e-7) == 2) / (20)
  
  hdshen_vars[i] <- sum((hdshen_dat - hdshen_result$U %*% t(hdshen_result$V))^2) / sum(hdshen_dat^2)
  
  hdshen_loss[i] <- hdshen_result$loss
  
  # ldhss // correlated components #
  set.seed(shenrandom[i])
  ld_z_cor <- mvrnorm(n = 100, mu = rep(0,10), Sigma = ldcov, empirical = TRUE)
  
  ld_z_cor <- normalize(ld_z_cor)
  
  ld_z_cor <- ld_z_cor*sqrt(100)
    
  ldcor_dat <- t(ld_eigenortho %*% diag(sqrt(ld_eigenvalues)) %*% t(ld_z_cor))
  
  ldcor_dat <- VAFcontrol(data = ldcor_dat, seed = shenrandom[i], VAFx = 0.7)
  
  # STANDARDIZE OR NOT #
  # ldshen_dat <- scaleData(ldshen_dat)
  
  ldcor_lasso <- findLasso(dat = ldcor_dat, zeros = c(4,4), R = 2, whichfunction = "rsvd_pca", init = 100, maxiterOut = 1000, maxiterIn = 1000)
  
  ldcor_result <- rsvd_spca(dat = ldcor_dat, R = 2, lambda = ldcor_lasso$lasso, penalty = "soft", inits = "SVD", maxiter = 1000)
  
  ldcor_tucker <- RegularizedSCA::TuckerCoef(ld_eigenortho[,1:2], ldcor_result$V)
  
  ldcor_tuckers[i] <- as.numeric(ldcor_tucker$tucker_value)
  
  Tldcor_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(ld_z_cor[,1:2], ldcor_result$U)$tucker_value)
  
  ldcor_nonhits[i] <- sum((abs(ldcor_result$V[,ldcor_tucker$perm]) > 1e-7) + (abs(ld_eigenortho[,1:2]) > 1e-7) == 2) / (12)
  
  ldcor_vars[i] <- sum((ldcor_dat - ldcor_result$U %*% t(ldcor_result$V))^2) / sum(ldcor_dat^2)
  
  ldcor_loss[i] <- ldcor_result$loss
  
  # hdlss // correlated components #
  set.seed(shenrandom[i])
  hd_z_cor <- rmvnorm(n = 50, mean = rep(0, 500), sigma = hdcov)
  
  hd_z_cor <- normalize(hd_z_cor)
  
  hd_z_cor <- hd_z_cor*sqrt(50)
  
  hdcor_dat <- t(hd_eigenortho %*% diag(sqrt(hd_eigenvalues)) %*% t(hd_z_cor))
  
  hdcor_dat <- VAFcontrol(data = hdcor_dat, seed = shenrandom[i], VAFx = 0.7)
  
  # STANDARDIZE OR NOT #
  # hdshen_dat <- scaleData(hdshen_dat)
  
  hdcor_lasso <- findLasso(dat = hdcor_dat, zeros = c(490, 490), R = 2, whichfunction = "rsvd_spca", init = 100, maxiterOut = 1000, maxiterIn = 1000)
  
  hdcor_result <- rsvd_spca(dat = hdcor_dat, R = 2, lambda = hdcor_lasso$lasso, penalty = "soft", maxiter = 1000, inits = "SVD")
  
  hdcor_tucker <- RegularizedSCA::TuckerCoef(hd_eigenortho[,1:2], hdcor_result$V)
  
  hdcor_tuckers[i] <- as.numeric(hdcor_tucker$tucker_value)
  
  Thdcor_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(hd_z_cor[,1:2], hdcor_result$U)$tucker_value)
  
  hdcor_nonhits[i] <- sum((abs(hdcor_result$V[,hdcor_tucker$perm]) > 1e-7) + (abs(hd_eigenortho[,1:2]) > 1e-7) == 2) / (20)
  
  hdcor_vars[i] <- sum((hdcor_dat - hdcor_result$U %*% t(hdcor_result$V))^2) / sum(hdcor_dat^2)
  
  hdcor_loss[i] <- hdcor_result$loss
  
  
  # MULTISTART // ldhss // correlated components #
  mult_ldcor_result_all  <- rsvd_spca(dat = ldcor_dat, R = 2, lambda = ldcor_lasso$lasso, penalty = "soft", inits = "multistart", maxiter = 1000, nrstart = 500)
  
  mult_ldcor_result <- mult_ldcor_result_all$bunch[[which.min(mult_ldcor_result_all$LOSS)]]
  
  mult_ldcor_tucker <- RegularizedSCA::TuckerCoef(ld_eigenortho[,1:2], mult_ldcor_result$V)
  
  mult_ldcor_tuckers[i] <- as.numeric(mult_ldcor_tucker$tucker_value)
  
  mult_Tldcor_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(ld_z_cor[,1:2], mult_ldcor_result$U)$tucker_value)
  
  mult_ldcor_nonhits[i] <- sum((abs(mult_ldcor_result$V[,mult_ldcor_tucker$perm]) > 1e-7) + (abs(ld_eigenortho[,1:2]) > 1e-7) == 2) / (12)
  
  mult_ldcor_vars[i] <- sum((ldcor_dat - mult_ldcor_result$U %*% t(mult_ldcor_result$V))^2) / sum(ldcor_dat^2)
  
  mult_ldcor_loss[i] <- mult_ldcor_result$loss
  
  
  # MULTISTART // hdlss // correlated components #
  mult_hdcor_result_all <- rsvd_spca(dat = hdcor_dat, R = 2, lambda = hdcor_lasso$lasso, penalty = "soft", maxiter = 1000, inits = "multistart", nrstart = 500)
  
  mult_hdcor_result <- mult_hdcor_result_all$bunch[[which.min(mult_hdcor_result_all$LOSS)]]
  
  mult_hdcor_tucker <- RegularizedSCA::TuckerCoef(hd_eigenortho[,1:2], mult_hdcor_result$V)
  
  mult_hdcor_tuckers[i] <- as.numeric(mult_hdcor_tucker$tucker_value)
  
  mult_Thdcor_tuckers[i] <- as.numeric(RegularizedSCA::TuckerCoef(hd_z_cor[,1:2], mult_hdcor_result$U)$tucker_value)
  
  mult_hdcor_nonhits[i] <- sum((abs(mult_hdcor_result$V[,mult_hdcor_tucker$perm]) > 1e-7) + (abs(hd_eigenortho[,1:2]) > 1e-7) == 2) / (20)
  
  mult_hdcor_vars[i] <- sum((hdcor_dat - mult_hdcor_result$U %*% t(mult_hdcor_result$V))^2) / sum(hdcor_dat^2)
  
  mult_hdcor_loss[i] <- mult_hdcor_result$loss
  
}
