# findLasso function #

# a function that figures out the lambda value for each component,
# given the number of zeros that you want for each component

# allows sparse PCA functions of:
# spca.defiled
# rsvd_spca

# i want to match the number of zeros for all my components
# i iterate the process, until all of the components have the right number of zeros

source("/home/soogs/Desktop/Rstudio Github/SPCA-W-P/spca-defiled function.R")
source("/home/soogs/Desktop/Rstudio Github/sparsePCA_optimism/rsvd_spca_function.R")

findLasso <- function(dat, zeros, R, whichfunction = c("spca.defiled", "rsvd_spca"), init, ridge = 1e-6, maxiterOut, maxiterIn){
  # zeros should be a vector
  
  estimatedzeros <- rep(0,R)
  lasso <- rep(init,R)
  
  converged <- FALSE
  
  iterOut <- 0
  
  while(abs(zeros - estimatedzeros) > 0 && iterOut <= maxiterOut ){
    iterOut <- iterOut + 1
    
    for (j in 1:R){
      iterIn <- 0
      
      up <- init
      down <- 0
      
      estimatedzero <- 0
      while(abs(zeros[j] - estimatedzero) > 0 && iterIn <= maxiterIn){
        iterIn <- iterIn + 1
        lasso[j] <- (down + up)/2
        
        # lasso[-j] <- 0
        
        if(whichfunction == "spca.defiled"){
          fit <- spca.defiled(x = dat, K = R, para = lasso, type = "predictor", sparse = "penalty", inits = "SVD", lambda = ridge)
          estimatedzero <- sum(abs(fit$Wraw[,j]) < 1e-06)
        } else {
          fit <- rsvd_spca(dat = dat, R = R, lambda = lasso, penalty = "soft", ridge = ridge, maxiter = 1000, inits = "SVD")
          estimatedzero <- sum(abs(fit$V[,j]) < 1e-06)
        }
        
        
        if(zeros[j] > estimatedzero){
          down  <- lasso[j]
          # if the estimated zeros are not enough,
          # pull up the 'down' 
        } else if (zeros[j] < estimatedzero){
          up  <- lasso[j]
          # if the estimated zeros are more than enough,
          # pull down the 'up'
        } 
        # else (don't do anything)
        
        print(round(lasso,10))
      }
      
      if(j == R){
        if(whichfunction == "spca.defiled"){
          estimatedzeros <- apply((abs(fit$Wraw) < 1e-06),2,sum)
        } else {
          estimatedzeros <- apply((abs(fit$V) < 1e-06),2,sum)
        }
      }
    }
    
  }
  
  if( iterOut < maxiterOut && iterIn < maxiterIn ){
    converged <- TRUE
  } 
  return(list(lasso = lasso, converged = converged))
}