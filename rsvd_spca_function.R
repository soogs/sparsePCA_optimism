# rSVD-sPCA implementation #

# the soft-thresholding or scad-threshodling can be done independently
# because this is a deflation method
# so this is univariate 

normalize <- function(MATRIX){
  norm_inluding_zeros <- function(x){
    if(sum(x==0) == length(x)){
      result <- 1
    } else {
      result <- norm(as.matrix(x),"F")
    }
    return(result)
  }
  apply(as.matrix(MATRIX),2,function(x){x/norm_inluding_zeros(x)})
}

softVec <- function(alpha, lambda){
  big <- alpha > lambda
  small <- alpha < (-1*lambda)
  between <- (-1*lambda) <= alpha & alpha <= lambda
  
  result <- alpha
  
  result[big] <- result[big] - lambda
  result[small] <- result[small] + lambda
  result[between] <- 0
  
  return (result)
}

scadVec <- function(alpha, lambda, gamma){
  
  # if abs(alpha) <= 2*lambda
  # softbig <- alpha > 2*lambda
  # softsmall <- alpha < (-1*2*lambda)
  softbetween <- (-2*lambda) <= alpha & alpha <= 2*lambda
  
  # if 2lambda < abs(alpha) <= gamma * lambda 
  scadpos <- (2*lambda < alpha & alpha  <= gamma * lambda)
  scadneg <- ((-2*lambda) > alpha & alpha  >= gamma * lambda)
  
  # if abs(alpha) > gamma * lambda
  scadbig <- abs(alpha) > gamma*lambda
  
  result <- alpha
  
  result[softbetween] <- 0
  result[scadpos] <- ((gamma-1)*result[scadpos] - (gamma*lambda)) / (gamma - 2)
  result[scadneg] <- ((gamma-1)*result[scadneg] + (gamma*lambda)) / (gamma - 2)
  result[scadbig] <- result[scadbig]
  
  return (result)
}

rsvd_spca <- function(dat, R, lambda, penalty = c("soft", "scad"), 
                      ridge = 1e-6, maxiter, scadgamma,
                      inits = c("SVD", "oracle", "multistart"),
                      nrstart = NULL, oracle = NULL){
    
  if((inits == "SVD" || inits == "oracle") && !is.null(nrstart)){
    stop("Cannot define nrstart while inits = SVD or oracle")
  }
  
  if(inits == "oracle" && is.null(oracle)){
    stop("Please specify the oracle initial values for the loadings")
  }
  
  # svd on dat
  svdobj<-svd(dat)
    
  if (inits == "SVD"){
    v <- svdobj$v[,1:R]
    nrstart <- 1
  } else if (inits == "oracle"){
    v <- oracle
    nrstart <- 1
  } 
  
  resultbunch <- list()
  LOSS <- c()
  
  Vresults <- matrix(0, nrow = ncol(dat), ncol = R)
  Uresults <- matrix(0, nrow = nrow(dat), ncol = R)
  
  for (nr in 1:nrstart){
    
    if (inits == "multistart"){
      p <- ncol(dat)
      v <- matrix(stats::runif(n = p*R, min = -5, max = 5), nrow = p, ncol = R)
    }
  
    alliter <- c()
  
    for (r in 1:R){
      
      if (r == 1){
        dat_iter <- dat
      }
      
      uold <- normalize(dat_iter %*% v)[,r]
      iter <- 0
      
      converge <- FALSE
      lossold <- sum(dat^2)
    
      while (!converge && iter < maxiter){
        iter <- iter + 1 
        
        vnew <- softVec(alpha = t(dat_iter) %*% uold, lambda = lambda[r])
        
        Xv <- dat_iter %*% vnew
        
        unew <- normalize(Xv)
        
        lossnew <- sum((dat_iter - uold %*% t(vnew))^2) + sum(abs(vnew))*lambda[r] + sum(vnew^2)*ridge
        
        lossdiff <- abs(lossold - lossnew)
        
        if (lossdiff > 0.00001){
          uold <- unew
          lossold <- lossnew
        } else {
          converge <- TRUE
        }
      }
      
      residual <- dat_iter - (unew %*% t(vnew))
      
      dat_iter <- residual
      
      Vresults[,r] <- vnew
      Uresults[,r] <- unew
      alliter[r] <- iter
    }
    
    loss <- sum((dat - Uresults %*% t(Vresults))^2) + (ridge * sum(Vresults^2)) + (sum(abs(Vresults) %*% diag(lambda)))
    
    result <- list(V = Vresults, U = Uresults, iter = alliter, loss = loss)
    # class(obj) <- "spca"
    
    resultbunch[[nr]] <- result
    LOSS[nr] <- loss
  }
  
  if(nr == 1){
    allresult <- resultbunch[[1]]
  } else {
    allresult <- list(bunch = resultbunch, LOSS = LOSS)
  }
  
  return (allresult)
  
}







# testing ####
P <- matrix(0,30,3)
P[1:10,1] <- 1
P[11:20,2] <- 1
P[11:30,3] <- 1

t(P) %*% P

Pblock <- P[11:20,2:3]
Pblock2 <- qr.Q(qr(Pblock))

P[11:20,2:3] <- Pblock2

t(P) %*% P



P <- normalize(P)

t(P) %*% P

# specifying the T matrix (norm 1: U matrix)
Tmat <- MASS::mvrnorm(n = 100, mu = rep(0,3), Sigma = diag(3), empirical=TRUE)

# orthogonalizing the Tmat
Tmat <- qr.Q(qr(Tmat))

t(Tmat) %*% Tmat

# X = UDV'
# X = TP'
# T = UD 
# D = diagonal matrix of PC vairance

D <- diag(c(50, 30, 20))

X<- Tmat %*% D %*% t(P)

n <- 100; p <- 30; propNoise <- 0.2
E <- matrix( rnorm(n*p,0,1), n, p )
g = sqrt((var(as.vector(X))*propNoise) / (var(as.vector(E)) * (1 - propNoise)))
X2 <- X + E*g

SStrue <- var(as.vector(X))
SSX <- var(as.vector(X2))

1 - (SStrue/SSX)

svd(X2)$v[,1:3]

spca1 <- elasticnet::spca(x = X2, K = 3, para = c(10, 10, 20), type = "predictor", sparse = "varnum")

uinit <- svd(X2)$u[,1:3]
vinit <- svd(X2)$v[,1:3]

shen1 <- rsvd_spca(dat = X2, R = 3, lambda = c(3, 3, 1), penalty = "soft", maxiter = 10000, inits = "SVD")
shen2 <- rsvd_spca(dat = X2, R = 3, lambda = c(3, 3, 1), penalty = "soft", maxiter = 10000, inits = "multistart", nrstart = 500)

shen1$loss
shen2$bunch[[which.min(shen2$LOSS)]]


# scad penalty for later ####
# shenhuang_scad <- function(dat, uinit, vinit, lambda, maxiter, R){
#   
#   V <- vinit
#   U <- uinit
#   alliter <- c()
#   
#   for (r in 1:R){
#     uold <- uinit[,r]
#     iter <- 0
#     
#     converge <- FALSE
#     lossold <- sum(dat^2)
#     
#     while (!converge && iter < maxiter){
#       iter <- iter + 1 
#       
#       vnew <- scadVec(alpha = t(dat) %*% uold, lambda = lambda[r], gamma = 3.7)
#       
#       unew <- (dat %*% vnew) / sqrt(sum((dat %*% vnew)^2))
#       
#       lossnew <- sum((dat - uold %*% t(vnew))^2)
#       
#       lossdiff <- lossold - lossnew
#       
#       if (lossdiff > 0.00001){
#         uold <- unew
#         lossold <- lossnew
#       } else {
#         converge <- TRUE
#       }
#     }
#     
#     residual <- dat - (unew %*% t(vnew))
#     
#     dat <- residual
#     
#     V[,r] <- vnew
#     U[,r] <- unew
#     alliter[r] <- iter
#   }
#   
#   result <- list(V = V, U = U, iter = alliter)
#   return(result)
# }
# 
