library(quantreg)

# screening procedure: CQ-SIS
CQSIS <- function(x, y, tau, cn) {
  n <- nrow(x); p <- ncol(x)
  psi <- function(x, tau) (tau - as.numeric(x<0))
  y.tau.hat <- as.numeric(stats::quantile(y, probs=tau))
  
  Margin_Utility <- function(k) {
    xk <- x[,k]
    dhat_func <- function(t) {
      mean(psi(y-y.tau.hat, tau) * as.numeric(xk<t))
    }
    rho.hat.temp <- sapply(xk, dhat_func)
    rho.hat <- mean(rho.hat.temp^2)
    return(rho.hat)
  }
  
  marginal.temp <- sapply(1:p, Margin_Utility)
  subset.ind <- order(marginal.temp, decreasing=TRUE)[1:cn]
  return(subset.ind)
}

# Testing procedure for H0:\beta_{\tau}=0
InferencePLQM <- function(x, z, y, tau, is.split=FALSE, is.screen=FALSE, seed.fix) {
  
  # the sample size and dimension
  n <- length(y)
  d <- ncol(x); q <- ncol(z); p <- d + q
  
  if(is.screen) {
    screening.set <- CQSIS(x=z, y=y, tau=tau, cn=ceiling(n/log(n)))
  }
  if(!is.screen) {
    screening.set <- 1:q
  }
  
  z <- z[,screening.set]
  z.tilde <- cbind(rep(1, n), z)
  
  # p-value without data-splitting
  if(is.split) {
    set.seed(seed.fix)
    
    # data-splitting
    n1 <- floor(n/2)
    ind1 <- sample(1:n, n1, replace=F)
    ind2 <- (1:n)[-ind1]
    
    # construction of standardized test statistics
    Sn_func <- function(ind1, ind2) {
      x1 <- x[ind1,]; x2 <- x[ind2,]
      z1 <- z[ind1,]; z2 <- z[ind2,]
      y1 <- y[ind1]; y2 <- y[ind2]
      z1.tilde <- z.tilde[ind1,]; z2.tilde <- z.tilde[ind2,]
      
      # estimate only gamma
      lasso.model <- rq(y1 ~ z1, tau=tau, method="lasso")
      mz.pred <- z2.tilde%*%coefficients(lasso.model)
      
      # construction of test statistics
      phi_func <- function(u, tau=tau) { tau - ifelse(u<0, 1, 0) }
      phi.vec <- phi_func(as.numeric(y2-mz.pred), tau=tau)
      phi.mat <- outer(phi.vec, phi.vec, "*")
      x.mat <- x2 %*% t(x2)
      Sn <- (sum(phi.mat*x.mat) - sum(diag(phi.mat*x.mat)))/n1
      
      tr.Sigx2.hat <- (sum(x.mat^2) - sum(diag(x.mat^2)))/(n1*(n1-1))
      Lambda.hat <- 2*(tau^2)*((1-tau)^2)*tr.Sigx2.hat
      
      Sn.std <- Sn/sqrt(Lambda.hat)
      return(Sn.std)
    }
    
    Sn1.std <- Sn_func(ind1, ind2)
    Sn2.std <- Sn_func(ind2, ind1)
    Sn.std <- (Sn1.std + Sn2.std)/sqrt(2)
    pval <- 1-pnorm(Sn.std)
  }
  
  if(!is.split) {
    # estimate only gamma
    lasso.model <- rq(y ~ z, tau=tau, method="lasso")
    mz.pred <- z.tilde%*%coefficients(lasso.model)
    
    ## construction of test statistics
    phi_func <- function(u, tau=tau) { tau - ifelse(u<0, 1, 0) }
    phi.vec <- phi_func(as.numeric(y-mz.pred), tau=tau)
    phi.mat <- outer(phi.vec, phi.vec, "*")
    x.mat <- x %*% t(x)
    Sn <- (sum(phi.mat*x.mat) - sum(diag(phi.mat*x.mat)))/n
    
    tr.Sigx2.hat <- (sum(x.mat^2) - sum(diag(x.mat^2)))/(n*(n-1))
    Lambda.hat <- 2*(tau^2)*((1-tau)^2)*tr.Sigx2.hat
    
    Sn.std <- Sn/sqrt(Lambda.hat)
    pval <- 1-pnorm(Sn.std)
  }
  return(pval)
}
