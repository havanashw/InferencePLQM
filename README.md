# InferencePLQM
R code for "Inference for partially linear quantile regression models in ultrahigh dimension".

## Usage

```{R}
InferencePLQM(x, z, y, tau, is.split=FALSE, is.screen=FALSE, seed.fix)
```

## Required Packages
- `quantreg`
- `MASS`

## Inputs
- `x`: A matrix of n*d, where n is the sample size and d is the dimension of interested covariates.
- `z`: A matrix of n*q, where q is the dimension of potential nuisance variables.
- `y`: The response variable, with length n.
- `tau`: The quantile level.
- `is.split`: A logical value indicating whether to perform data splitting for testing.
- `is.screen`: A logical value indicating whether to apply the CQ-SIS screening procedure implemented by function `CQSIS`, where `cn` is the screening number
- `seed.fix`: A fixed seed value for reproducibility when data splitting is used.

## Examples
```{R}
library(quantreg)
library(parallel)
source("InferencePLQM.R")

# setting
tau <- 0.75; h0 <- 0
n <- 300; p <- 1000; d <- 500; q <- p - d
varsigma <- 0; sx <- 10
sz.signal <- 0.5; sz <- 20
Sig <- toeplitz(0.5^seq(0, p-1))

# generate x, z and y
v <- mvrnorm(n, mu=rep(0, p), Sigma=Sig)
x <- v[,1:d]; z <- v[,-(1:d)]
x.tilde <- cbind(rep(1, n), x)
b0.beta0 <- c(1, rep(varsigma/sqrt(sx), sx), rep(0, d-sx))
gamma0 <- c(rep(sz.signal, sz), rep(0, q-sz))
error <- rnorm(n) - qnorm(tau)
mz <- z%*%gamma0 
y <- x.tilde%*%b0.beta0 + mz + (1+h0*x[,1])*error

# p-value
InferencePLQM(x, z, y, tau, is.split=TRUE, is.screen=TRUE, seed.fix=1)
```
