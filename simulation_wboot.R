library(MASS)

library(Rcpp)
sourceCpp("tau-a.cpp")
sourceCpp("fun-tau_n.cpp")


# for this script to run
# variables n, delta, rho, k need to be
# assigned values
# see quoted block codes for an example
# n <- 50
# delta <- 1
# rho <- 0.2
# k <- 1

OneSim <- function(loop) {
    print(loop)
set.seed(loop)
# specify model
theta0 <- c(log(1), -1, log(2), 0.5)
p <- length(theta0)
funH <- function(x, theta) {
      theta1 <- theta[1]
      theta2 <- theta[2]
      theta3 <- theta[3]
      theta4 <- theta[4]
      return(min(x[1], theta2 + exp(theta1) * x[2], theta4 + exp(theta3) * x[3]))
}

# generate training sample
tempmat <- matrix(10 * rho, 3, 3)
diag(tempmat) <- 10
X <- mvrnorm(n * 10, mu = rep(0, 3), Sigma = tempmat)

H <- apply(X, 1, function(x) funH(x, theta0))

Y <- sapply(H, function(x) {
      return(rbinom(1, 1, pnorm(x, -2, delta)))
})

X <- rbind(X[which(Y == 0)[1:(n/2)], ],
           X[which(Y == 1)[1:(n/2)], ])
Y <- c(Y[which(Y == 0)[1:(n/2)]], Y[which(Y == 1)[1:(n/2)]])


peMRC <- function(X, Y, theta_start = theta0) { 
      # point estimation of theta and auc given data
      rc <- function(theta) {
            # calculates rank correlation for parameter theta
            # given data X and Y
            H <- apply(X, 1, function(x) funH(x, theta))
            result <- taua(H, Y)
            return(result)
      }
      
      # compute estimate
      start <- matrix(mvrnorm(n = 10, mu = theta_start, Sigma = diag(1, p)), 
                      ncol = p)
      start <- rbind(start, theta_start)
      comp.sol <- do.call(rbind, lapply(1:nrow(start), function(i) {
            fit <- optim(start[i, ], rc, control = list(fnscale = -1))
            temp <- c(fit$par, fit$value)
      }))
      
    sol <- comp.sol[comp.sol[, p+1] == max(comp.sol[, p+1]), 1:p, drop = F]
      if (nrow(sol) > 1) {
        dis <- apply(sol, 1, function(x) sum((x - theta0)^2))
        out <- sol[dis == min(dis), , drop = F][1, ]
      } else {
        out <- sol
      }
      
      theta_est <- out[1:p]  # estimate
      return(theta_est)
}

pe <- peMRC(X, Y)

# boostrap to compute variance
bootn <- 10000  # set number of bootstrapped samples to sample size
smp <- cbind(matrix(sample(1:(n/2), bootn * n/2, replace = T), nrow = bootn), 
             matrix(sample(1:(n/2) + n/2, bootn*n/2, replace = T), nrow = bootn))
pe_boot <- lapply(1:nrow(smp), function(x) {
      temp <- peMRC(X[smp[x, ], ], Y[smp[x, ]], pe)
      return(temp)
})
return(list(pe = pe, boot = pe_boot)) 
}

# array version
res <- OneSim(k)

main.dir <- getwd()
sub.dir <- paste0("mrc-linear7-d", delta, "-r", rho, "-n", n)
if(file.exists(sub.dir)) {
    setwd(file.path(main.dir, sub.dir))
} else {
    dir.create(file.path(main.dir, sub.dir))
    setwd(file.path(main.dir, sub.dir))
}
save.image(file = paste0(sub.dir, "-k", k, ".rda"))
setwd(file.path(main.dir))
