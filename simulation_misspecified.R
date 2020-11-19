library(MASS)
library(rpart)
library(partykit)
library(Rcpp)
library(pROC)
sourceCpp("tau-a.cpp")
sourceCpp("fun-tau_n.cpp")

# for this script to run
# variables n, delta, rho, tree, k need to be
# assigned values
# see quoted block codes for an example
# n <- 100 # 200, 500
# delta <- 1 # 3
# rho <- 0.2 # 0.5, 0.8
# tree <- "tree1" # "tree2", "tree3", "tree4"
# k <- 1 # 1, ... , 1000

fc <- "misfcpoly2"
print(k)
set.seed(k)

# second-order polynomial
h <- function(x, beta) { # beta is the un-transformed parameter
      out <- as.numeric(x <= 0) * (beta[1] * x^2 + beta[2] * x + beta[3]) +
            as.numeric(x > 0) * (beta[4] * x^2 + beta[5] *x + beta[3])
      return(out)
}
h1 <- function(x) {
      out <- as.numeric(x <= 0) * (-x^2/2 +x/20 - 0.5) +
            as.numeric(x > 0) * (x^2/5 + x/20 - 0.5)
      
      return(out)
}
h2 <- function(x) {
      out <- as.numeric(x <= 0) * (-x^2/10 + x/5 - 2) +
            as.numeric(x > 0) * (x^2 + x/3 - 2)
      
      return(out)
}

if (tree == "tree1") {
      # tree1: X1 > x1 & X2 > x2 & X3 > x3
      Htree <- function(x1, x2, x3) {
            n <- length(x1)
            out <- sapply(1:n, function(i) {
                  return(min(x1[i], x2[i], x3[i]))
            })
            return(out)
      }
} else if (tree == "tree2") {
      # tree2: X1 > x1 | X2 > x2 | X3 > x3
      Htree <- function(x1, x2, x3) {
            n <- length(x1)
            out <- sapply(1:n, function(i) {
                  return(max(x1[i], x2[i], x3[i]))
            })
            return(out)
      }
} else if (tree == "tree3") {
      # tree3: X1 > x1 & (X2 > x2 | X3 > x3)
      Htree <- function(x1, x2, x3) {
            n <- length(x1)
            out <- sapply(1:n, function(i) {
                  return(min(x1[i], max(x2[i], x3[i])))
            })
            return(out)
      }
} else if (tree == "tree4") {
      # tree4: (X1 > x1 & X2 > x2) | X3 > x3
      Htree <- function(x1, x2, x3) {
            n <- length(x1)
            out <- sapply(1:n, function(i) {
                  return(max(x1[i], min(x2[i], x3[i])))
            })
            return(out)
      }
} else {
      stop("Tree structure not coded.")
}

# model function
funH <- function(x, theta) { # theta is the transformed parameter
      beta <- theta
      beta[1] <- -exp(theta[1])
      beta[2] <- exp(theta[2])
      beta[4] <- exp(theta[4])
      beta[5] <- exp(theta[5])
      
      beta[6] <- -exp(theta[6])
      beta[7] <- exp(theta[7])
      beta[9] <- exp(theta[9])
      beta[10] <- exp(theta[10])
      
      
      out <- Htree(h(x[1], beta[1:5]), h(x[2], beta[6:10]), x[3])
      return(out)
}
# function to check if estimate is in parameter space
funCheck <- function(theta) { # theta is the transformed parameter
      return(TRUE)
}


# true parameter value under model
theta0 <- c(log(1/2), log(1/20), -0.5, log(1/5), log(1/20),
            log(1/10), log(1/5), -2, log(1), log(1/3))


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

# create testing data
X_temp <- mvrnorm(n * 10, mu = rep(0, 3), Sigma = tempmat)

H_temp <- apply(X_temp, 1, function(x) funH(x, theta0))

Y_temp <- sapply(H_temp, function(x) {
      return(rbinom(1, 1, pnorm(x, -2, delta)))
})

X_test <- rbind(X_temp[which(Y_temp == 0)[1:(n/2)], ],
                X_temp[which(Y_temp == 1)[1:(n/2)], ])
Y_test <- c(Y_temp[which(Y_temp == 0)[1:(n/2)]], 
            Y_temp[which(Y_temp == 1)[1:(n/2)]])


peMRC <- function(X, Y, theta_start = theta0, nstart = 50) { 
      # point estimation of theta and auc given data
      rc <- function(theta) {
            # calculates rank correlation for parameter theta
            # given data X and Y
            H <- apply(X, 1, function(x) funH(x, theta))
            result <- taua(H, Y)
            return(result)
      }
      p <- length(theta0)
      # compute estimate
      start <- matrix(mvrnorm(n = nstart, mu = theta_start, Sigma = diag(1, p)), 
                      ncol = p)
      
      start <- rbind(start, theta_start)
      comp.sol <- do.call(rbind, lapply(1:nrow(start), function(i) {
            fit <- optim(start[i, ], rc, control = list(fnscale = -1))
            temp <- c(fit$par, fit$value)
      }))
      
      check <- apply(comp.sol, 1, funCheck) # only take the first p elements
      
      comp.sol <- comp.sol[check, ]
      
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

H0_train <- apply(X, 1, function(x) funH(x, theta0)) # under true parameter
H_pe_train <- apply(X, 1, function(x) funH(x, pe))
H_pe_test <- apply(X_test, 1, function(x) funH(x, pe))




# load estimation function
pwl <- function(m, theta, knot = NULL) { 
      # piece-wise linear function
      # if (length(theta) != length(knot) + 2) {
      #       stop("length of parameter and number of knots given are not compatible")
      # }
      theta[-1] <- exp(theta[-1])
      if (is.null(knot)) {
            mt <- c(1, m)
      } else {
            knot <- sort(knot)
            mt <- unlist(c(1, m, sapply(knot, function(mk) {
                  out <- (m - mk) * as.numeric(m >= mk)
                  return(out)
            })))
      }
      return(sum(mt * theta))
}

rc_est_cart <- function(var1, var2, var3, outcome, nk) {
      rc <- function(theta1, theta2, nknot = NULL, knot = NULL) {
            if (nk != 0) {
                  lower1 <- range(var1)[1]
                  span1 <- range(var1)[2] - range(var1)[1]
                  knot1 <- lower1 + 1:nk / (nk + 1) * span1
                  
                  knot1r <- -rev(knot1)
                  
                  lower2 <- range(var2)[1]
                  span2 <- range(var2)[2] - range(var2)[1]
                  knot2 <- lower2 + 1:nk / (nk + 1) * span2
                  
                  knot2r <- -rev(knot2)
            } else {
                  knot1 <- NULL
                  knot2 <- NULL
                  knot1r <- NULL
                  knot2r <- NULL
            }
            
            H <- as.numeric(sapply(1:length(var1), function(i) {
                  out <- Htree(pwl(var1[i], theta1, knot1), 
                               pwl(var2[i], theta2, knot2), var3)
                  return(out)
            }))
            result <- taua(H, outcome)
      }
      
      p <- nk + 2
      optimfn <- function(theta) {
            return(rc(theta1 = theta[1:p], 
                      theta2 = theta[1:p+p],
                      nknot = nk))
      }
      
      start <- matrix(mvrnorm(n = 500, mu = rep(0, 2*p), 
                              Sigma = diag(2, 2*p)), ncol = 2*p)
      comp.sol <- do.call(rbind, lapply(1:nrow(start), function(i) {
            fit <- optim(start[i, ], optimfn, control = list(fnscale = -1))
            temp <- c(fit$par, fit$value)
      }))
      sol <- matrix(comp.sol[comp.sol[, 2*p+1] == max(comp.sol[, 2*p+1]), 1:(2*p)], ncol = 2*p)[1, ]
      return(sol)
}

rcfit_wcv_cart <- function(var1, var2, var3, outcome,
                           nkvec = 0:2, nfold = 5, seed = 37) {
      set.seed(seed)
      folds0 <- createFolds(which(outcome == 0), nfold)
      folds1 <- createFolds(which(outcome == 1), nfold)
      folds <- list()
      for (l in 1:nfold) {
            folds[[l]] <- c(which(outcome == 0)[folds0[[l]]], 
                            which(outcome == 1)[folds1[[l]]])
      }
      cvrc <- NULL
      for (nk in nkvec) {
            print(nk)
            p <- nk + 2
            if (nk != 0) {
                  lower1 <- range(var1)[1]
                  span1 <- range(var1)[2] - range(var1)[1]
                  knot1 <- lower1 + 1:nk/(nk+1)*span1
                  knot1r <- -rev(knot1)
                  
                  lower2 <- range(var2)[1]
                  span2 <- range(var2)[2] - range(var2)[1]
                  knot2 <- lower2 + 1:nk/(nk+1)*span2
                  
                  knot2r <- -rev(knot2)
            } else {
                  knot1 <- NULL
                  knot2 <- NULL
                  knot1r <- NULL
                  knot2r <- NULL
            }
            tempauc <- NULL
            for (i in 1:nfold) {
                  # print(i)
                  est <- rc_est_cart(var1 = var1[-folds[[i]]], 
                                     var2 = var2[-folds[[i]]], 
                                     var3 = var3[-folds[[i]]],
                                     outcome = outcome[-folds[[i]]], nk)
                  
                  H <- as.numeric(sapply(folds[[i]], function(i) {
                        out <- Htree(pwl(var1[i], est[1:p], knot1), 
                                     pwl(var2[i], est[1:p+p], knot2), var3)
                        return(out)
                  }))
                  
                  tempauc <- c(tempauc, taua(H, outcome[folds[[i]]]))
                  # evarauc <- aucvar(H, dat$prog[folds[[i]]])
            }
            print(tempauc)
            print(mean(tempauc))
            cvrc <- rbind(cvrc, c(nk, mean(tempauc)))
      }
      return(cvrc)
}


# fit tree with linear knot = 0, 1, 2, 
# or through cross-validation

# nk = 0
fit_nk0 <- rc_est_cart(X[, 1], X[, 2], X[, 3], Y, nk = 0)
nk <- 0
p <- nk+2
if (nk != 0) {
      lower1 <- range(X[, 1])[1]
      span1 <- range(X[, 1])[2] - range(X[, 1])[1]
      knot1 <- lower1 + 1:nk/(nk+1)*span1
      knot1r <- -rev(knot1)
      
      lower2 <- range(X[, 2])[1]
      span2 <- range(X[, 2])[2] - range(X[, 2])[1]
      knot2 <- lower2 + 1:nk/(nk+1)*span2
      knot2r <- -rev(knot2)
} else {
      knot1 <- NULL
      knot2 <- NULL
      knot1r <- NULL
      knot2r <- NULL
}
H_train_nk0 <- sapply(1:n, function(i) {
      out <- Htree(pwl(X[i, 1], fit_nk0[1:p], knot1), 
                   pwl(X[i, 2], fit_nk0[1:p+p], knot2), var3)
})
H_test_nk0 <- sapply(1:n, function(i) {
      out <- Htree(pwl(X_test[i, 1], fit_nk0[1:p], knot1), 
                   pwl(X_test[i, 2], fit_nk0[1:p+p], knot2), var3)
})



# nk = 1
fit_nk1 <- rc_est_cart(X[, 1], X[, 2], X[, 3], Y, nk = 1)
nk <- 1
p <- nk+2
if (nk != 0) {
      lower1 <- range(X[, 1])[1]
      span1 <- range(X[, 1])[2] - range(X[, 1])[1]
      knot1 <- lower1 + 1:nk/(nk+1)*span1
      knot1r <- -rev(knot1)
      
      lower2 <- range(X[, 2])[1]
      span2 <- range(X[, 2])[2] - range(X[, 2])[1]
      knot2 <- lower2 + 1:nk/(nk+1)*span2
      knot2r <- -rev(knot2)
} else {
      knot1 <- NULL
      knot2 <- NULL
      knot1r <- NULL
      knot2r <- NULL
}
H_train_nk1 <- sapply(1:n, function(i) {
      out <- Htree(pwl(X[i, 1], fit_nk1[1:p], knot1), 
                   pwl(X[i, 2], fit_nk1[1:p+p], knot2), var3)
})
H_test_nk1 <- sapply(1:n, function(i) {
      out <- Htree(pwl(X_test[i, 1], fit_nk1[1:p], knot1), 
                   pwl(X_test[i, 2], fit_nk1[1:p+p], knot2), var3)
})


# nk = 2
fit_nk2 <- rc_est_cart(X[, 1], X[, 2], X[, 3], Y, nk = 2)
nk <- 2
p <- nk+2
if (nk != 0) {
      lower1 <- range(X[, 1])[1]
      span1 <- range(X[, 1])[2] - range(X[, 1])[1]
      knot1 <- lower1 + 1:nk/(nk+1)*span1
      knot1r <- -rev(knot1)
      
      lower2 <- range(X[, 2])[1]
      span2 <- range(X[, 2])[2] - range(X[, 2])[1]
      knot2 <- lower2 + 1:nk/(nk+1)*span2
      knot2r <- -rev(knot2)
} else {
      knot1 <- NULL
      knot2 <- NULL
      knot1r <- NULL
      knot2r <- NULL
}
H_train_nk2 <- sapply(1:n, function(i) {
      out <- Htree(pwl(X[i, 1], fit_nk2[1:p], knot1), 
                   pwl(X[i, 2], fit_nk2[1:p+p], knot2), var3)
})
H_test_nk2 <- sapply(1:n, function(i) {
      out <- Htree(pwl(X_test[i, 1], fit_nk2[1:p], knot1), 
                   pwl(X_test[i, 2], fit_nk2[1:p+p], knot2), var3)
})

main.dir <- getwd()
sub.dir <- paste0("mrcsim-", fc, tree,"-r", rho*10, "-d", delta, "-n", n)
if(file.exists(sub.dir)) {
      setwd(file.path(main.dir, sub.dir))
} else {
      dir.create(file.path(main.dir, sub.dir))
      setwd(file.path(main.dir, sub.dir))
}
save.image(file = paste0(sub.dir, "-k", k, ".rda"))
setwd(file.path(main.dir))