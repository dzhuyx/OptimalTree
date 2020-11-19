library(MASS)
library(rpart)
library(partykit)
library(Rcpp)
library(pROC)
sourceCpp("tau-a.cpp")
sourceCpp("fun-tau_n.cpp")

# for this script to run
# variables n, delta, rho, tree, fc, k need to be
# assigned values
# see quoted block codes for an example
# n <- 100 # 200, 500
# delta <- 1 # 3
# rho <- 0.2 # 0.5, 0.8
# tree <- "tree1" # "tree2", "tree3", "tree4"
# fc <- "linear0" # "linear1", "linear2", "poly2"
# k <- 1 # 1, ... , 1000
print(k)
set.seed(k)

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

if (fc == "linear0") {
      # linear h with 0 knot
      h <- function(x, beta) { # beta is the un-transformed parameter
            out <- beta[1] * x + beta[2]
            return(out)
      }
      h1 <- function(x) {
            return(x - 1)
      }
      h2 <- function(x) {
            return(2 * x + 0.5)
      }
      # model function
      funH <- function(x, theta) { # theta is the transformed parameter
            beta <- theta
            beta[1] <- exp(theta[1])
            beta[3] <- exp(theta[3])
            
            out <- Htree(h(x[1], beta[1:2]), h(x[2], beta[3:4]), x[3])
            return(out)
      }
      # function to check if estimate is in parameter space
      funCheck <- function(theta) { # theta is the transformed parameter
            return(TRUE)
      }
      
      
      # true parameter value under model
      theta0 <- c(log(1), -1, log(2), 0.5)
} else if (fc == "linear1") {
      # linear h with 1 knot
      h <- function(x, beta) { # beta is the un-transformed parameter
            out <- as.numeric(x <= 0) * (beta[1] * x + beta[2]) + 
                  as.numeric(x > 0) * (beta[3] * x + beta[2])
            return(out)
      }
      # linear h with 1 knot (at median 0)
      h1 <- function(x) {
            out <- as.numeric(x <= 0) * (x - 1) + 
                  as.numeric(x > 0) * (3 * x - 1)
            return(out)
      }
      h2 <- function(x) {
            out <- as.numeric(x <= 0) * (2 * x + 0.5) + 
                  as.numeric(x > 0) * (0.5 * x + 0.5)
            return(out)
      }
      # model function
      funH <- function(x, theta) { # theta is the transformed parameter
            beta <- theta
            beta[1] <- exp(theta[1])
            beta[4] <- exp(theta[4])
            
            out <- Htree(h(x[1], beta[1:3]), h(x[2], beta[4:6]), x[3])
            return(out)
      }
      # function to check if estimate is in parameter space
      funCheck <- function(theta) { # theta is the transformed parameter
            return(TRUE)
      }
      
      
      # true parameter value under model
      theta0 <- c(log(1), -1, log(3),
                  log(2), 0.5, log(0.5))
} else if (fc == "linear2") {
      # linear h with 2 knots
      h <- function(x, beta) { # beta is the un-transformed parameter
            out <- as.numeric(x <= -1) * (beta[1] * x + beta[2]) + 
                  as.numeric(x > -1 & x <= 1) * (beta[3] * x + beta[2] + beta[3] - beta[1]) + 
                  as.numeric(x > 1) * (beta[4] * x + 2 * beta[3] + beta[2] - beta[1] - beta[4])
            return(out)
      }
      # linear h with 2 knots
      h1 <- function(x) {
            out <- as.numeric(x <= -1) * (0.5 * x - 1.5) + 
                  as.numeric(x > -1 & x <= 1) * (x - 1) + 
                  as.numeric(x > 1) * (3 * x - 3)
            return(out)
      }
      h2 <- function(x) {
            out <- as.numeric(x <= -1) * (x - 0.5) + 
                  as.numeric(x > -1 & x <= 1) * (2 * x + 0.5) + 
                  as.numeric(x > 1) * (x/2 + 2)
      }
      # model function
      funH <- function(x, theta) { # theta is the transformed parameter
            beta <- theta
            beta[1] <- exp(theta[1])
            beta[3] <- exp(theta[3])
            beta[4] <- exp(theta[4])
            
            out <- Htree(h(x[1], beta[1:4]), h(x[2], beta[5:8]), x[3])
            return(out)
      }
      # function to check if estimate is in parameter space
      funCheck <- function(theta) { # theta is the transformed parameter
            return(TRUE)
      }
      
      
      # true parameter value under model
      theta0 <- c(log(0.5), -1.5, log(1), log(3),
                  log(1), -0.5, log(2), log(0.5))
} else if (fc == "poly2") {
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
} else {
      stop("Function class not coded.")
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

# create training and testing data.frame
dat_train <- data.frame(Y = Y, X1 = X[, 1], X2 = X[, 2], X3 = X[, 3])
dat_test <- data.frame(Y = Y_test, X1 = X_test[, 1],
                       X2 = X_test[, 2], X3 = X_test[, 3])

# fit CART and CIT
set.seed(37)
fit_cart <- rpart(Y ~ X1 + X2 + X3, data = dat_train, method = "class",
                  control = rpart.control(xval = 5)) # 5-fold cross-validation
cptab_cart <- printcp(fit_cart)
fit_cart_pruned <- prune(fit_cart, 
                         cp = cptab_cart[cptab_cart[, 4] == min(cptab_cart[, 4]), 1][1])
fit_cit <- ctree(Y ~ X1 + X2 + X3, data = dat_train)

H_cart_train <- predict(fit_cart_pruned, newdata = dat_train)[, 1]
H_cit_train <- predict(fit_cit, newdata = dat_train)

H_cart_test <- predict(fit_cart_pruned, newdata = dat_test)[, 1]
H_cit_test <- predict(fit_cit, newdata = dat_test)

# calculate ROCs
roc_pe_train <- roc(factor(dat_train$Y), H_pe_train)
roc_cart_train <- roc(factor(dat_train$Y), H_cart_train)
roc_cit_train <- roc(factor(dat_train$Y), H_cit_train)

roc_pe_test <- roc(factor(dat_test$Y), H_pe_test)
roc_cart_test <- roc(factor(dat_test$Y), H_cart_test)
roc_cit_test <- roc(factor(dat_test$Y), H_cit_test)

# training set AUC
auc_pe_train <- auc(roc_pe_train)[1]
auc_cart_train <- auc(roc_cart_train)[1]
auc_cit_train <- auc(roc_cit_train)[1]
# testing set AUC
auc_pe_test <- auc(roc_pe_test)[1]
auc_cart_test <- auc(roc_cart_test)[1]
auc_cit_test <- auc(roc_cit_test)[1]

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
