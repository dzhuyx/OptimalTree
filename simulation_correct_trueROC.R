library(MASS)
library(pROC)
n <- 1000000
dat_trueROC <- data.frame()
for (fc in c("linear0", "linear1", "linear2", "poly2")) {
      for (tree in paste0("tree", 1:4)) {
            for (rho in c(0.2, 0.5, 0.8)) {
                  for (delta in c(1, 3)) {
                        Sigma <- diag(1, 3) * 10 * (1 - rho) + matrix(10 * rho, 3, 3)
                        
                        if (fc == "linear0") {
                              # linear h with 0 knot
                              h1 <- function(x) {
                                    return(x - 1)
                              }
                              h2 <- function(x) {
                                    return(2 * x + 0.5)
                              }
                        } else if (fc == "linear1") {
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
                        } else if (fc == "linear2") {
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
                        } else {
                              # second-order polynomial with 1 knot (at median 0)
                              h1 <- function(x) {
                                    out <- as.numeric(x > 0) * (x^2/5 + x/20 - 0.5) + 
                                          as.numeric(x <= 0) * (-x^2/2 +x/20 - 0.5)
                                    return(out)
                              }
                              h2 <- function(x) {
                                    out <- as.numeric(x > 0) * (x^2 + x/3 - 2) + 
                                          as.numeric(x <= 0) * (-x^2/10 + x/5 - 2)
                                    return(out)
                              }
                        }
                        
                        if (tree == "tree1") {
                              # tree1: X1 > x1 & X2 > x2 & X3 > x3
                              H <- function(x1, x2, x3) {
                                    n <- length(x1)
                                    out <- sapply(1:n, function(i) {
                                          return(min(x1[i], x2[i], x3[i]))
                                    })
                                    return(out)
                              }
                        } else if (tree == "tree2") {
                              # tree2: X1 > x1 | X2 > x2 | X3 > x3 
                              H <- function(x1, x2, x3) {
                                    n <- length(x1)
                                    out <- sapply(1:n, function(i) {
                                          return(max(x1[i], x2[i], x3[i]))
                                    })
                                    return(out)
                              }
                        } else if (tree == "tree3") {
                              # tree3: X1 > x1 & (X2 > x2 | X3 > x3)
                              H <- function(x1, x2, x3) {
                                    n <- length(x1)
                                    out <- sapply(1:n, function(i) {
                                          return(min(x1[i], max(x2[i], x3[i])))
                                    })
                                    return(out)
                              }
                        } else {
                              # tree4: X1 > x1 | (X2 > x2 & X3 > x3)
                              H <- function(x1, x2, x3) {
                                    n <- length(x1)
                                    out <- sapply(1:n, function(i) {
                                          return(max(x1[i], min(x2[i], x3[i])))
                                    })
                                    return(out)
                              }
                        }
                        
                        # calculate real ROC
                        set.seed(37)
                        X <- mvrnorm(n*100, mu = rep(0, 3), Sigma = Sigma)
                        Hvec <- H(h1(X[, 1]), h2(X[, 2]), X[, 3])
                        U <- rnorm(n, 2, delta)
                        D <- as.numeric(Hvec + U > 0)
                        sampleid0 <- sample(which(D == 0), n, replace = F)
                        sampleid1 <- sample(which(D == 1), n, replace = F)
                        Hvec <- Hvec[c(sampleid0, sampleid1)]
                        D <- D[c(sampleid0, sampleid1)]
                        X <- X[c(sampleid0, sampleid1), ]
                        sum(h1(X[, 1]) == Hvec)
                        sum(h2(X[, 2]) == Hvec)
                        sum(X[, 3] == Hvec)
                        roc_true <- roc(factor(D), Hvec)
                        p <- length(roc_true$sensitivities)
                        dat_trueROC <- rbind(dat_trueROC,
                                             data.frame(fc = rep(fc, p), 
                                                        tree = rep(tree, p), 
                                                        rho = rep(rho, p), 
                                                        delta = rep(delta, p),
                                                        sens = roc_true$sensitivities,
                                                        spec = roc_true$specificities,
                                                        auc = rep(auc(roc_true), p)))
                  }
            }
      }
}

save(list = "dat_trueROC", file = "trueROC.rda")