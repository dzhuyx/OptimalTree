rm(list = ls())
library(pROC)
library(caret)
library(Rcpp)
library(MASS)
library(beepr)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(partykit)
sourceCpp("tau-a.cpp")
# read data
dat <- read.csv(file = "prostate_cancer_challenge_data_training/CoreTable_training.csv", 
                na.strings = c(".", " "))

# define outcome variable 
dat$event <- NA
dat$event[dat$LKADT_P < 500 & dat$DEATH == "YES"] <- 1
dat$event[dat$LKADT_P >= 500] <- 0
dat <- dat[!is.na(dat$event), ]

# convert variables to correct type
# code NA values
dat$AGEGRP <- as.character(dat$AGEGRP)
dat$AGEGRP[dat$AGEGRP == ">=85"] <- 85
dat$AGEGRP <- as.numeric(dat$AGEGRP)
dat$AGEGRP2 <- as.character(dat$AGEGRP2)
dat$AGEGRP2[dat$AGEGRP2 == ">=75"] <- 3
dat$AGEGRP2[dat$AGEGRP2 == "18-64"] <- 1
dat$AGEGRP2[dat$AGEGRP2 == "65-74"] <- 2

dato <- dat

# estimate AUCs of variables 
auc_analysis <- matrix(NA, nrow = 120, ncol = 2)
for (j in 12:131) {
      nobs <- sum(!is.na(dat[, j]))
      if (is.numeric(dat[, j]) & nobs > 900) {
            temp <- roc(dat$event, dat[, j])
            auc_est <- paste0(ci.auc(temp), collapse = ",")
            if (auc_est < 0.5) {
                  auc_est <- paste0(ci.auc(temp), collapse = ",")
                  dat[, j] <- -dat[, j]
            }
      } else {
            auc_est <- NA
      }
      
      auc_analysis[j-11, ] <- c(nobs, auc_est)
}
row.names(auc_analysis) <- colnames(dat)[12:131]

auc_analysis <- auc_analysis[!is.na(auc_analysis[, 2]) & auc_analysis[, 1] > 900, ]
auc_analysis <- auc_analysis[order(auc_analysis[, 2], decreasing = T), ]
View(auc_analysis)

marker_cov_pearson <- matrix(NA, nrow = nrow(auc_analysis), ncol = nrow(auc_analysis))
marker_cov_kendall <- matrix(NA, nrow = nrow(auc_analysis), ncol = nrow(auc_analysis))
for (i in 1:nrow(auc_analysis)) {
      for (j in 1:nrow(auc_analysis)) {
            marker_cov_pearson[i, j] <- cor.test(dat[, rownames(auc_analysis)[i]], 
                                                 dat[, rownames(auc_analysis)[j]],
                                                 method = "pearson")$estimate
            marker_cov_kendall[i, j] <- cor.test(dat[, rownames(auc_analysis)[i]], 
                                                 dat[, rownames(auc_analysis)[j]],
                                                 method = "kendall")$estimate
      }
}
rownames(marker_cov_pearson) <- rownames(auc_analysis)
colnames(marker_cov_pearson) <- rownames(auc_analysis)
rownames(marker_cov_kendall) <- rownames(auc_analysis)
colnames(marker_cov_kendall) <- rownames(auc_analysis)

# generic functions for estimation
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

# functions for "and" tree
rc_est_var1andvar2 <- function(var1, var2, outcome, nk) {
      rc <- function(theta1, nknot = NULL, knot = NULL) {
            if (nk != 0) {
                  lower <- range(var2)[1]
                  span <- range(var2)[2] - range(var2)[1]
                  knot <- lower + 1:nk / (nk + 1) * span
            } else {
                  knot <- NULL
            }
            
            H <- as.numeric(sapply(1:length(var1), function(i) {
                  out <- min(pwl(as.numeric(var2[i]), theta1, knot = knot), 
                             as.numeric(var1[i]))
                  return(out)
            }))
            result <- taua(H, outcome)
      }
      
      p <- nk + 2
      optimfn <- function(theta) {
            return(rc(theta, nknot = nk))
      }
      
      start <- matrix(mvrnorm(n = 500, mu = rep(0, p), 
                              Sigma = diag(2, p)), ncol = p)
      comp.sol <- do.call(rbind, lapply(1:nrow(start), function(i) {
            fit <- optim(start[i, ], optimfn, control = list(fnscale = -1))
            temp <- c(fit$par, fit$value)
      }))
      sol <- matrix(comp.sol[comp.sol[, p+1] == max(comp.sol[, p+1]), 1:p], ncol = p)[1, ]
      return(sol)
}
# k-fold cross validation
rcfit_wcv_var1andvar2 <- function(var1_name, var2_name, outcome_name, dat, 
                                  nkvec = 3, nfold = 5, seed = 37) {
      set.seed(seed)
      folds0 <- createFolds(which(dat[, outcome_name] == 0), nfold)
      folds1 <- createFolds(which(dat[, outcome_name] == 1), nfold)
      folds <- list()
      for (l in 1:nfold) {
            folds[[l]] <- c(which(dat[, outcome_name] == 0)[folds0[[l]]], 
                            which(dat[, outcome_name] == 1)[folds1[[l]]])
      }
      cvrc <- NULL
      for (nk in 0:nkvec) {
            print(nk)
            p <- nk + 2
            if (nk != 0) {
                  lower <- range(dat[, var2_name])[1]
                  span <- range(dat[, var2_name])[2] - range(dat[, var2_name])[1]
                  knot <- lower + 1:nk/(nk+1)*span
            } else {
                  knot <- NULL
            }
            tempauc <- NULL
            for (i in 1:nfold) {
                  # print(i)
                  est <- rc_est_var1andvar2(var1 = dat[-folds[[i]], var1_name], 
                                            var2 = dat[-folds[[i]], var2_name], 
                                            outcome = dat[-folds[[i]], outcome_name], nk)
                  H <- as.numeric(apply(dat[folds[[i]], ], 1, function(x) {
                        out <- min(pwl(as.numeric(x[var2_name]), est, knot = knot), 
                                   as.numeric(x[var1_name]))
                        return(out)
                  }))
                  tempauc <- c(tempauc, taua(H, dat[, outcome_name][folds[[i]]]))
                  # evarauc <- aucvar(H, dat$prog[folds[[i]]])
            }
            print(tempauc)
            print(mean(tempauc))
            cvrc <- rbind(cvrc, c(nk, mean(tempauc)))
      }
      return(cvrc)
}

# functions for "or" tree
rc_est_var1orvar2 <- function(var1, var2, outcome, nk) {
      rc <- function(theta1, nknot = NULL, knot = NULL) {
            if (nk != 0) {
                  lower <- range(var2)[1]
                  span <- range(var2)[2] - range(var2)[1]
                  knot <- lower + 1:nk / (nk + 1) * span
            } else {
                  knot <- NULL
            }
            
            H <- as.numeric(sapply(1:length(var1), function(i) {
                  out <- max(pwl(as.numeric(var2[i]), theta1, knot = knot), 
                             as.numeric(var1[i]))
                  return(out)
            }))
            result <- taua(H, outcome)
      }
      
      p <- nk + 2
      optimfn <- function(theta) {
            return(rc(theta, nknot = nk))
      }
      
      start <- matrix(mvrnorm(n = 500, mu = rep(0, p), 
                              Sigma = diag(2, p)), ncol = p)
      comp.sol <- do.call(rbind, lapply(1:nrow(start), function(i) {
            fit <- optim(start[i, ], optimfn, control = list(fnscale = -1))
            temp <- c(fit$par, fit$value)
      }))
      sol <- matrix(comp.sol[comp.sol[, p+1] == max(comp.sol[, p+1]), 1:p], ncol = p)[1, ]
      return(sol)
}

# k-fold cross validation
rcfit_wcv_var1orvar2 <- function(var1_name, var2_name, outcome_name, dat, 
                                  nkvec = 3, nfold = 5, seed = 37) {
      set.seed(seed)
      folds0 <- createFolds(which(dat[, outcome_name] == 0), nfold)
      folds1 <- createFolds(which(dat[, outcome_name] == 1), nfold)
      folds <- list()
      for (l in 1:nfold) {
            folds[[l]] <- c(which(dat[, outcome_name] == 0)[folds0[[l]]], 
                            which(dat[, outcome_name] == 1)[folds1[[l]]])
      }
      cvrc <- NULL
      for (nk in 0:nkvec) {
            print(nk)
            p <- nk + 2
            if (nk != 0) {
                  lower <- range(dat[, var2_name])[1]
                  span <- range(dat[, var2_name])[2] - range(dat[, var2_name])[1]
                  knot <- lower + 1:nk/(nk+1)*span
            } else {
                  knot <- NULL
            }
            tempauc <- NULL
            for (i in 1:nfold) {
                  # print(i)
                  est <- rc_est_var1orvar2(var1 = dat[-folds[[i]], var1_name], 
                                            var2 = dat[-folds[[i]], var2_name], 
                                            outcome = dat[-folds[[i]], outcome_name], nk)
                  H <- as.numeric(apply(dat[folds[[i]], ], 1, function(x) {
                        out <- max(pwl(as.numeric(x[var2_name]), est, knot = knot), 
                                   as.numeric(x[var1_name]))
                        return(out)
                  }))
                  tempauc <- c(tempauc, taua(H, dat[, outcome_name][folds[[i]]]))
                  # evarauc <- aucvar(H, dat$prog[folds[[i]]])
            }
            print(tempauc)
            print(mean(tempauc))
            cvrc <- rbind(cvrc, c(nk, mean(tempauc)))
      }
      return(cvrc)
}


# analysis
# choose top 8 as markers and use complete data
marker_name  <- rownames(auc_analysis[order(auc_analysis[, 2], decreasing = T), ])[1:8]
dat <- dat[complete.cases(dat[, c(marker_name, "event")]), ]
dato <- dato[complete.cases(dato[, c(marker_name, "event")]), ]
# standardize variables to use
# record mean and sd
marker_mean <- colMeans(dat[, marker_name])
marker_sd <- apply(dat[, marker_name], 2, sd)
dat[, marker_name] <- scale(dat[, marker_name])

# split dat into training and testing sets
set.seed(37)
folds_control <- createFolds(which(dat$event == 0), 2)
folds_case <- createFolds(which(dat$event == 1), 2)
controlid <- which(dat$event == 0)
caseid <- which(dat$event == 1)
dat_train <- dat[c(controlid[folds_control[[1]]], caseid[folds_case[[1]]]), ]
dat_test <- dat[c(controlid[folds_control[[2]]], caseid[folds_case[[2]]]), ]
dato_train <- dato[c(controlid[folds_control[[1]]], caseid[folds_case[[1]]]), ]
dato_test <- dato[c(controlid[folds_control[[2]]], caseid[folds_case[[2]]]), ]

# "and" and "or" tree
out_andtree <- list()
out_ortree <- list()
count <- 0
for (i in 1:5) {
      for (j in 1:5) {
            if (i < j) {
                  print(paste(marker_name[i], marker_name[j]))
                  count <- count + 1
                  
                  fit_cv1 <- rcfit_wcv_var1andvar2(var1_name = marker_name[i],
                                                   var2_name = marker_name[j], 
                                                   outcome_name = "event", 
                                                   dat = dat_train)
                  nk <- fit_cv1[fit_cv1[, 2] == max(fit_cv1[, 2]), 1][1]
                  p <- nk+2
                  fit1 <- rc_est_var1andvar2(var1 = dat_train[, marker_name[i]], 
                                             var2 = dat_train[, marker_name[j]],
                                             outcome = dat_train$event, nk)
                  if (nk != 0) {
                        lower <- range(dat_train[, marker_name[j]])[1]
                        span <- range(dat_train[, marker_name[j]])[2] - 
                              range(dat_train[, marker_name[j]])[1]
                        knot <- lower + 1:nk/(nk+1)*span
                  } else {
                        knot <- NULL
                        
                  }
                  H_est1 <- as.numeric(apply(dat_test, 1, function(x) {
                        out <- min(pwl(as.numeric(x[marker_name[j]]), 
                                       fit1, knot = knot), as.numeric(x[marker_name[i]]))
                        return(out)
                  }))
                  eauc <- taua(H_est1, dat_test$event)
                  
                  H_est1_train <- as.numeric(apply(dat_train, 1, function(x) {
                        out <- min(pwl(as.numeric(x[marker_name[j]]), 
                                       fit1, knot = knot), as.numeric(x[marker_name[i]]))
                        return(out)
                  }))
                  eauc_train <- taua(H_est1_train, dat_train$event)
                  
                  out_andtree[[count]] <- list(var1 = marker_name[i],
                                               var2 = marker_name[j],
                                               fit_cv = fit_cv1,
                                               nk = nk, fit = fit1,
                                               eauc_test = eauc, 
                                               eauc_train = eauc_train)
                  
                  
                  
                  # or tree
                  fit_cv2 <- rcfit_wcv_var1orvar2(var1_name = marker_name[i],
                                                  var2_name = marker_name[j], 
                                                  outcome_name = "event", dat = dat_train)
                  nk <- fit_cv2[fit_cv2[, 2] == max(fit_cv2[, 2]), 1][1]
                  p <- nk + 2
                  fit2 <- rc_est_var1orvar2(var1 = dat_train[, marker_name[i]], 
                                            var2 = dat_train[, marker_name[j]],
                                            outcome = dat_train$event, nk)
                  if (nk != 0) {
                        lower <- range(dat_train[, marker_name[j]])[1]
                        span <- range(dat_train[, marker_name[j]])[2] - 
                              range(dat_train[, marker_name[j]])[1]
                        knot <- lower + 1:nk/(nk+1)*span
                  } else {
                        knot <- NULL
                  }
                  H_est2 <- as.numeric(apply(dat_test, 1, function(x) {
                        out <- max(pwl(as.numeric(x[marker_name[j]]), 
                                       fit2, knot = knot), as.numeric(x[marker_name[i]]))
                        return(out)
                  }))
                  eauc <- taua(H_est2, dat_test$event)
                  
                  H_est2_train <- as.numeric(apply(dat_train, 1, function(x) {
                        out <- max(pwl(as.numeric(x[marker_name[j]]), 
                                       fit2, knot = knot), as.numeric(x[marker_name[i]]))
                        return(out)
                  }))
                  eauc_train <- taua(H_est2_train, dat_train$event)
                  
                  out_ortree[[count]] <- list(var1 = marker_name[i],
                                              var2 = marker_name[j],
                                              fit_cv = fit_cv2,
                                              nk = nk, fit = fit2,
                                              eauc_test = eauc, 
                                              eauc_train = eauc_train)
            }
      }
}

# select trees that provide improvement over individual marker
good_tree <- list()
count <- 0
options(warn=2)
for(i in 1:length(out_andtree)) {
      # and tree
      nk <- out_andtree[[i]]$nk
      p <- nk + 2
      fit1 <- out_andtree[[i]]$fit
      if (length(fit1) != p) print(i)
      if (nk != 0) {
            lower <- range(dat_train[, out_andtree[[i]]$var1])[1]
            span <- range(dat_train[, out_andtree[[i]]$var2])[2] - 
                  range(dat_train[, out_andtree[[i]]$var2])[1]
            knot <- lower + 1:nk/(nk+1)*span
      } else {
            knot <- NULL
            
      }
      H_est1 <- as.numeric(apply(dat_test, 1, function(x) {
            out <- min(pwl(as.numeric(x[out_andtree[[i]]$var2]), 
                           fit1, knot = knot), as.numeric(x[out_andtree[[i]]$var1]))
            return(out)
      }))
      eauc <- taua(H_est1, dat_test$event)
      temp <- roc(dat_test$event, H_est1)
      
      H_est1_train <- as.numeric(apply(dat_train, 1, function(x) {
            out <- min(pwl(as.numeric(x[out_andtree[[i]]$var2]), 
                           fit1, knot = knot), as.numeric(x[out_andtree[[i]]$var1]))
            return(out)
      }))
      eauc_train <- taua(H_est1_train, dat_train$event)
      
      auc_marginal <- c(taua(dat_test[, out_andtree[[i]]$var1], dat_test$event), 
                        taua(dat_test[, out_andtree[[i]]$var2], dat_test$event))
      if (TRUE) {
            count <- count + 1
            good_tree[[count]] <- list(type = "and",
                                       fit = out_andtree[[i]],
                                       auc_marginal = auc_marginal, 
                                       eauc = paste0(ci.auc(dat_test$event, H_est1), collapse = ","),
                                       TP = temp$sensitivities[temp$sensitivities + temp$specificities == max(temp$sensitivities + temp$specificities)], 
                                       FP = 1- temp$specificities[temp$sensitivities + temp$specificities == max(temp$sensitivities + temp$specificities)])
      }
      
      # or tree
      nk <- out_ortree[[i]]$nk
      p <- nk + 2
      fit2 <- out_ortree[[i]]$fit
      if (length(fit2) != p) print(i)
      if (nk != 0) {
            lower <- range(dat_train[, out_ortree[[i]]$var1], na.rm = T)[1]
            span <- range(dat_train[, out_ortree[[i]]$var2], na.rm = T)[2] - 
                  range(dat_train[, out_ortree[[i]]$var2], na.rm = T)[1]
            knot <- lower + 1:nk/(nk+1)*span
      } else {
            knot <- NULL
      }
      H_est2 <- as.numeric(apply(dat_test, 1, function(x) {
            out <- max(pwl(as.numeric(x[out_ortree[[i]]$var2]), 
                           fit2, knot = knot), as.numeric(x[out_ortree[[i]]$var1]))
            return(out)
      }))
      eauc <- taua(H_est2, dat_test$event)
      temp <- roc(dat_test$event, H_est2)
      
      H_est2_train <- as.numeric(apply(dat_train, 1, function(x) {
            out <- max(pwl(as.numeric(x[out_ortree[[i]]$var2]), 
                           fit2, knot = knot), as.numeric(x[out_ortree[[i]]$var1]))
            return(out)
      }))
      eauc_train <- taua(H_est2_train, dat_train$event)
      
      if (TRUE) {
            count <- count + 1
            good_tree[[count]] <- list(type = "or", 
                                       fit = out_ortree[[i]], 
                                       auc_marginal = auc_marginal, 
                                       eauc = paste0(ci.auc(dat_test$event, H_est2), collapse = ","),
                                       TP = temp$sensitivities[temp$sensitivities + temp$specificities == max(temp$sensitivities + temp$specificities)], 
                                       FP = 1 - temp$specificities[temp$sensitivities + temp$specificities == max(temp$sensitivities + temp$specificities)])
      }
}



# fit CART
set.seed(37)
fit_cart <- rpart(as.formula(paste0("event~", paste0(marker_name, collapse = "+"))),
                  method = "class", data = dato_train, 
                  control = rpart.control(minsplit = 20, xval = 10))
cptab_cart <- printcp(fit_cart)
plotcp(fit_cart)
fit_cart_pruned <- prune(fit_cart, cp = cptab_cart[cptab_cart[, 4] == min(cptab_cart[-c(1, 2), 4]), 1][1])
rpart.plot(fit_cart_pruned, box.palette = 0, type = 0, xflip = F)
auc(roc(factor(dato_train$event), predict(fit_cart)[, 1]))
auc(roc(dato_test$event, predict(fit_cart, dato_test)[, 1]))
ci.auc(roc(dato_test$event, predict(fit_cart, dato_test)[, 1]))

roc_cart <- roc(factor(dato_test$event), predict(fit_cart, dato_test)[, 1])
roc_cart$sensitivities[which(roc_cart$sensitivities + roc_cart$specificities == 
                                       max(roc_cart$sensitivities + roc_cart$specificities))]
1-roc_cart$specificities[which(roc_cart$sensitivities + roc_cart$specificities == 
                                         max(roc_cart$sensitivities + roc_cart$specificities))]

set.seed(37)
rc_est_cart <- function(var1, var2, var3, var4, outcome, nk) {
      rc <- function(theta1, theta2, theta3, nknot = NULL, knot = NULL) {
            if (nk != 0) {
                  lower1 <- range(var2)[1]
                  span1 <- range(var2)[2] - range(var2)[1]
                  knot1 <- lower1 + 1:nk / (nk + 1) * span1
                  
                  knot1r <- -rev(knot1)
                  
                  lower2 <- range(var3)[1]
                  span2 <- range(var3)[2] - range(var3)[1]
                  knot2 <- lower2 + 1:nk / (nk + 1) * span2
                  
                  knot2r <- -rev(knot2)
                  
                  lower3 <- range(var4)[1]
                  span3 <- range(var4)[2] - range(var4)[1]
                  knot3 <- lower3 + 1:nk / (nk + 1) * span3
                  
                  knot3r <- -rev(knot3)
            } else {
                  knot1 <- NULL
                  knot2 <- NULL
                  knot3 <- NULL
                  knot1r <- NULL
                  knot2r <- NULL
                  knot3r <- NULL
            }
            
            H <- as.numeric(sapply(1:length(var1), function(i) {
                  out <- max(min(pwl(-as.numeric(var2[i]), theta1, knot = knot1r), 
                             as.numeric(var1[i])),
                             min(as.numeric(var1[i]),
                                 pwl(as.numeric(var3[i]), theta2, knot = knot2),
                                 pwl(-as.numeric(var4[i]), theta3, knot = knot3r)))
                  return(out)
            }))
            result <- taua(H, outcome)
      }
      
      p <- nk + 2
      optimfn <- function(theta) {
            return(rc(theta1 = theta[1:p], theta2 = theta[1:p+p],
                      theta3 = theta[1:p+2*p], nknot = nk))
      }
      
      start <- matrix(mvrnorm(n = 500, mu = rep(0, 3*p), 
                              Sigma = diag(2, 3*p)), ncol = 3*p)
      comp.sol <- do.call(rbind, lapply(1:nrow(start), function(i) {
            fit <- optim(start[i, ], optimfn, control = list(fnscale = -1))
            temp <- c(fit$par, fit$value)
      }))
      sol <- matrix(comp.sol[comp.sol[, 3*p+1] == max(comp.sol[, 3*p+1]), 1:(3*p)], ncol = 3*p)[1, ]
      return(sol)
}

rcfit_wcv_cart <- function(var1_name, var2_name, var3_name, var4_name, outcome_name, dat, 
                                  nkvec = 2, nfold = 5, seed = 37) {
      set.seed(seed)
      folds0 <- createFolds(which(dat[, outcome_name] == 0), nfold)
      folds1 <- createFolds(which(dat[, outcome_name] == 1), nfold)
      folds <- list()
      for (l in 1:nfold) {
            folds[[l]] <- c(which(dat[, outcome_name] == 0)[folds0[[l]]], 
                            which(dat[, outcome_name] == 1)[folds1[[l]]])
      }
      cvrc <- NULL
      for (nk in 0:nkvec) {
            print(nk)
            p <- nk + 2
            if (nk != 0) {
                  lower1 <- range(dat[, var2_name])[1]
                  span1 <- range(dat[, var2_name])[2] - range(dat[, var2_name])[1]
                  knot1 <- lower1 + 1:nk/(nk+1)*span1
                  knot1r <- -rev(knot1)
                  
                  lower2 <- range(dat[, var3_name])[1]
                  span2 <- range(dat[, var3_name])[2] - range(dat[, var3_name])[1]
                  knot2 <- lower2 + 1:nk/(nk+1)*span2
                  
                  knot2r <- -rev(knot2)
                  
                  lower3 <- range(dat[, var4_name])[1]
                  span3 <- range(dat[, var4_name])[2] - range(dat[, var4_name])[1]
                  knot3 <- lower3 + 1:nk/(nk+1)*span3
                  
                  knot3r <- -rev(knot3)
            } else {
                  knot1 <- NULL
                  knot2 <- NULL
                  knot3 <- NULL
                  knot1r <- NULL
                  knot2r <- NULL
                  knot3r <- NULL
            }
            tempauc <- NULL
            for (i in 1:nfold) {
                  # print(i)
                  est <- rc_est_cart(var1 = dat[-folds[[i]], var1_name], 
                                            var2 = dat[-folds[[i]], var2_name], 
                                            var3 = dat[-folds[[i]], var3_name],
                                     var4 = dat[-folds[[i]], var4_name],
                                            outcome = dat[-folds[[i]], outcome_name], nk)
                  H <- as.numeric(apply(dat[folds[[i]], ], 1, function(x) {
                        out <- max(min(pwl(-as.numeric(x[var2_name]), est[1:p], knot = knot1r), 
                                   as.numeric(x[var1_name])), 
                                   min(as.numeric(x[var1_name]),
                                       pwl(as.numeric(x[var3_name]), est[1:p+p], knot = knot2),
                                       pwl(-as.numeric(x[var4_name]), est[1:p+2*p], knot = knot3r)))
                        return(out)
                  }))
                  tempauc <- c(tempauc, taua(H, dat[, outcome_name][folds[[i]]]))
                  # evarauc <- aucvar(H, dat$prog[folds[[i]]])
            }
            print(tempauc)
            print(mean(tempauc))
            cvrc <- rbind(cvrc, c(nk, mean(tempauc)))
      }
      return(cvrc)
}

fit_cv1 <- rcfit_wcv_cart(var1_name = "ALP",
                          var2_name = "HB", 
                          var3_name = "PSA", 
                          var4_name = "BMI",
                          outcome_name = "event",
                          dat = dat_train, nkvec = 2)

nk <- fit_cv1[fit_cv1[, 2] == max(fit_cv1[, 2]), 1][1] # nk = 1
p <- nk+2


fit1 <- rc_est_cart(var1 = dat_train[, "ALP"], 
                    var2 = dat_train[, "HB"],
                    var3 = dat_train[, "PSA"],
                    var4 = dat_train[, "BMI"],
                    outcome = dat_train$event, nk)

if (nk != 0) {
      lower1 <- range(dat_train$HB)[1]
      span1 <- range(dat_train$HB)[2] - range(dat_train$HB)[1]
      knot1 <- lower1 + 1:nk/(nk+1)*span1
      knot1r <- -rev(knot1)
      
      lower2 <- range(dat_train$PSA)[1]
      span2 <- range(dat_train$PSA)[2] - range(dat_train$PSA)[1]
      knot2 <- lower2 + 1:nk/(nk+1)*span2
      knot2r <- -rev(knot2)
      
      lower3 <- range(dat_train$BMI)[1]
      span3 <- range(dat_train$BMI)[2] - range(dat_train$BMI)[1]
      knot3 <- lower3 + 1:nk/(nk+1)*span3
      knot3r <- -rev(knot3)
} else {
      knot1 <- NULL
      knot2 <- NULL
      knot3 <- NULL
      knot1r <- NULL
      knot2r <- NULL
      knot3r <- NULL
}
H_est1 <- as.numeric(apply(dat_test, 1, function(x) {
      out <- max(min(pwl(-as.numeric(x["HB"]), fit1[1:p], knot = knot1r), 
                 as.numeric(x["ALP"])), 
                 min(as.numeric(x["ALP"]),
                     pwl(as.numeric(x["PSA"]), fit1[1:p+p], knot = knot2),
                     pwl(-as.numeric(x["BMI"]), fit1[1:p+2*p], knot = knot3r)))
      return(out)
}))
eauc <- taua(H_est1, dat_test$event)
auc(roc(dat_test$event, H_est1))
ci.auc(roc(dat_test$event, H_est1))
temp <- roc(dat_test$event, H_est1)
TP <- temp$sensitivities[temp$sensitivities + temp$specificities == 
                              max(temp$sensitivities + temp$specificities)]
FP <- 1- temp$specificities[temp$sensitivities + temp$specificities == 
                                 max(temp$sensitivities + temp$specificities)]

H_est1_train <- as.numeric(apply(dat_train, 1, function(x) {
      out <- max(min(pwl(-as.numeric(x["HB"]), fit1[1:p], knot = knot1r), 
                 as.numeric(x["ALP"])), 
                 min(as.numeric(x["ALP"]), 
                     pwl(as.numeric(x["PSA"]), fit1[1:p+p], knot = knot2),
                     pwl(-as.numeric(x["BMI"]), fit1[1:p+2*p], knot = knot3r)))
      return(out)
})) # mrc estimation based on CART
eauc_train <- taua(H_est1_train, dat_train$event)

# save.image(file = "DREAM_tree_analysis_new.rda")


mrccart_roc <- roc(dat_test$event, H_est1)


c <- mrccart_roc$thresholds[mrccart_roc$sensitivities + mrccart_roc$specificities == 
                                        max(mrccart_roc$sensitivities + mrccart_roc$specificities)]
mrccart_roc$sensitivities[mrccart_roc$sensitivities + mrccart_roc$specificities == 
                                      max(mrccart_roc$sensitivities + mrccart_roc$specificities)]
1 - mrccart_roc$specificities[mrccart_roc$sensitivities + mrccart_roc$specificities == 
                                          max(mrccart_roc$sensitivities + mrccart_roc$specificities)]

postscript(file = "CART_roc.eps", width = 5, height = 5)
plot(1 - roc_cart$specificities, roc_cart$sensitivities, 
     type = "l", lty = "dashed", 
     xlab = "False positive rate", ylab = "True positive rate")
lines(1 - mrccart_roc$specificities, mrccart_roc$sensitivities,
      type = "l")
legend(0.5, 0.3 , legend = c("optimal CART", "CART"),
       lty = c("solid", "dashed"), cex = 0.8)
dev.off()

c * marker_sd[1] + marker_mean[1]

library(GoFKernel)
# HB
f <- function(x) {
   out <- pwl(-x, fit1[1:p], knot = knot1r)
   return(out)
}
inverse(f)(c) * marker_sd[2] + marker_mean[2]
# PSA
f <- function(x) {
   out <-  pwl(x, fit1[1:p+p], knot = knot2)
   return(out)
}
inverse(f)(c) * marker_sd[3] + marker_mean[3]
# BMI
f <- function(x) {
   out <- pwl(-x, fit1[1:p+2*p], knot = knot3r)
   return(out)
}
inverse(f)(c) * marker_sd[6] + marker_mean[6]

# CIT
set.seed(37)
fit_cit <- ctree(as.formula(paste0("event~", paste0(marker_name, collapse = "+"))),
                 data = dato_train)
auc(roc(factor(dato_train$event), predict(fit_cit, newdata = dato_train)))
auc(roc(factor(dato_test$event), predict(fit_cit, newdata = dato_test)))
ci.auc(roc(factor(dato_test$event), predict(fit_cit, newdata = dato_test)))

# make ctree plot
dato_train$event <- dato_train$event - 1
fit_cit <- ctree(as.formula(paste0("event~", paste0(marker_name, collapse = "+"))),
                 data = dato_train)
plot(as.simpleparty(fit_cit))
dato_train$event <- dato_train$event + 1


# rank correlation estimation using ctree structure
rc_est_cit <- function(var1, var2, var3, outcome, nk) {
   rc <- function(theta1, theta2, nknot = NULL, knot = NULL) {
      if (nk != 0) {
         lower1 <- range(var2)[1]
         span1 <- range(var2)[2] - range(var2)[1]
         knot1 <- lower1 + 1:nk / (nk + 1) * span1
         
         knot1r <- -rev(knot1)
         
         lower2 <- range(var3)[1]
         span2 <- range(var3)[2] - range(var3)[1]
         knot2 <- lower2 + 1:nk / (nk + 1) * span2
         
         knot2r <- -rev(knot2)
         
         
      } else {
         knot1 <- NULL
         knot2 <- NULL
         
         knot1r <- NULL
         knot2r <- NULL
         
      }
      
      H <- as.numeric(sapply(1:length(var1), function(i) {
         out <- min(as.numeric(var1[i]), 
                    pwl(-as.numeric(var2[i]), theta1, knot = knot1r),
                    pwl(-as.numeric(var3[i]), theta2, knot = knot2r))
         return(out)
      }))
      result <- taua(H, outcome)
   }
   
   p <- nk + 2
   optimfn <- function(theta) {
      return(rc(theta1 = theta[1:p], theta2 = theta[1:p+p],
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

rcfit_wcv_cit <- function(var1_name, var2_name, var3_name, outcome_name, dat, 
                           nkvec = 2, nfold = 5, seed = 37) {
   set.seed(seed)
   folds0 <- createFolds(which(dat[, outcome_name] == 0), nfold)
   folds1 <- createFolds(which(dat[, outcome_name] == 1), nfold)
   folds <- list()
   for (l in 1:nfold) {
      folds[[l]] <- c(which(dat[, outcome_name] == 0)[folds0[[l]]], 
                      which(dat[, outcome_name] == 1)[folds1[[l]]])
   }
   cvrc <- NULL
   for (nk in 0:nkvec) {
      print(nk)
      p <- nk + 2
      if (nk != 0) {
         lower1 <- range(dat[, var2_name])[1]
         span1 <- range(dat[, var2_name])[2] - range(dat[, var2_name])[1]
         knot1 <- lower1 + 1:nk/(nk+1)*span1
         knot1r <- -rev(knot1)
         
         lower2 <- range(dat[, var3_name])[1]
         span2 <- range(dat[, var3_name])[2] - range(dat[, var3_name])[1]
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
         est <- rc_est_cit(var1 = dat[-folds[[i]], var1_name], 
                            var2 = dat[-folds[[i]], var2_name], 
                            var3 = dat[-folds[[i]], var3_name],
                            outcome = dat[-folds[[i]], outcome_name], nk)
         H <- as.numeric(apply(dat[folds[[i]], ], 1, function(x) {
            out <- min(as.numeric(x[var1_name]),
                       pwl(-as.numeric(x[var2_name]), est[1:p], knot = knot1r),
                       pwl(-as.numeric(x[var3_name]), est[1:p+p], knot = knot2r))
            return(out)
         }))
         tempauc <- c(tempauc, taua(H, dat[, outcome_name][folds[[i]]]))
         # evarauc <- aucvar(H, dat$prog[folds[[i]]])
      }
      print(tempauc)
      print(mean(tempauc))
      cvrc <- rbind(cvrc, c(nk, mean(tempauc)))
   }
   return(cvrc)
}

dat_train$eventrev <- 1 - dat_train$event
dat_test$eventrev <- 1 - dat_test$event
set.seed(37)
fit_cit_cv1 <- rcfit_wcv_cit(var1_name = "HB",
                          var2_name = "AST", 
                          var3_name = "PSA", 
                          outcome_name = "eventrev",
                          dat = dat_train, nkvec = 2)
nk <- fit_cit_cv1[fit_cit_cv1[, 2] == max(fit_cit_cv1[, 2]), 1][1] # nk = 1
p <- nk+2 #


fit_cit_mrc <- rc_est_cit(var1 = dat_train[, "HB"], 
                    var2 = dat_train[, "AST"],
                    var3 = dat_train[, "PSA"],
                    outcome = dat_train$eventrev, nk)

if (nk != 0) {
   lower1 <- range(dat_train$AST)[1]
   span1 <- range(dat_train$AST)[2] - range(dat_train$AST)[1]
   knot1 <- lower1 + 1:nk/(nk+1)*span1
   knot1r <- -rev(knot1)
   
   lower2 <- range(dat_train$PSA)[1]
   span2 <- range(dat_train$PSA)[2] - range(dat_train$PSA)[1]
   knot2 <- lower2 + 1:nk/(nk+1)*span2
   knot2r <- -rev(knot2)
} else {
   knot1 <- NULL
   knot2 <- NULL
   knot1r <- NULL
   knot2r <- NULL
}
H_cit_mrc_est <- as.numeric(apply(dat_test, 1, function(x) {
   out <- min(as.numeric(x["HB"]),
              pwl(-as.numeric(x["AST"]), fit_cit_mrc[1:p], knot = knot1r),
              pwl(-as.numeric(x["PSA"]), fit_cit_mrc[1:p+p], knot = knot2r))
   return(out)
}))


eauc <- taua(H_cit_mrc_est, dat_test$eventrev)
auc(roc(dat_test$eventrev, H_cit_mrc_est))
ci.auc(roc(dat_test$eventrev, H_cit_mrc_est))
temp <- roc(dat_test$eventrev, H_cit_mrc_est)
TP <- temp$sensitivities[temp$sensitivities + temp$specificities == 
                            max(temp$sensitivities + temp$specificities)]
FP <- 1- temp$specificities[temp$sensitivities + temp$specificities == 
                               max(temp$sensitivities + temp$specificities)]

H_cit_mrc_train <- as.numeric(apply(dat_train, 1, function(x) {
   
   
   out <- max(min(pwl(-as.numeric(x["HB"]), fit1[1:p], knot = knot1r), 
                  as.numeric(x["ALP"])), 
              min(as.numeric(x["ALP"]), 
                  pwl(as.numeric(x["PSA"]), fit1[1:p+p], knot = knot2),
                  pwl(-as.numeric(x["BMI"]), fit1[1:p+2*p], knot = knot3r)))
   return(out)
}))
eauc_cit_mrc_train <- taua(H_cit_mrc_train, dat_train$eventrev)


roc_cit <- roc(factor(dato_test$event), predict(fit_cit, dato_test))
roc_cit$sensitivities[which(roc_cit$sensitivities + roc_cit$specificities == 
                                    max(roc_cit$sensitivities + roc_cit$specificities))]
1-roc_cit$specificities[which(roc_cit$sensitivities + roc_cit$specificities == 
                                      max(roc_cit$sensitivities + roc_cit$specificities))]


mrccit_roc <- roc(factor(1-dato_test$event), H_cit_mrc_est)
mrccit_roc$sensitivities[which(mrccit_roc$sensitivities + mrccit_roc$specificities == 
                                   max(mrccit_roc$sensitivities + mrccit_roc$specificities))]
1-mrccit_roc$specificities[which(mrccit_roc$sensitivities + mrccit_roc$specificities == 
                                     max(mrccit_roc$sensitivities + mrccit_roc$specificities))]
# find cutoffs
postscript(file = "CIT_roc.eps", width = 5, height = 5)
plot(1 - roc_cit$specificities, roc_cit$sensitivities, 
     type = "l", lty = "dashed",
     xlab = "False positive rate", ylab = "True positive rate")
lines(1 - mrccit_roc$specificities, mrccit_roc$sensitivities,
      type = "l")
legend(0.5, 0.3, legend = c("optimal CIT", "CIT"), 
       lty = c("solid", "dashed"), cex = 0.8)
dev.off()

c <- mrccit_roc$thresholds[mrccit_roc$sensitivities + mrccit_roc$specificities == 
                               max(mrccit_roc$sensitivities + mrccit_roc$specificities)]
c * marker_sd[2] + marker_mean[2] # HB
# AST
f <- function(x) {
   out <- pwl(-x, fit_cit_mrc[1:p], knot = knot1r)
   return(out)
}
inverse(f)(c) * marker_sd[5] + marker_mean[5]
# PSA
f <- function(x) {
   out <- pwl(-x, fit_cit_mrc[1:p+p], knot = knot2r)
   return(out)
}
inverse(f)(c) * marker_sd[3] + marker_mean[3]
