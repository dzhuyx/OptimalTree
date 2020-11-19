rm(list = ls())
library(rpart)
library(partykit)
library(pROC)
load("trueROC.rda") # true ROC calculated using `simulation_correct_trueROC.R`

# summarize rda files from simulations with 
# correctly specified models into data.frame
# function for calculating maximum absolute 
# difference in two ROC curves
myrocdiff <- function(roc_obj, sens_sub, spec_sub) {
      # sens is increasing
      # spec is decreasing
      out <- sapply(1:(length(sens_sub)-1), function(i) {
            senstemp <- sort(roc_obj$sensitivities)
            spectemp <- sort(roc_obj$specificities, decreasing = T)
            pin <- max(which(senstemp <= sens_sub[i]))
            return(spectemp[pin] + (spectemp[pin + 1] - spectemp[pin]) / 
                         (senstemp[pin + 1] - senstemp[pin]) * 
                         (sens_sub[i] - senstemp[pin]))
      })
      return(max(abs(out - spec_sub[-length(spec_sub)])))
}
# n <- 100
# delta <- 1
# rho <- 0.2
# tree <- "tree1"
# fc <- "linear2"

# summarize point estimate and ROC results
# pe is of different length for different fc

summaryresult_pe0 <- data.frame()
summaryresult_pe1 <- data.frame()
summaryresult_pe2 <- data.frame()
summaryresult_pe3 <- data.frame()
summaryresult_roc <- data.frame()
for (fc in c("linear0", "linear1", "linear2", "poly2")) {
      for (tree in paste0("tree", 1:4)) {
            for (rho in c(0.2, 0.5, 0.8)) {
                  for (delta in c(1, 3)) {
                        ROC_true <- dat_trueROC[dat_trueROC$fc == fc & 
                                                      dat_trueROC$tree == tree &
                                                      dat_trueROC$delta == delta &
                                                      dat_trueROC$rho == rho, ]
                        # calculate true AUC
                        sens <- sort(ROC_true$sens) # increasing
                        spec <- rev(sort(ROC_true$spec)) # decreasing
                        auc_true <- sum((sens[-1] + sens[-length(sens)]) * 
                                              diff(1-spec) / 2)
                        
                        sens_sub <- 0.001 * 0:1000
                        spec_sub <- sapply(sens_sub, function(x) {
                              return(spec[max(which(sens <= x))])
                        })
                        
                        
                        for (n in c(100, 200, 500)) {
                              print(paste(fc, tree, rho, delta, n))
                              # read simulation result .rda files
                              main.dir <- getwd()
                              sub.dir <- paste0("mrcsim-", fc, tree,"-r", rho*10, "-d", delta, "-n", n)
                              filenames <- list.files(path = sub.dir)
                              load(paste0(sub.dir, "/", filenames[1]))
                              p <- length(theta0)
                              
                              result <- do.call(rbind, lapply(filenames, function(x) {
                                    load(paste0(sub.dir, "/", x))
                                    
                                    # max difference from true ROC
                                    rocdiff_pe_train <- myrocdiff(roc_pe_train, 
                                                                  sens_sub = sens_sub,
                                                                  spec_sub = spec_sub)
                                    rocdiff_cart_train <- myrocdiff(roc_cart_train, 
                                                                    sens_sub = sens_sub,
                                                                    spec_sub = spec_sub)
                                    rocdiff_cit_train <- myrocdiff(roc_cit_train, 
                                                                   sens_sub = sens_sub,
                                                                   spec_sub = spec_sub)
                                    
                                    rocdiff_pe_test <- myrocdiff(roc_pe_test,
                                                                 sens_sub = sens_sub,
                                                                 spec_sub = spec_sub)
                                    rocdiff_cart_test <- myrocdiff(roc_cart_test,
                                                                   sens_sub = sens_sub,
                                                                   spec_sub = spec_sub)
                                    rocdiff_cit_test <- myrocdiff(roc_cit_test,
                                                                  sens_sub = sens_sub,
                                                                  spec_sub = spec_sub)
                                    
                                    out <- c(theta0, pe, auc_pe_train = auc_pe_train,
                                             auc_cart_train = auc_cart_train, 
                                             auc_cit_train = auc_cit_train,
                                             auc_pe_test = auc_pe_test,
                                             auc_cart_test = auc_cart_test,
                                             auc_cit_test = auc_cit_test,
                                             rocdiff_pe_train = rocdiff_pe_train,
                                             rocdiff_cart_train = rocdiff_cart_train,
                                             rocdiff_cit_train = rocdiff_cit_train,
                                             rocdiff_pe_test = rocdiff_pe_test,
                                             rocdiff_cart_test = rocdiff_cart_test,
                                             rocdiff_cit_test = rocdiff_cit_test)
                                    return(out)
                              }))
                              
                              
                              # summarize results
                              # point estimate results
                              pe_bias <- round(colMeans(result[, 1:p + p]), 3)
                              pe_ese <- round(apply(result[, 1:p + p], 2, sd), 3)
                              pe_mse <- round(apply(result[, 1:p+p] - result[, 1:p], 2, function(x) {
                                    return(sqrt(mean(x^2)))
                              }), 3) # is actually the sqrt of mse
                              pe_summary <- sapply(1:p, function(j) {
                                    return(paste0(pe_bias[j], " (", 
                                                  pe_ese[j], ", ", 
                                                  pe_mse[j], ")"))
                              })
                              
                              # roc results
                              auc_bias <- round(colMeans(result[, (1 + 2*p):(6+2*p)]) - auc_true, 3)
                              auc_ese <- round(apply(result[, (1 + 2*p):(6+2*p)], 2, sd), 3)
                              auc_mse <- round(apply(result[, (1 + 2*p):(6+2*p)], 2, function(x) {
                                    return(sqrt(mean((x - auc_true)^2)))
                              }), 3)
                              auc_summary <- sapply(1:6, function(j) {
                                    return(paste0(auc_bias[j], " (", 
                                                  auc_ese[j], ", ", 
                                                  auc_mse[j], ")"))
                              })
                              rocdiff_bias <- round(colMeans(result[, (7+2*p):(12+2*p)]), 3)
                              rocdiff_ese <- round(apply(result[, (7+2*p):(12+2*p)], 2, sd), 3)
                              rocdiff_summary <- sapply(1:6, function(j) {
                                    return(paste0(rocdiff_bias[j], " (", 
                                                  rocdiff_ese[j], ")"))
                              })
                              
                              # store into result data.frame
                              summaryresult_roc <- rbind(summaryresult_roc, 
                                                         cbind(data.frame(fc = fc, 
                                                                          tree = tree,
                                                                          rho = rho,
                                                                          delta = delta,
                                                                          n = n),
                                                               t(auc_summary),
                                                               t(rocdiff_summary)))
                              if (fc == "linear0") {
                                    summaryresult_pe0 <- rbind(summaryresult_pe0, 
                                                               cbind(data.frame(fc = fc, 
                                                                                tree = tree,
                                                                                rho = rho,
                                                                                delta = delta, 
                                                                                n = n),
                                                                     t(pe_summary)))
                              } else if (fc == "linear1") {
                                    summaryresult_pe1 <- rbind(summaryresult_pe1, 
                                                               cbind(data.frame(fc = fc, 
                                                                                tree = tree,
                                                                                rho = rho,
                                                                                delta = delta, 
                                                                                n = n),
                                                                     t(pe_summary)))
                              } else if (fc == "linear2") {
                                    summaryresult_pe2 <- rbind(summaryresult_pe2, 
                                                               cbind(data.frame(fc = fc, 
                                                                                tree = tree,
                                                                                rho = rho,
                                                                                delta = delta, 
                                                                                n = n),
                                                                     t(pe_summary)))
                              } else if (fc == "poly2") {
                                    summaryresult_pe3 <- rbind(summaryresult_pe3, 
                                                               cbind(data.frame(fc = fc, 
                                                                                tree = tree,
                                                                                rho = rho,
                                                                                delta = delta, 
                                                                                n = n),
                                                                     t(pe_summary)))
                              } else {
                                print("Function class not coded.")
                              }
                        }
                  }
            }
      }
}

colnames(summaryresult_roc) <- c("fc", "tree", "rho", "delta", "n",
                                 "auc_pe_train", "auc_cart_train", "auc_cit_train",
                                 "auc_pe_test", "auc_cart_test", "auc_cit_test",
                                 "rocdiff_pe_train", "rocdiff_cart_train", "rocdiff_cit_train",
                                 "rocdiff_pe_test", "rocdiff_cart_test", "rocdiff_cit_test")


# reformat summarized data.frame into manuscript tables
# linear0
theta0 <- c(log(1), -1, log(2), 0.5)
tempest <- apply(summaryresult_pe0[, 6:9], 1:2, function(x) {
      return(as.numeric(as.character(strsplit(x, " (", fixed = T)[[1]][1])))
})
tempse1 <- apply(summaryresult_pe0[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][1]))
})
tempse2 <- apply(summaryresult_pe0[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][2], ")", fixed = T)[[1]][1]))
})
# format to three decimals
tempse1 <- formatC(tempse1, format = "f", digits = 3)
tempse2 <- formatC(tempse2, format = "f", digits = 3)
newresult <- cbind(summaryresult_pe0[, 1:5], sapply(1:length(theta0), function(j) {
      out <- paste0(formatC(tempest[, j] - theta0[j], 
                            format = "f", digits = 3), 
                    " (", tempse1[, j], ", ", 
                    tempse2[, j], ")")
      return(out)
}))
write.csv(newresult[newresult$tree == "tree1", ], file = "summaryresult_pe0_tree1.csv")
write.csv(newresult[newresult$tree == "tree2", ], file = "summaryresult_pe0_tree2.csv")
write.csv(newresult[newresult$tree == "tree3", ], file = "summaryresult_pe0_tree3.csv")
write.csv(newresult[newresult$tree == "tree4", ], file = "summaryresult_pe0_tree4.csv")

# linear1
theta0 <- c(log(1), -1, log(3), log(2), 0.5, log(0.5))
tempest <- apply(summaryresult_pe1[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(as.character(strsplit(x, " (", fixed = T)[[1]][1])))
})
tempse1 <- apply(summaryresult_pe1[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][1]))
})
tempse2 <- apply(summaryresult_pe1[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][2], ")", fixed = T)[[1]][1]))
})
# format to three decimals
tempse1 <- formatC(tempse1, format = "f", digits = 3)
tempse2 <- formatC(tempse2, format = "f", digits = 3)
newresult <- cbind(summaryresult_pe1[, 1:5], sapply(1:length(theta0), function(j) {
      out <- paste0(formatC(tempest[, j] - theta0[j], 
                            format = "f", digits = 3), 
                    "(", tempse1[, j], ", ", 
                    tempse2[, j], ")")
      return(out)
}))
write.csv(newresult[newresult$tree == "tree1", ], file = "summaryresult_pe1_tree1.csv")
write.csv(newresult[newresult$tree == "tree2", ], file = "summaryresult_pe1_tree2.csv")
write.csv(newresult[newresult$tree == "tree3", ], file = "summaryresult_pe1_tree3.csv")
write.csv(newresult[newresult$tree == "tree4", ], file = "summaryresult_pe1_tree4.csv")


# linear2
theta0 <- c(log(0.5), -1.5, log(1), log(3),
            log(1), -0.5, log(2), log(0.5))
tempest <- apply(summaryresult_pe2[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(as.character(strsplit(x, " (", fixed = T)[[1]][1])))
})
tempse1 <- apply(summaryresult_pe2[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][1]))
})
tempse2 <- apply(summaryresult_pe2[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][2], ")", fixed = T)[[1]][1]))
})
# format to three decimals
tempse1 <- formatC(tempse1, format = "f", digits = 3)
tempse2 <- formatC(tempse2, format = "f", digits = 3)
newresult <- cbind(summaryresult_pe2[, 1:5], sapply(1:length(theta0), function(j) {
      out <- paste0(formatC(tempest[, j] - theta0[j], 
                            format = "f", digits = 3), 
                    "(", tempse1[, j], ", ", 
                    tempse2[, j], ")")
      return(out)
}))
write.csv(newresult[newresult$tree == "tree1", ], file = "summaryresult_pe2_tree1.csv")
write.csv(newresult[newresult$tree == "tree2", ], file = "summaryresult_pe2_tree2.csv")
write.csv(newresult[newresult$tree == "tree3", ], file = "summaryresult_pe2_tree3.csv")
write.csv(newresult[newresult$tree == "tree4", ], file = "summaryresult_pe2_tree4.csv")

# poly2
load("mrc_summary.rda")
theta0 <- c(log(1/2), log(1/20), -0.5, log(1/5), log(1/20),
            log(1/10), log(1/5), -2, log(1), log(1/3))
tempest <- apply(summaryresult_pe3[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(as.character(strsplit(x, " (", fixed = T)[[1]][1])))
})
tempse1 <- apply(summaryresult_pe3[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][1]))
})
tempse2 <- apply(summaryresult_pe3[, 6:(5+length(theta0))], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][2], ")", fixed = T)[[1]][1]))
})
# format to three decimals
tempse1 <- formatC(tempse1, format = "f", digits = 3)
tempse2 <- formatC(tempse2, format = "f", digits = 3)
newresult <- cbind(summaryresult_pe3[, 1:5], sapply(1:length(theta0), function(j) {
      out <- paste0(formatC(tempest[, j] - theta0[j], 
                            format = "f", digits = 3), 
                    "(", tempse1[, j], ", ", 
                    tempse2[, j], ")")
      return(out)
}))
write.csv(newresult[newresult$tree == "tree1", ], file = "summaryresult_pe3_tree1.csv")
write.csv(newresult[newresult$tree == "tree2", ], file = "summaryresult_pe3_tree2.csv")
write.csv(newresult[newresult$tree == "tree3", ], file = "summaryresult_pe3_tree3.csv")
write.csv(newresult[newresult$tree == "tree4", ], file = "summaryresult_pe3_tree4.csv")


# roc
colnames(summaryresult_roc) <- c("fc", "tree", "rho", "delta", "n",
                                 "auc_pe_train", "auc_cart_train", "auc_cit_train",
                                 "auc_pe_test", "auc_cart_test", "auc_cit_test",
                                 "rocdiff_pe_train", "rocdiff_cart_train", "rocdiff_cit_train",
                                 "rocdiff_pe_test", "rocdiff_cart_test", "rocdiff_cit_test")
tempest_auc <- apply(summaryresult_roc[, 6:11], 1:2, function(x) {
      return(as.numeric(as.character(strsplit(x, " (", fixed = T)[[1]][1])))
})
tempse1_auc <- apply(summaryresult_roc[, 6:11], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][1]))
})
tempse2_auc <- apply(summaryresult_roc[, 6:11], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][2], ")", fixed = T)[[1]][1]))
})

tempest_rocdiff <- apply(summaryresult_roc[, 12:17], 1:2, function(x) {
      return(as.numeric(as.character(strsplit(x, " (", fixed = T)[[1]][1])))
})
tempse_rocdiff <- apply(summaryresult_roc[, 12:17], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ")", fixed = T)[[1]][1]))
})
# format to three decimals
tempest_auc <- formatC(tempest_auc, format = "f", digits = 3)
tempse1_auc <- formatC(tempse1_auc, format = "f", digits = 3)
tempse2_auc <- formatC(tempse2_auc, format = "f", digits = 3)
tempest_rocdiff <- formatC(tempest_rocdiff, format = "f", digits = 3)
tempse_rocdiff <- formatC(tempse_rocdiff, format = "f", digits = 3)

tempresults <- cbind(summaryresult_roc[, 1:5], sapply(c(6:8, 12:14, 9:11, 15:17), function(j) {
      if (j %in% c(6:8, 9:11)) {
            out <- paste0(tempest_auc[, j-5], " (",
                          tempse1_auc[, j-5], ", ",
                          tempse2_auc[, j-5], ")")
      } else {
            out <- paste0(tempest_rocdiff[, j-11], " (",
                          tempse_rocdiff[, j-11], ")")
      }
      return(out)
}))

newresult <- data.frame()
for (fc in c(paste0("linear", 0:2), "poly2")) {
      for (tree in paste0("tree", 1:4)) {
            for (rho in c(0.2, 0.5, 0.8)) {
                  for (delta in c(1, 3)) {
                        for (n in c(100, 200, 500)) {
                              tempdat <- tempresults[tempresults$fc == fc & 
                                                           tempresults$tree == tree &
                                                           tempresults$rho == rho &
                                                           tempresults$delta == delta &
                                                           tempresults$n == n, ]
                              newresult <- rbind(newresult, 
                                                 cbind(rbind(tempdat[1, 1:5], tempdat[1, 1:5]),
                                                 rbind(unlist(c("train", tempdat[1, 6:11])),
                                                 unlist(c("test", tempdat[1, 12:17])))))
                        }
                  }
            }
      }
}
newresult <- data.frame(newre)

fc <- "poly2"
write.csv(newresult[newresult$fc == fc & newresult$tree == "tree1", ], 
          file = paste0("summaryresult_roc", "_", fc, "_tree1.csv"))
write.csv(newresult[newresult$fc == fc & newresult$tree == "tree2", ], 
          file = paste0("summaryresult_roc", "_", fc, "_tree2.csv"))
write.csv(newresult[newresult$fc == fc & newresult$tree == "tree3", ], 
          file = paste0("summaryresult_roc", "_", fc, "_tree3.csv"))
write.csv(newresult[newresult$fc == fc & newresult$tree == "tree4", ], 
          file = paste0("summaryresult_roc", "_", fc, "_tree4.csv"))