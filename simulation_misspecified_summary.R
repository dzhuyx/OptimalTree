rm(list = ls())
library(rpart)
library(partykit)
library(pROC)
load("trueROC.rda")

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

n <- 100
fc <- "misfcpoly2"

summaryresult_misfc <- data.frame()
for (tree in paste0("tree", 1:4)) {
      for (rho in c(0.2, 0.5, 0.8)) {
            for (delta in c(1, 3)) {
                  ROC_true <- dat_trueROC[dat_trueROC$fc == "poly2" & 
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
                  
                  print(paste(fc, tree, rho, delta, n))
                  # read simulation result .rda files
                  main.dir <- getwd()
                  sub.dir <- paste0("mrcsim-", fc, tree,"-r", rho*10, "-d", delta, "-n", n)
                  filenames <- list.files(path = sub.dir)
                  # load(paste0(sub.dir, "/", filenames[1]))
                  # p <- length(theta0)
                  
                  result <- do.call(rbind, lapply(filenames, function(x) {
                        load(paste0(sub.dir, "/", x))
                        
                        # calculate ROC
                        roc_train_nk0 <- roc(factor(Y), H_train_nk0)
                        roc_train_nk1 <- roc(factor(Y), H_train_nk1)
                        roc_train_nk2 <- roc(factor(Y), H_train_nk2)
                        
                        roc_test_nk0 <- roc(factor(Y), H_test_nk0)
                        roc_test_nk1 <- roc(factor(Y), H_test_nk1)
                        roc_test_nk2 <- roc(factor(Y), H_test_nk2)
                        
                        # calculate AUC
                        auc_train_nk0 <- auc(roc_train_nk0)
                        auc_train_nk1 <- auc(roc_train_nk1)
                        auc_train_nk2 <- auc(roc_train_nk2)
                        
                        auc_test_nk0 <- auc(roc_test_nk0)
                        auc_test_nk1 <- auc(roc_test_nk1)
                        auc_test_nk2 <- auc(roc_test_nk2)
                        
                        # calculate max abs difference in ROC
                        rocdiff_train_nk0 <- myrocdiff(roc_train_nk0, 
                                                       sens_sub =  sens_sub,
                                                       spec_sub = spec_sub)
                        rocdiff_train_nk1 <- myrocdiff(roc_train_nk1, 
                                                       sens_sub =  sens_sub,
                                                       spec_sub = spec_sub)
                        rocdiff_train_nk2 <- myrocdiff(roc_train_nk2, 
                                                       sens_sub =  sens_sub,
                                                       spec_sub = spec_sub)
                        
                        rocdiff_test_nk0 <- myrocdiff(roc_test_nk0, 
                                                       sens_sub =  sens_sub,
                                                       spec_sub = spec_sub)
                        rocdiff_test_nk1 <- myrocdiff(roc_test_nk1, 
                                                       sens_sub =  sens_sub,
                                                       spec_sub = spec_sub)
                        rocdiff_test_nk2 <- myrocdiff(roc_test_nk2, 
                                                       sens_sub =  sens_sub,
                                                       spec_sub = spec_sub)
                  
                        return(c(auc_train_nk0, 
                                 auc_train_nk1, 
                                 auc_train_nk2,
                                 auc_test_nk0, 
                                 auc_test_nk1, 
                                 auc_test_nk2,
                                 rocdiff_train_nk0, 
                                 rocdiff_train_nk1, 
                                 rocdiff_train_nk2,
                                 rocdiff_test_nk0, 
                                 rocdiff_test_nk1, 
                                 rocdiff_test_nk2))
                  }))
                  
                  # summarize results
                  auc_bias <- round(colMeans(result[, 1:6]) - auc_true, 3)
                  auc_ese <- round(apply(result[, 1:6], 2, sd), 3)
                  auc_mse <- round(apply(result[, 1:6], 2, function(x) {
                        return(sqrt(mean((x - auc_true)^2)))
                  }), 3)
                  auc_summary <- sapply(1:6, function(j) {
                        return(paste0(auc_bias[j], " (",
                                      auc_ese[j], ", ",
                                      auc_mse[j], ")"))
                  })
                  
                  rocdiff_bias <- round(colMeans(result[, 7:12]), 3)
                  rocdiff_ese <- round(apply(result[, 7:12], 2, sd), 3)
                  rocdiff_summary <- sapply(1:6, function(j) {
                        return(paste0(rocdiff_bias[j], " (",
                                      rocdiff_ese[j], ")"))
                  })
                  
                  summaryresult_misfc <- rbind(summaryresult_misfc, 
                                               data.frame(cbind(t(auc_summary),
                                                                t(rocdiff_summary))))
            }
      }
}


# reformat summarized data.frame into manuscript tables
for (tree in paste0("tree", 1:4)) {
      for (rho in c(0.2, 0.5, 0.8)) {
            for (delta in c(1, 3)) { 
                  tempdat <- rbind(tempdat, data.frame(tree = tree, 
                                                       rho = rho, 
                                                       delta = delta))
            }
      }
}
summaryresult_misfc <- cbind(tempdat, summaryresult_misfc)
tempbias_auc <- apply(summaryresult_misfc[, 4:13], 1:2, function(x) {
      return(as.numeric(as.character(strsplit(x, " (", fixed = T)[[1]][1])))
})
tempese_auc <-  apply(summaryresult_misfc[, 4:13], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][1]))
})
tempmse_auc <-  apply(summaryresult_misfc[, 4:13], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ",", fixed = T)[[1]][2], ")", fixed = T)[[1]][1]))
})
tempbias_rocdiff <- apply(summaryresult_misfc[, 14:23], 1:2, function(x) {
      return(as.numeric(as.character(strsplit(x, " (", fixed = T)[[1]][1])))
})
tempese_rocdiff <-  apply(summaryresult_misfc[, 14:23], 1:2, function(x) {
      return(as.numeric(strsplit(strsplit(x, " (", fixed = T)[[1]][2], ")", fixed = T)[[1]][1]))
})

tempbias_auc <- formatC(tempbias_auc, format = "f", digits = 3)
tempese_auc <- formatC(tempese_auc, format = "f", digits = 3)
tempmse_auc <- formatC(tempmse_auc, format = "f", digits = 3)
tempbias_rocdiff <- formatC(tempbias_rocdiff, format = "f", digits = 3)
tempese_rocdiff <- formatC(tempese_rocdiff, format = "f", digits = 3)

tempresults <- cbind(summaryresult_misfc[, 1:3], sapply(c(4:23), function(j) {
      if (j %in% 4:13) {
            out <- paste0(tempbias_auc[, j-3], " (",
                          tempese_auc[, j-3], ", ",
                          tempmse_auc[, j-3], ")")
      } else {
            out <- paste0(tempbias_rocdiff[, j-13], " (",
                          tempese_rocdiff[, j-13], ")")
      }
      return(out)
}))



newresult <- data.frame()
for (tree in paste0("tree", 1:4)) {
      for (rho in c(0.2, 0.5, 0.8)) {
            for (delta in c(1, 3)) {
                  for (K in 0:2) {
                        tempdat <- tempresults[tempresults$tree == tree &
                                                             tempresults$rho == rho &
                                                             tempresults$delta == delta, ]
                        newresult <- rbind(newresult, 
                                           c(unlist(tempdat[1, 1:3]),
                                                 unlist(c(K, tempdat[1, 4+K], tempdat[1, 9+K], 
                                                              tempdat[1, 14+K], tempdat[1, 19+K]))))
                  }
            }
      }
}
write.csv(newresult[newresult[, 1] == "tree1", ], 
          file = paste0("summaryresult_misfc_tree1.csv"))
write.csv(newresult[newresult[, 1] == "tree2", ], 
          file = paste0("summaryresult_misfc_tree2.csv"))
write.csv(newresult[newresult[, 1] == "tree3", ], 
          file = paste0("summaryresult_misfc_tree3.csv"))
write.csv(newresult[newresult[, 1] == "tree4", ], 
          file = paste0("summaryresult_misfc_tree4.csv"))
