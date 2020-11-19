# setwd("/users/yzhu/treeROC")

# for this script to run
# variables n, delta, rho need to be
# assigned values
# see quoted block codes for an example
# n <- 50
# delta <- 1
# rho <- 0.2


filenames <- list.files(path = paste0("mrc-linear7-d", delta, "-r", rho, "-n", n))
cwd <- getwd()
setwd(paste0("mrc-linear7-d", delta, "-r", rho, "-n", n))

simres <- t(sapply(filenames, function(x) {
	load(x)
	tpe <- res[["pe"]]
	tpe[c(1, 3)] <- exp(tpe[c(1, 3)])
	tboot <- do.call(rbind, lapply(res[["boot"]], function(x) {
		return(x)
	}))
	# tboot[, c(1, 3)] <- exp(tboot[, c(1, 3)])
	tse <- sqrt(diag(var(tboot)))
	tse[c(1, 3)] <- tse[c(1, 3)] * tpe[c(1, 3)]
	tcp <- (tpe - 1.96 * tse <= c(1, -1, 2, 0.5)) & (tpe + 1.96 * tse >= c(1, -1, 2, 0.5))
	# tcp <- (c(log(1), -1, log(2), 0.5) <= quantile(tboot, 0.975)) &
	#    	   (c(log(1), -1, log(2), 0.5) >= quantile(tboot, 0.025))
	return(c(tpe, tse, tcp))
}))

# simres <- simres[simres[, 3] < 25, ]
bias <- colMeans(simres[, 1:4]) - c(1, -1, 2, 0.5)
ese <- sqrt(diag(var(simres[, 1:4])))
mse <- colMeans(simres[, 5:8])
cp <- colMeans(simres[, 9:12])


length(filenames)

bias 
ese
mse
cp

# summary(simres)

setwd(cwd)
