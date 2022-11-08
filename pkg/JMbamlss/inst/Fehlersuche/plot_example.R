source("R/JM_example.R")


par(mfrow = c(1,2))
curve(1.4*log((120*x + 10)/1000), from = 0, to = 1)

pred_data <- data.frame(seq(0, 1, length.out = 100))
colnames(pred_data) <- b$x$lambda$smooth.construct[[1]]$term

k <- b$x$lambda$smooth.construct[[1]]$bs.dim
b_it <- b$x$lambda$smooth.construct[[1]]$state$parameters[-k]

pred_mat <- PredictMat(b$x$lambda$smooth.construct[[1]], pred_data)
plot(pred_data[, 1], pred_mat%*%b$parameters$lambda$s$`s(survtime)`[1:9])
