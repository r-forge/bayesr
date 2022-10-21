sim_results <- function(result_list, dat_list, name = "est") {
  
  mapply(function (est, sim) {
    
    # Bias, MSE, Coverage
    eval_lambga <- data.frame(
      type = c("Bias", "MSE", "Coverage"),
      model = name,
      predictor = "lambga",
      marker = "all",
      t = "all",
      value = c(mean(-(est$lambga$Mean - sim$lambga)),
                mean((est$lambga$Mean - sim$lambga)^2),
                mean(est$lambga[, 1] < sim$lambga & 
                       est$lambga[, 3] > sim$lambga)))
    
    eval_alpha = data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 6),
      model = name,
      predictor = "alpha",
      marker = paste0("m", 1:6),
      t = "all",
      value = c(-(est$alpha$Mean - sim$alpha),
                (est$alpha$Mean - sim$alpha)^2,
                as.numeric(est$alpha[, 1] < sim$alpha & 
                             est$alpha[, 3] > sim$alpha)))
    
    eval_sigma = data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 6),
      model = name,
      predictor = "sigma",
      marker = paste0("m", 1:6),
      t = "all",
      value = c(-(est$sigma$Mean - sim$sigma),
                (est$sigma$Mean - sim$sigma)^2,
                as.numeric(est$sigma[, 1] < sim$sigma & 
                             est$sigma[, 3] > sim$sigma)))
    
    eval_mu <- data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 6),
      model =  name,
      predictor = "mu",
      marker = paste0("m", 1:6),
      t = "all",
      value = c(mapply(function (e, s) {
                  mean(s - e)
                }, e = split(est$mu$Mean, est$mu$marker), 
                s = split(sim$mu$mu, sim$mu$marker)),
                mapply(function (e, s) {
                  mean((e - s)^2)
                }, e = split(est$mu$Mean, est$mu$marker), 
                s = split(sim$mu$mu, sim$mu$marker)),
                mapply(function (l, u, s) {
                  mean(as.numeric(l < s & u > s))
                }, l = split(est$mu[, 1], est$mu$marker), 
                u = split(est$mu[, 3], est$mu$marker),
                s = split(sim$mu$mu, sim$mu$marker))))
    
    sim_marker <- split(sim$mu_long, est$mu_long$marker)
    eval_mu_long <- data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 6*101),
      model =  name,
      predictor = "mu_long",
      marker = paste0("m", 1:6),
      t = rep(rep(seq(0, 1, by = 0.01), each = 6), times = 3),
      value = c(c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    mean(-(e$Mean[same_t] - s[same_t]))
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                })),
                c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    mean((e$Mean[same_t] - s[same_t])^2)
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                })),
                c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    mean(e[same_t, 1] < s[same_t] & e[same_t, 3] > s[same_t])
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                }))))
    
    rbind(eval_lambga, eval_alpha, eval_sigma, eval_mu, eval_mu_long)
    
  }, est = result_list, sim = dat_list, SIMPLIFY = FALSE)
  
}
