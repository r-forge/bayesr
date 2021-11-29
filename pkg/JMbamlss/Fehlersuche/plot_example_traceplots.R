source("R/JM_example.R")
library(tidyverse)

extract_fun <- function(param_ = 1, partial_ = 1, part_ = "s",
                        pred_ = "lambda", it_param_ = it_param) {
  sapply(it_param_, function(x) {
    if(part_ == "s") {
      x[[pred_]][[part_]][[partial_]][param_]
    } else {
      x[[pred_]][["p"]][param_]
    }
  })
}

plot_fun <- function(pred = "lambda", part = "s", partial = 1, 
                     param_list = it_param) {
  
  pcre <- pred == "mu" & part == "s" & partial == 3
  seq <- if (pcre) 100 else if (part == "s") 9 else length(it_param[[1]][[pred]][["p"]])
  dat <- data.frame(lapply(seq_len(seq), extract_fun, partial_ = partial,
                           part_ = part, pred_ = pred))
  names(dat) <- paste0("b", seq_len(seq))
  dat$it <- seq_len(length(param_list))
  dat <- pivot_longer(data = dat, seq_len(seq))
  
  if (pcre) {
    p1 <- ggplot(dat %>% filter(name %in% paste0("b", 1:25)),
                 aes(x = it, y = value)) +
      geom_point() +
      facet_wrap(~ name, scales = "free") +
      ggtitle(paste(pred, part, partial, "1/4"))
    p2 <- ggplot(dat %>% filter(name %in% paste0("b", 26:50)),
                 aes(x = it, y = value)) +
      geom_point() +
      facet_wrap(~ name, scales = "free") +
      ggtitle(paste(pred, part, partial, "2/4"))
    p3 <- ggplot(dat %>% filter(name %in% paste0("b", 51:75)),
                 aes(x = it, y = value)) +
      geom_point() +
      facet_wrap(~ name, scales = "free") +
      ggtitle(paste(pred, part, partial, "3/4"))
    p4 <- ggplot(dat %>% filter(name %in% paste0("b", 76:98)),
                 aes(x = it, y = value)) +
      geom_point() +
      facet_wrap(~ name, scales = "free") +
      ggtitle(paste(pred, part, partial, "3/4"))
    list(p1, p2, p3, p4)
  } else{
    ggplot(dat, aes(x = it, y = value)) +
      geom_point() +
      facet_wrap(~ name, scales = "free") +
      ggtitle(paste(pred, part, partial))
    
  }

}
plot_fun("lambda", "s", 1)
plot_fun("gamma", "p")
plot_fun("mu", "p")
plot_fun("mu", "s", 1)
plot_fun("mu", "s", 2)
pcre_plots <-plot_fun("mu", "s", 3)
pcre_plots[[1]]; pcre_plots[[2]]; pcre_plots[[3]]; pcre_plots[[4]];
plot_fun("sigma", "p")
plot_fun("alpha", "s", 1)
plot_fun("alpha", "p")

plot(1:200, sapply(it_param, function (x) {
  x$lambda$s[[1]][1]
}))
