# Simulation Description --------------------------------------------------

library(tidyverse)

# Read in all data sets ---------------------------------------------------

file_location <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                        "berlin.de,share=volkmana.hub/JMbamlss/simulation/scen",
                        "_I_130922/data")
files <- list.files(path = file_location, pattern = "Rdata")
all_dat <- lapply(files, function (dat) {
  load(paste0(file_location, "/", dat))
  out <- d_rirs$data
  rm(d_rirs)
  out
})
all_dat_id <- lapply(files, function (dat) {
  load(paste0(file_location, "/", dat))
  out <- d_rirs$data_short
  rm(d_rirs)
  out
})


# Number of observations --------------------------------------------------

# Total
n_obs <- sapply(all_dat, nrow)
summary(n_obs)

# Number of observations per subject over all data sets
n_obs_id <- lapply(all_dat, function(x) {
  table(x$id)
})
summary(do.call(c, n_obs_id))

# Summary of summaries
n_obs_id <- sapply(all_dat, function(x) {
  summary(as.numeric(table(x$id)))
})
rowMeans(n_obs_id)


# Number of events --------------------------------------------------------

# Mean of number of events
rowMeans(sapply(all_dat_id, function(x) {
  table(x$event)/6
}))

# Mean event rate
summary(sapply(all_dat_id, function(x) {
  mean(x$event[1:150])
}))

# Mean summary of survival times over all observations
rowMeans(sapply(all_dat_id, function (x) {
  summary(x$survtime[1:150])
}))

# Mean summary of survival times of events
rowMeans(sapply(all_dat_id, function (x) {
  summary(x$survtime[1:150][x$event[1:150] == 1])
}))

# Mean over all simulated data
summary(do.call(c, lapply(all_dat_id, function (x) {
  x$survtime[1:150][x$event[1:150] == 1]
})))
