# New version of survint_gq with FOR loops --------------------------------

survint_gqFOR <- function (pred = c("lambda", "gamma", "long"), pre_fac,
                           pre_vec = NULL, omega, int_fac = NULL,
                           int_vec = NULL, weights, survtime) {
  
  # Helpful quantities
  nw <- length(weights)
  
  # Set up the matrices for score and hess
  nsubj <- length(survtime)
  p <- if (pred == "gamma") ncol(pre_vec) else ncol(int_vec)
  score_int <- matrix(0, nrow = nsubj, ncol = p)
  hess_int <- matrix(0, nrow = nsubj, ncol = p^2)
  
  switch(pred, 
         "lambda" = {
           
           # Iterate over subjects
           for (ni in seq_len(nsubj)) {
             
             score_i <- 0
             hess_i <- 0
             
             # Iterate over GQ evaluation points
             for (wi in seq_along(weights)) {
               
               tmp <- weights[wi]*omega[(ni-1)*nw + wi]
               hess <- vector(length = p^2)
               # Create crossproduct for each row
               for (pi in seq_len(p)) {
                 for (pi_cross in seq_len(pi)){
                   hess[(pi-1)*p+pi_cross] <- hess[(pi_cross-1)*p+pi] <- 
                     int_vec[(ni-1)*nw + wi, pi] *
                     int_vec[(ni-1)*nw + wi, pi_cross]
                 }
               }
               score_i <- score_i + tmp*int_vec[(ni-1)*nw + wi, ]
               hess_i <- hess_i + tmp*hess
               
             }
             
             # Individual weights
             tmp <- survtime[ni] / 2 * pre_fac[ni]
             score_int[ni,] <- tmp * score_i
             hess_int[ni,] <- tmp * hess_i
             
           }
         },
         "gamma" = {
           
           # Iterate over individuals
           for (ni in seq_len(nsubj)) {
             
             # Integral is a scalar value
             survint_i <- 0
             for (wi in seq_len(nw)) {
               survint_i <- survint_i +
                 omega[(ni - 1)*nw + wi] * weights[wi]
             }
             tmp <-  survtime[ni] / 2 * survint_i * pre_fac[ni] 

             score_int[ni,] <- tmp * pre_vec[ni, ]
             hess_int[ni,] <- tmp * tcrossprod(pre_vec[ni, ])
             
           }
         },
         "long" = {
           nmarker <- nrow(int_vec)/(nsubj*nw)
           if (nmarker %% 1 != 0) {
             stop("Dimensions of longitudinal design matrix do not match.")
           }
           
           # Iterate over subjects
           for (ni in seq_len(nsubj)) {
             
             score_i <- 0
             hess_i <- 0
             
             # Iterate over GQ evaluation points
             for (wi in seq_along(weights)) {
               
               tmp <- weights[wi]*omega[(ni-1)*nw + wi]
              
               score_vec_i <- vector(length = p)
               hess_vec_i <- vector(length = p^2)
               
               # Iterate over markers
               for(mi in seq_len(nmarker)) {
                 
                 score_vec_i <- score_vec_i +
                   int_fac[(mi-1)*nsubj*nw + (ni-1)*nw + wi] * 
                   int_vec[(mi-1)*nsubj*nw + (ni-1)*nw + wi, ]
                 
                 # Create crossproduct for each row
                 for (pi in seq_len(p)) {
                   for (pi_cross in seq_len(pi)){
                     hess_vec_i[(pi-1)*p+pi_cross] <- 
                       hess_vec_i[(pi_cross-1)*p+pi] <- 
                       int_fac[(mi-1)*nsubj*nw + (ni-1)*nw + wi]^2 * 
                       int_vec[(mi-1)*nsubj*nw + (ni-1)*nw + wi, pi] *
                       int_vec[(mi-1)*nsubj*nw + (ni-1)*nw + wi, pi_cross]
                   }
                 }
                 
               }
               score_i <- score_i + tmp*score_vec_i
               hess_i <- hess_i + tmp*hess_vec_i
               
             }
             
             # Individual weights
             tmp <- survtime[ni] / 2 * pre_fac[ni]
             score_int[ni,] <- tmp * score_i
             hess_int[ni,] <- tmp * hess_i
           }
         })
  
  list(score_int = score_int, hess_int = hess_int)
}
