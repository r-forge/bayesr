psi_mat_crossprod <- function(Psi, R = NULL, y = NULL) {
  
  X <- Psi$X[Psi$cp_info[[1]], ]
  ni_obs <- Psi$cp_info[[2]]
  
  if (!is.null(R)) {
    
    # Hesse
    R <- R[Psi$cp_info[[1]]]
    .Call("psi_mat_multiplication", X, ni_obs, R)
    
  } else {
    
    # Score longitudinal part
    y <- y[Psi$cp_info[[1]]]
    .Call("psi_vec_multiplication", X, ni_obs, y)
    
  }
  # Score survival part does not seem to profit from C
  
}

psi_mat_multiplication <- function (X, ni_obs, diags){
  p <- ncol(X)
  n <- nrow(X)
  out <- vector(length = p)
  col_it <- 0
  for(i in seq_len(p)) {
    out_i <- 0
    for (j in seq_len(ni_obs[i])) {
      out_i <- out_i + X[(i-1)*n + col_it + j]^2 * diags[col_it + j]
    }
    col_it <- col_it + ni_obs[i]
    out[i] <- out_i
  }
  out
}

psi_vec_multiplication <- function (X, ni_obs, y){
  p <- ncol(X)
  n <- nrow(X)
  out <- vector(length = p)
  col_it <- 0
  for(i in seq_len(p)) {
    out_i <- 0
    for (j in seq_len(ni_obs[i])) {
      out_i <- out_i + X[(i-1)*n + col_it + j] * y[col_it + j]
    }
    col_it <- col_it + ni_obs[i]
    out[i] <- out_i
  }
  out
}