
# Custom Function to Use MFPCA --------------------------------------------

Predict.matrix.pcre2.random.effect <- function(object, data)
{
  if(is.null(object$xt$mfpc))
    stop("need mfpa object!")
  X <- eval_mfpc(object$xt$mfpc, data[[object$timevar]])
  if(ncol(X) != (length(object$term) - 2))
    stop("check M argument in MFPCA()!")
  X <- data.frame(data[[object$term[1]]], X)
  colnames(X) <- object$term[-length(object$term)]
  
  # Muss man dann dieses X vielleicht noch als data argument in die
  # Predict.matrix.pcre.random.effect()
  # übergeben?
  object$term <- object$term[-length(object$term)]
  X <- refund:::Predict.matrix.pcre.random.effect(object, X)
  return(X)
}