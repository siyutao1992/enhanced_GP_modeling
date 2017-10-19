#' @title The Function for calculating the Negative Log-Likelehood in \code{GPM} Package
#' @description Calculates the negative log-likelihood (excluding all the constant terms) as described in \code{reference 1}.
#'
#' @param w The vector storing all the (hyper)parameters of the covariance function. The length of \code{w} depends on the \code{CovType}. See \code{reference 1}.
#' @param X Matrix containing the training (aka design or input) data points. The rows and columns of \code{X} denote individual observation settings and input dimension, respectively.
#' @param Y Matrix containing the output (aka response) data points. The rows and columns of \code{Y} denote individual observation responses and output dimension, respectively.
#' @param CovType The covariance function of the GP model. Choices include \code{'G'} (default), \code{'PE'}, \code{'LBG'}, and \code{'LB'}. See the \code{references} for the details.
#' @param MinEig The smallest eigen value that the covariance matrix is allowed to have, which in return determines the appraopriate nugget that should be added to the covariance matrix.
#' @param Fn A matrix on \code{1}'s with \code{nrow(X)} rows and \code{1} column. See \code{reference 1}.
#' @param k Number of observations, \code{nrow(X)}.
#' @param q Number of responses, \code{ncol(Y)}.
#'
#' @return nlogl The negative log-likelihood (excluding all the constant terms). See the \code{references}.
#'
#' @details \code{\link[GPM]{Fit}} calls this function with \emph{scaled} \code{X} and \code{Y}. That is, when the user fits a GP model by calling \code{Fit(X, Y)}, \code{X} and \code{Y} are mapped to the \code{[0, 1]} region and then passed to this function.
#'
#' @note This function is \strong{NOT} exported once GPM package is loaded.
#'
#' @references
#' \enumerate{
#' \item T. Kearney, S. Tao, R. Bostanabad, D.W. Apley, and W. Chen (2017). GPM: An R Package for Gaussian Process Modeling of Deterministic Computer Simulations. The Journal of Statistical Software.
#' \item M. Plumlee, D.W. Apley (2016). Lifted Brownian kriging models, Technometrics.
#' }
#' @seealso
#' \code{\link[GPM]{Fit}} to see how a GP model can be fitted to a training dataset.\cr
#' \code{\link[GPM]{Predict}} to use the fitted GP model for prediction.\cr
#' \code{\link[GPM]{Draw}} to plot the response via the fitted model.
#' @examples
#' # see the examples in \code{\link[GPM]{Fit}}

NLogL <-  function(w, X, Y, CovType, MinEig, Fn, k, q){

R = CovMat(X, X, CovType, w)

Raw_MinEig = sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values)[1]
if (Raw_MinEig < MinEig){
  R = R + diag(x = 1, k, k)*(MinEig - Raw_MinEig)
  #print(MinEig - Raw_MinEig)
}

L = t(chol(R))
if (CovType == 'PE' || CovType=='G'){
  FnTRinvFn = t(Fn)%*%solve(t(L), solve(L, Fn))
  B = t(Fn)%*%solve(t(L), solve(L, Y))/FnTRinvFn[1]
  temp = Y - Fn%*%B
  Sigma2 = t(temp)%*%solve(t(L), solve(L, temp))/k
  nlogl = k*log(det(Sigma2)) + q*2*sum(log(diag(L)))
}else{
  Alpha = t(Y)%*%solve(t(L), solve(L, Y))/k
  nlogl = 2*(log(det(Alpha)) + q*2*sum(log(diag(L)))/k)
}

return(nlogl)
}
