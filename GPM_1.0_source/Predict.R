#' @title The Prediction Function of \code{GPM} Package
#' @description Predicts the reponse(s), associated prediction uncertainties, and gradient(s) of a deterministic computer simulator via the GP model fitted by \code{\link[GPM]{Fit}}.
#'
#' @param XF Matrix containing the locations (settings) where the predictions are desired. The rows and columns of \code{XF} denote individual observation settings and input dimension, respectively.
#' @param Model The GP model fitted by \code{\link[GPM]{Fit}}.
#' @param MSE_on Flag (a scalar) indicating whether the uncertainty (i.e., mean squared error \code{MSE}) associated with prediction of the response(s) should be calculated. Set to a non-zero value to calculate \code{MSE}.
#' @param YgF_on Flag (a scalar) indicating whether the gradient(s) of the response(s) are desired. Set to a non-zero value to calculate the gradient(s). See \code{note} below.
#' @param grad_dim A binary vector of length \code{ncol(XF)}. The gradient of the response(s) will be calculated along the dimensions where the corresponding element of \code{grad_dim} is \code{1}. \code{grad_dim} is ignored if \code{YgF_on!=1}.
#'
#' @return Output A list containing the following components:
#' \itemize{
#' \item{\code{YF}} {A matrix with \code{n} rows (the number of prediction points) and \code{q} columns (the number of responses).}
#' \item{\code{MSE}} {A matrix with \code{n} rows and \code{q} columns where each element represents the prediction uncertainty (i.e., the expected value of the squared difference between the prediction and the true response) associated with the corresponding element in \code{YF}.}
#' \item{\code{YgF}} {An array of size \code{n} by \code{sum{grad_dim}} by \code{p}.}
#' }
#'
#' @export
#' @note
#' \enumerate{
#' \item The gradient(s) can be calculated if \code{CovType='G'} or \code{CovType='LBG'}. If \code{CovType='PE'} or \code{CovType='LB'}, the gradient(s) can only be calculated if \code{Power=2} and \code{Gamma=1}, respectively.
#' \item For efficiency, make sure the inputs are vecotrized and then passed to \code{\link[GPM]{Predict}}. Avoid passing inputs individually in a \code{for} loop.
#' }
#' @references
#' \enumerate{
#' \item T. Kearney, S. Tao, R. Bostanabad, D.W. Apley, and W. Chen (2017). GPM: An R Package for Gaussian Process Modeling of Deterministic Computer Simulations. The Journal of Statistical Software.
#' \item M. Plumlee, D.W. Apley (2016). Lifted Brownian kriging models, Technometrics.
#' }
#' @seealso
#' \code{\link[GPM]{Fit}} to see how a GP model can be fitted to a training dataset.\cr
#' \code{\link[GPM]{Draw}} to plot the response via the fitted model.
#' @examples
#' # see the examples in \code{\link[GPM]{Fit}}

Predict <-  function(XF, Model, MSE_on = 0, YgF_on = 0, grad_dim = rep(1, ncol(XF))){
  #if (is.GPM(Model) == FALSE){
  #  stop('The second input should be a model of class GMP built by Fit.')
  #}
  XN = Model$Data$XN
  p = ncol(XN)
  MSE = NULL
  YgF = NULL
  if (ncol(XF) != p){
    stop('The dimension of XF is not correct!')
  }
  if (length(MSE_on)!=1 || length(YgF_on)!=1){
    stop('MSE_on and YgF_on should be scalar flags. Set them to 1 to turn them "on".')
  }
  CovType = Model$CovFun$CovType
  Ymin = Model$Data$Ymin
  Yrange = Model$Data$Yrange
  k = Model$Data$k
  Xmin = Model$Data$Xmin
  Xmax = Model$Data$Xmax
  q = Model$Data$q
  Fn = Model$Details$Fn
  L = Model$Details$L

  m = nrow(XF)
  Fm = matrix(1, m, 1)
  XFN = t((t(XF)-Xmin)/(Xmax-Xmin))

  if (CovType == 'PE' || CovType=='G'){
    Theta = Model$CovFun$Parameters$Theta
    Power = Model$CovFun$Parameters$Power
    B = Model$CovFun$Parameters$B
    Rinv_YN = Model$CovFun$Parameters$Rinv_YN
    Rinv_Fn = Model$CovFun$Parameters$Rinv_Fn
    FnTRinvFn = Model$CovFun$Parameters$FnTRinvFn
    Sigma2 = Model$CovFun$Parameters$Sigma2
    if (CovType == 'PE'){
      Rxf = CovMat(XN, XFN, CovType, c(Theta, Power))
      if (MSE_on){
        Rf = CovMat(XFN, XFN, CovType, c(Theta, Power))
      }
    } else {
      Rxf = CovMat(XN, XFN, CovType, Theta)
      if (MSE_on){
        Rf = CovMat(XFN, XFN, CovType, Theta)
      }
    }
    YFN = Fm%*%B + t(Rxf)%*%(Rinv_YN - Rinv_Fn%*%B)
    YF = t(t(YFN)*Yrange + Ymin)
    if (MSE_on){
      MSE = Rf - t(Rxf)%*%(solve(t(L), solve(L, Rxf))) +
        t(t(Fm) - t(Fn)%*%(solve(t(L), solve(L, Rxf))))%*%
        (t(Fm) - t(Fn)%*%(solve(t(L), solve(L, Rxf))))/FnTRinvFn
      MSE = diag(kronecker(Sigma2, MSE))
      MSE = t(matrix(MSE, q, m)*Yrange^2)
    }
    if (YgF_on){
      if (Power != 2){
        stop('The gradient can be calculated only if Power == 2.')
      }
      if (!is.vector(grad_dim)) grad_dim = as.vector(grad_dim)
      if (length(grad_dim) != p) {
        stop(paste('grad_dim should be a vector of size (1, ', toString(p), ').'))
      }
      if (any(grad_dim<0) || any(grad_dim>1) || any((grad_dim>0)*(grad_dim<1))){
        stop('The elements of grad_dim should be either 1 or 0.')
      }
      YgF = array(0, c(m, sum(grad_dim), q))
      jj = 1
      for (d in 1:p){
        if (grad_dim[i] > 0){
          XFNd = XFN[, d]
          XNd = XN[, d]
          RxfD = (Power*10^Theta[d])/(Xmax[d] - Xmin[d])*(replicate(m, XNd)-t(replicate(k, XFNd)))*Rxf
          YFN_der = t(RxfD)%*%(Rinv_YN - Rinv_Fn%*%B)
          YgF[, jj, ] = t(t(YFN_der)*Yrange)
          jj = jj + 1
        }
      }
    }
  } else {
    XN0 = Model$CovFun$Parameters$XN0
    XFN = t(t(XFN) - XN0)
    YN0 = Model$CovFun$Parameters$YN0
    Alpha = Model$CovFun$Parameters$Alpha
    Rinv_YN = Model$CovFun$Parameters$Rinv_YN
    A = Model$CovFun$Parameters$A
    Beta = Model$CovFun$Parameters$Beta
    Gamma = Model$CovFun$Parameters$Gamma
    if (CovType == 'LB'){
      Rxf = CovMat(XN, XFN, CovType, c(A, Beta, Gamma))
      if (MSE_on){
        Rf = CovMat(XFN, XFN, CovType, c(A, Beta, Gamma))
      }
    } else {
      Rxf = CovMat(XN, XFN, CovType, c(A, Beta))
      if (MSE_on){
        Rf = CovMat(XFN, XFN, CovType, c(A, Beta))
      }
    }
    YFN = Fm%*%YN0 + t(Rxf)%*%Rinv_YN
    YF = t(t(YFN)*Yrange + Ymin)
    if (MSE_on){
      MSE = Rf - t(Rxf)%*%(solve(t(L), solve(L, Rxf)))
      MSE = diag(kronecker(Alpha, MSE))
      MSE = t(matrix(MSE, q, m)*Yrange^2)
    }
    if (YgF_on){
      A = 10^A - 1
      if (Gamma != 1){
        stop('The gradient can be calculated only if Gamma == 1.')
      }
      YgF = array(0, c(m, sum(grad_dim), q))
      jj = 1
      for (d in 1:p){
        if (grad_dim[d] > 0){
          XFNd = XFN[, d]
          XNd = XN[, d]
          RxfD = matrix(0, k, m)
          if (k >= m){
            for (i in 1:m){
              RxfD[, i] = 2*A[d]*Beta/(Xmax[d] - Xmin[d])*(XFNd[i]*(1 + sum(XFN[i, ]^2*A))^(Beta - 1) -
                          (XFNd[i] - XNd)*(1 + colSums((XFN[i, ] - t(XN))^2*A))^(Beta - 1))
            }
          } else{
            for (i in 1:k){
              RxfD[i, ] = 2*A[d]*Beta/(Xmax[d] - Xmin[d])*(XFNd*(1 + colSums(t(XFN^2)*A))^(Beta - 1) -
                          (XFNd - XNd[i])*(1 + colSums((t(XFN) - XN[i, ])^2*A))^(Beta - 1))
            }
          }

          YFN_der = t(RxfD)%*%Rinv_YN
          YgF[, jj, ] = t(t(YFN_der)*Yrange)
          jj = jj + 1
        }
      }
    }
  }
  Output = list(YF = YF, MSE = MSE, YgF = YgF)
  return(Output)
}
