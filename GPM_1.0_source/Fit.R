#' @title The Fitting Function of \code{GPM} Package
#' @description Fits a Gaussian process (GP) to a set of simulation data as described in \code{reference 1}. Both the inputs and outputs can be multi-dimensional.
#' @param X Matrix containing the training (aka design or input) data points. The rows and columns of \code{X} denote individual observation settings and input dimension, respectively.
#' @param Y Matrix containing the output (aka response) data points. The rows and columns of \code{Y} denote individual observation responses and output dimension, respectively.
#' @param CovType The covariance function of the GP model. Choices include \code{'G'} (default), \code{'PE'}, \code{'LBG'}, and \code{'LB'}. See the \code{references} for the details.
#' @param MinEig The smallest eigen value that the covariance matrix is allowed to have, which in return determines the appraopriate nugget that should be added to the diagonal elements of the covariance matrix. The default depends on the machine precision and is set to \code{10^(round(log10(.Machine$double.eps))/2)}.
#' @param N_int_opt The number of times the log-likelihood function is optimized with an overly restrictive constraint on the smallest eigen value that the covariance matrix is allowed to have.
#' @param N_final_opt The number of times the optimization problem (i.e., maximizing the log-likelihood) is solved with the \code{L-BFGS-B} algorithm to estimate the (hyper)parameters of the covariance function of the GP model. As \code{L-BFGS-B} is gradient based, and to ensure global optimality, for each time the optimization is started from a different initial point. If much variation is seen in the optimization history (stored in \code{Model$Opt_History}), increase \code{N_final_opt}. Variations usually happen when either the dataset or the number of dimensions are large.
#' @param TraceIt Non-negative integer. If positive, tracing information on the progress of the optimization is \strong{printed}. There are six levels of tracing (see \code{\link{optim}}) and higher values will produce more tracing information.
#' @param MaxIter Maximum number of iterations allowed for each optimization (see \code{\link{optim}}).
#' @param Seed An integer for the random number generator. Use this to make the results reproducible.
#' @param w_min,w_max,Refine_Range To estimate the (hyper)parameters of the covariance function, the feasible range should be defined. \code{w_min} and \code{w_max} are vectors determining, resepectively, the lower and upper bounds and their length depends on the parametric form of the covariance function (see \code{reference 1} for the details). \code{Refine_Range} is a scalar \code{Flag} which determines if the provided (or default) feasible range is refined (as explained in \code{reference 1}) before solving the optimization problem \code{N_final_opt} times. If \code{Refine_Range != 0}, the refinement is done.
#' @param LargeMinEig Overly restrictive \code{MinEig} used for refining the provided/default feasible range. See \code{reference 1} for more details.
#' @param Progress \code{Flag} indicating if the fitting process should be summarized. Set it to \code{!= 0} to turn it on.
#' @import lhs
#' @import randtoolbox
#'
#' @return Model A list containing the following components:
#' \itemize{
#' \item{\code{CovFunc}} {A list containing the type and estimated parameters of the covariance function.}
#' \item{\code{Data}} {A list storing the original (but scaled) data.}
#' \item{\code{Details}} {A list of some parameters (used in prediction) as well as some values reporting the total run-time (\code{cost}), leave-one-out cross-validation error (\code{LOO_CV}), the added nugget (\code{Nug_opt}) for satisfying the constraint (i.e., \code{MinEig}) on the smallest eigen value of the covariance matrix.}
#' \item{\code{Opt_History}} {The optimization history.}
#' \item{\code{Setting}} {The default/provided settings for running the code.}
#' }
#'
#' @references
#' \enumerate{
#' \item T. Kearney, S. Tao, R. Bostanabad, D.W. Apley, and W. Chen (2017). GPM: An R Package for Gaussian Process Modeling of Deterministic Computer Simulations. The Journal of Statistical Software.
#' \item M. Plumlee, D.W. Apley (2016). Lifted Brownian kriging models, Technometrics.
#' }
#' @export
#' @seealso
#' \code{\link[stats]{optim}} for the details on \code{L-BFGS-B} used in optimization.\cr
#' \code{\link[GPM]{Predict}} to use the fitted GP model for prediction.\cr
#' \code{\link[GPM]{Draw}} to plot the response via the fitted model.
#' @examples
#' # 1D example: Fit a model (with default settings) and evaluate the performance
#' # by computing the mean squared error (MSE) in prediction.
#' library(lhs)
#' X <- 5*maximinLHS(10, 1)
#' Y <- 2*sin(2*X) + log(X+1)
#' M <- Fit(X, Y)
#' XF <- matrix(seq(0, 5, length.out = 1000), 1000, 1)
#' YF <- Predict(XF, M)
#' MSE <- mean((YF$YF - (2*sin(2*XF) + log(XF+1)))^2)
#'
#' # 1D example: Fit a model, evaluate the performance, and plot the response
#' # along with 95% prediction interval
#' X <- 10*maximinLHS(25, 1) - 5
#' Y <- X*cos(pi*X)
#' M <- Fit(X, Y)
#' XF <- matrix(seq(-5, 5, length.out = 1000), 1000, 1)
#' YF <- Predict(XF, M)
#' MSE <- mean((YF$YF - (XF*cos(pi*XF)))^2)
#' Draw(M, 1, res = 100)
#'
#' # 2D example: Fit a model, evaluate the performance, and plot the response
#' # surface along with 95% prediction interval
#' X <- 2*maximinLHS(20, 2) - 1
#' Y <- X[, 1]^2 + X[, 2]^2
#' M <- Fit(X, Y, CovType = "PE")
#' XF <- 2*maximinLHS(1000, 2) - 1
#' YF <- Predict(XF, M)
#' MSE <- mean((YF$YF - (XF[, 1]^2 + XF[, 2]^2))^2)
#' library(lattice)
#' Draw(M, c(1, 1), res = 30, PI95=1)
#'
#' # 2D example: Plot the previous model wrt X1 in the [-2, 2]
#' # interval with X2=1
#' Draw(M, c(1, 0), LB = -2, UB = 2, res = 30, PI95=1)
#'
#' \dontrun{
#' # 3D example: Compare the performance of Gaussian ("G") and lifted Browninan
#' # with Gamma=1 ("LBG")
#' X <- 2*maximinLHS(50, 3) - 1
#' Y <- cos(X[, 1]^2) + 2*sin(X[, 2]^2) + X[, 3]^2
#' M_G <- Fit(X, Y)
#' M_LBG <- Fit(X, Y, CovType = "LBG")
#' XF <- 2*maximinLHS(1000, 3) - 1
#' YF_G <- Predict(XF, M_G)
#' YF_LBG <- Predict(XF, M_LBG)
#' MSE_G <- mean((YF_G$YF - (cos(XF[, 1]^2) + 2*sin(XF[, 2]^2) + XF[, 3]^2))^2)
#' MSE_LBG <- mean((YF_LBG$YF - (cos(XF[, 1]^2) + 2*sin(XF[, 2]^2) + XF[, 3]^2))^2)
#'
#' # 3D example: Draw the response in 2D using the M_G model when X3=0
#' Draw(M_G, c(1, 1, 0), PI95 = 0, Values = 0, X1Label = 'Input 1', X2Label = 'Input 2')
#' }

Fit <-  function(X, Y, CovType = 'G', Eps = 10^(seq(-2, -12)),N_opt = 5,TraceIt = 0,
                 MaxIter = 100, Seed = 1, w_min = NULL, Criterion = 'LOO_CV',
                 w_max = NULL, Refine_Range = 1, Progress = 0) {

  set.seed(Seed)

  ## Check the data
  if (Progress != 0){
    cat('*** Checking the inputs ...\n')
  }
  if (missing(X)){
    stop('     X must be provided.')
  }
  if (missing(Y)){
    stop('     Y must be provided.')
  }
  if (is.matrix(X) == FALSE) {
    X <- as.matrix(X)
  }
  if (is.matrix(Y) == FALSE) {
    Y <- as.matrix(Y)
  }
  if (nrow(X)!= nrow(Y)){
    stop('     The number of rows (i.e., observations) in X and Y should match!')
  }
  if (is.character(CovType)){
    CovType <- toupper(CovType)
  } else {
    stop('CovType should be a character (options are: "G", "PE", "LB", and "LBG")')
  }
  if (CovType!= 'G' && CovType!= 'PE' && CovType!= 'LBG' && CovType!='LB'){
    stop('     The type of the covariance function is not determined correctly. Supported functions types are "G", "PE", "LB", and "LBG".')
  }
  if (any(Eps < (10^(round(log10(.Machine$double.eps))+3)))){
    stop(paste('     Increase the smallest member of Eps. The minimum allowable is ', toString(10^(round(log10(.Machine$double.eps))+3))))
  }
  if (is.character(Criterion)){
    if(Criterion != 'LOO_CV' && Criterion != 'NLogL'){
      stop('Criterion should be either LOO_CV or NLogL.')
    }
  } else{
    stop('Criterion should be a string (either LOO_CV or NLogL).')
  }

  k <- nrow(X); p <- ncol(X); q <- ncol(Y);


  w_min_default <- -2
  w_max_default <- 2
  Beta_min_default <- 1e-6
  Beta_max_default <- 1 - 1e-4
  Gamma_min_default <- 1e-6
  Gamma_max_default <- 1

  if (is.null(w_min)){
    if (Progress != 0){
      cat('     The default values are used for initializing the lower bounds (w_min).\n')
    }
    if (CovType == 'G') {w_min = rep(x = w_min_default, times = p)}
    else if (CovType == 'PE') {w_min = c(rep(x = w_min_default, times = p), 1)}
    else if (CovType == 'LBG') {w_min = c(rep(x = w_min_default, times = p), Beta_min_default)}
    else if (CovType == 'LB') {w_min = c(rep(x = w_min_default, times = p), Beta_min_default, Gamma_min_default)}
  } else {
    if (!is.vector(w_min)) w_min = as.vector(w_min)
    if (CovType == 'G' && length(w_min) != p){
      stop(paste('     The size of the provided w_min is incorrect! It should be', toString(p)))}
    else if (CovType == 'PE' && length(w_min) != (1 + p)){
      stop(paste('     The size of the provided w_min is incorrect! It should be', toString(p + 1)))}
    else if (CovType == 'LBG' && length(w_min) != (1 + p)){
      stop(paste('     The size of the provided w_min is incorrect! It should be', toString(p + 1)))}
    else if (CovType == 'LB' && length(w_min) != (2 + p)){
      stop(paste('     The size of the provided w_min is incorrect! It should be', toString(p + 2)))}
    if (any(w_min[1:p] < -10)){
      warning(paste('The first ', toString(p), ' values in w_min are recommended to be larger than -10.'))
    }
  }

  if (is.null(w_max)){
    if (Progress != 0){
      cat('     The default values are used for initializing the upper bounds (w_max).\n')
    }
    if (CovType == 'G') {w_max = rep(x = w_max_default, times = p)}
    else if (CovType == 'PE') {w_max = c(rep(x = w_max_default, times = p), 2)}
    else if (CovType == 'LBG') {w_max = c(rep(x = w_max_default, times = p), Beta_max_default)}
    else if (CovType == 'LB') {w_max = c(rep(x = w_max_default, times = p), Beta_max_default, Gamma_max_default)}
  } else {
    if (!is.vector(w_max)) w_max = as.vector(w_max)
    if (CovType == 'G' && length(w_max) != p){
      stop(paste('     The size of the provided w_max is incorrect! It should be', toString(p)))}
    else if (CovType == 'PE' && length(w_max) != (1 + p)){
      stop(paste('     The size of the provided w_max is incorrect! It should be', toString(p + 1)))}
    else if (CovType == 'LBG' && length(w_max) != (1 + p)){
      stop(paste('     The size of the provided w_max is incorrect! It should be', toString(p + 1)))}
    else if (CovType == 'LB' && length(w_max) != (2 + p)){
      stop(paste('     The size of the provided w_max is incorrect! It should be', toString(p + 2)))}
  }
  if (CovType == 'PE' && (w_min[p + 1] < 1 || w_max[p+ 1] > 2)){
    stop(paste('     The exponent in PE should be between ', toString(c(1, 2)), ' (inclusive).'))
  } else if (CovType == 'LBG' && (w_min[p + 1] <= 0 || w_max[p+ 1] > Beta_max_default)){
    stop(paste('     Beta in LBG should be between ', toString(c(0, 1)), ' (exclusive).'))
  } else if (CovType == 'LB' && (w_min[p+1]<=0 || w_max[p+1]>Beta_max_default || w_min[p+2]<Gamma_min_default || w_max[p+2]>Gamma_max_default)){
    stop('     The provided range for Beta and/or Gamma is not acceptable.')
  }


  Setting <- list(CovType = CovType, Eps = Eps, MaxIter = MaxIter,
                 Seed = Seed, Refine_Range = Refine_Range, N_opt = N_opt,
                 w_min = w_min, w_max = w_max, TraceIt = TraceIt)

  ## Normalize the data
  N_hyperC <- length(w_min)
  Xmin <- apply(X, 2, min)
  Xmax <- apply(X, 2, max)

  XN <- t((t(X)-Xmin)/(Xmax-Xmin))
  Ymin <- apply(Y, 2, min)
  Yrange <- apply(Y, 2, max) - Ymin
  YN <- t((t(Y)-Ymin)/Yrange)
  if (CovType == 'LBG' || CovType == 'LB'){
    XN0 <- XN[1, ]
    YN0 <- YN[1, ]
    XN <- t(t(XN) - XN0)[2:k, ]
    YN <- t(t(YN) - YN0)[2:k, ]
    k <- k - 1
  }
  Fn <- matrix(data = 1, nrow = k, ncol = 1)

  ptm <- proc.time()


  Tries <- length(Eps)
  N_unique_local_opt <- matrix(0, Tries, 1)
  History <- NULL
  #A <- sobol(N_opt, dim = N_hyperC, scrambling = 3, seed = Seed)
  CTRL <- c(trace=0, maxit = MaxIter,  REPORT = 1, lmm = 15, factr = 1e6, pgtol = 1e-8)
  for (i in 1:Tries){
    if (i == 1){
      A <- maximinLHS(N_opt, N_hyperC)
      w_int_i <- t(t(A)*(w_max - w_min) + w_min)
      N_opt_i = N_opt
      Eps_i <- Eps[i]
      w_i <- matrix(0, N_opt_i, N_hyperC)
      NLogL_i <- matrix(0, N_opt_i, 1)
      lb <- w_min
      ub <- w_max
      for (j in 1: N_opt_i){
        temp <- optim(w_int_i[j, ], NLogL, gr = NULL, XN, YN, CovType, Eps_i, Fn, k, q,
                      method = 'L-BFGS-B', lower = lb, upper = ub, control = CTRL)
        w_i[j, ] <- temp$par
        NLogL_i[j] <- temp$value
      }
    }else{
      N_opt_i <- N_unique_local_opt[i-1]
      w_int_i <- History$w[[i-1]]
      Eps_i <- Eps[i]
      w_i <- matrix(0, N_opt_i, N_hyperC)
      NLogL_i <- matrix(0, N_opt_i, 1)
      lb[] <- -10
      ub[] <- 4
      for (j in 1: N_opt_i){
        temp <- optim(w_int_i[j, ], NLogL, gr = NULL, XN, YN, CovType, Eps_i, Fn, k, q,
                      method = 'L-BFGS-B', lower = lb, upper = ub, control = CTRL)
        w_i[j, ] <- temp$par
        NLogL_i[j] <- temp$value
      }
    }


    ID = sort(NLogL_i, index.return=TRUE)$ix
    NLogL_i <- as.matrix(NLogL_i[ID, drop = FALSE])
    w_i <- w_i[ID, , drop = FALSE]


    if (Progress != 0){
      Table <- cbind(apply(round(w_i, 3), 1, toString), round(NLogL_i, 2))
      colnames(Table)<- c('Estimated Hyperparameters', ' NLogL')
      Table <- as.table(Table)
      rownames(Table) <- 1:N_opt_i
      cat('\n\t\tOptimization Results with Eps = ', toString(Eps[i]), ' :\n')
      print(Table, row.names = FALSE)
    }

    rownames(NLogL_i) <- 1:N_opt_i
    NLogL_unique  <-  as.matrix(unique(round(NLogL_i, 2)))
    N_unique_local_opt[i] <- nrow(NLogL_unique)
    w_unique <- as.matrix(as.data.frame(w_i)[c(rownames(NLogL_unique)), ])

    if (Progress != 0){
      Table <- cbind(apply(round(w_unique, 3), 1, toString), round(NLogL_unique, 2))
      colnames(Table)<- c('Estimated Hyperparameters', ' NLogL')
      Table <- as.table(Table)
      rownames(Table) <- 1:nrow(NLogL_unique)
      cat('\n\t\tThe unique local optima are:\n')
      print(Table, row.names = FALSE)
    }

    R <- CovMat(XN, XN, CovType, w_unique[1, ])
    Raw_MinEig <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values)[1]
    if (Raw_MinEig < Eps[i]){
      Nug_opt <- Eps[i] - Raw_MinEig;
      R <- R + diag(x = 1, k, k)*(Eps[i] - Raw_MinEig)
    } else {
      Nug_opt <- 0
    }
    L <- t(chol(R))
    if (CovType == 'PE' || CovType=='G'){
      Rinv_Fn <- solve(t(L), solve(L, Fn))
      Rinv_YN <- solve(t(L), solve(L, YN))
      FnTRinvFn <- t(Fn)%*%Rinv_Fn
      B <- t(Fn)%*%Rinv_YN/FnTRinvFn[1]
      LOO_CV <- matrix(0, q, 1)
      inv_L <- solve(L)
      for (ii in 1:q){
        temp <- solve(diag(diag(t(inv_L)%*%inv_L)), Rinv_YN[, ii] - Rinv_Fn%*%B[ii])
        LOO_CV[ii] <- mean(abs(temp))
      }
    } else {
      Rinv_YN <- solve(t(L), solve(L, YN))
      LOO_CV <- matrix(0, q, 1)
      inv_L <- solve(L)
      Rinv_YN <- as.matrix(Rinv_YN)
      for (ii in 1:q){
        temp <- solve(diag(diag(t(inv_L)%*%inv_L)), Rinv_YN[, ii])
        LOO_CV[ii] <- mean(abs(temp))
      }
    }

    History$Nug_opt[[i]] <- Nug_opt
    History$w[[i]] <- w_unique
    History$LOO_CV[[i]] <- LOO_CV
    History$NLogL[[i]] <- NLogL_i[c(rownames(NLogL_unique)), ]
    History$MinNLogL[i] <- NLogL_i[1]
    if (i >= 2){
      if ((History$Nug_opt[[i]] == History$Nug_opt[[i-1]]) ||
          ( (History$MinNLogL[[i]]>History$MinNLogL[[i-1]])&&(History$LOO_CV[[i]]>History$LOO_CV[[i-1]]) ) ){
        break
      }
    }

  }

  Cost <- proc.time()[3] - ptm[3]
  # Organize the optimization results
  #cat('*** Summary of the optimization results (sorted and rounded):\n')
  if (Criterion == 'LOO_CV'){
    ID <- which.min(History$LOO_CV)
  }else{
    ID <- which.min(History$MinNLogL)
  }
  w <- as.matrix(History$w[[ID]][1, ])

  ## Post-processing
  R <- CovMat(XN, XN, CovType, w)
  Raw_MinEig <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values)[1]
  if (Raw_MinEig < Eps[ID]){
    Nug_opt <- Eps[ID] - Raw_MinEig;
    if (Progress != 0){
      cat('\n*** The best model is fitted via Eps = ', toString(Eps[ID]), '. A nugget of ', Nug_opt, ' has been used.')
    }
    R <- R + diag(x = 1, k, k)*(Eps[ID] - Raw_MinEig)
  } else {
    Nug_opt <- 0
  }

  L <- t(chol(R))
  if (CovType == 'PE' || CovType=='G'){
    Rinv_Fn <- solve(t(L), solve(L, Fn))
    Rinv_YN <- solve(t(L), solve(L, YN))
    FnTRinvFn <- t(Fn)%*%Rinv_Fn
    B <- t(Fn)%*%Rinv_YN/FnTRinvFn[1]
    temp <- YN - Fn%*%B
    Sigma2 <- t(temp)%*%solve(t(L), solve(L, temp))/k
    if (CovType=='G'){
      Theta <- w; Power <- 2;
    } else {
      Theta <- w[1:p]; Power <- w[p+1]
      if (Power > 1.999){
        Power <- 2
      }
    }
    Parameters <- list(FnTRinvFn = FnTRinvFn[1], B = B, Sigma2 = Sigma2, Rinv_Fn = Rinv_Fn,
                      Rinv_YN = Rinv_YN, Theta = Theta, Power = Power)
    LOO_CV <- matrix(0, q, 1)
    inv_L <- solve(L)
    for (i in 1:q){
      temp <- solve(diag(diag(t(inv_L)%*%inv_L)), Rinv_YN[, i] - Rinv_Fn%*%B[i])
      LOO_CV[i] <- mean(abs(temp))
    }
  } else {
    Rinv_YN <- solve(t(L), solve(L, YN))
    Alpha <- t(YN)%*%Rinv_YN/k
    if (CovType == 'LBG'){
      A <- w[1:p]; Beta <- w[p+1]; Gamma <- 1
    } else {
      A <- w[1:p]; Beta <- w[p+1]; Gamma <- w[p+2]
      if (Gamma>0.999) Gamma <- 1
    }
    Parameters <- list(XN0 = XN0, YN0 = YN0, Alpha = Alpha, Rinv_YN = Rinv_YN,
                      A = A, Beta = Beta, Gamma = Gamma)
    LOO_CV <- matrix(0, q, 1)
    inv_L <- solve(L)
    Rinv_YN <- as.matrix(Rinv_YN)
    for (i in 1:q){
      temp <- solve(diag(diag(t(inv_L)%*%inv_L)), Rinv_YN[, i])
      LOO_CV[i] <- mean(abs(temp))
    }
  }
  ## Save the results
  Model <- NULL
  Model$CovFunc <- list('CovType' = CovType, 'Parameters' = Parameters, 'w_min' = w_min, 'w_max' = w_max)
  Model$Data <- list('XN' = XN, 'k' = k, 'Xmin' = Xmin, 'Xmax' = Xmax,  'YN' = YN, 'q' = q, 'Yrange' = Yrange, 'Ymin' = Ymin)
  Model$Details <- list('Fn' = Fn, 'L' = L, 'Raw_MinEig' = Raw_MinEig, 'Nug_opt' = Nug_opt, 'LOO_CV' = LOO_CV, 'Cost' = Cost, 'MinNLogL' = History$MinNLogL[ID])
  Model$History <- History
  Model$Setting <- Setting
  class(Model) <- 'GPM'

  return(Model)
}
