#' @title The Plotting Function of \code{GPM} Package
#' @description Plots the predicted response along with the assocaited uncertainty via the GP model fitted by \code{\link[GPM]{Fit}}. Accepts multi-input and multi-output models. See \code{Arguments} for more details on the options.
#'
#' @param Model The GP model fitted by \code{\link[GPM]{Fit}}.
#' @param Plot_wrt A binary vector of length p where p is the dimension of the inputs in \code{Model}. A maximum (minimum) of \code{2} (\code{1}) elements can be \code{1}. The elemenets set to \code{1}, would correspond to the plotting axes.
#' @param LB,UB Vectors of length \code{sum(Plot_wrt)} indicating the lower and upper bounds used for plotting. The first (second) element corresponds to the first (second) non-zero element of \code{Plot_wrt}.
#' @param Values A vector of length p-\code{sum(Plot_wrt)}. The values are assigned to the variables NOT used in plotting and correspond to the zeros in \code{Plot_wrt}.
#' @param Response_ID A positive integer indicating the response that should be plotted if \code{Model} is multi-response.
#' @param res A positive integer indicating the number of points used in plotting. Higher values will result into smoother plots.
#' @param X1Label A string for the label of axis \code{1}.
#' @param X2Label A string for the label of axis \code{2}, if plotting a surface.
#' @param YLabel A string for the label of the response axis.
#' @param Title A string for the title of the plot.
#' @param PI95 Flag (a scalar) indicating whether the \code{95\%} prediction interval should be plotted. Set it to a non-zero value to turn the flag "on".
#'
#' @return NULL
#'
#' @import lattice
#' @export
#' @references
#' \enumerate{
#' \item T. Kearney, S. Tao, R. Bostanabad, D.W. Apley, and W. Chen (2017). GPM: An R Package for Gaussian Process Modeling of Deterministic Computer Simulations. The Journal of Statistical Software.
#' \item M. Plumlee, D.W. Apley (2016). Lifted Brownian kriging models, Technometrics.
#' }
#' @seealso
#' \code{\link[GPM]{Fit}} to see how a GP model can be fitted to a training dataset.\cr
#' \code{\link[GPM]{Predict}} to use the fitted GP model for prediction.\cr
#' @examples
#' # see the examples in \code{\link[GPM]{Fit}}

Draw <-  function(Model, Plot_wrt, LB = NULL, UB = NULL, Values = NULL,
                Response_ID = NULL, res = 20, X1Label = NULL, X2Label = NULL,
                YLabel = NULL, Title = NULL, PI95 = NULL){

  #if (is.GPM(Model) == FALSE){
  #  stop('The 1st input should be a model of class GMP built by Fit.')
  #}
  p <- ncol(Model$Data$XN)
  q <- ncol(Model$Data$YN)

  if (!is.vector(Plot_wrt)) Plot_wrt <- as.vector(Plot_wrt)
  if (length(Plot_wrt) != p){
    stop(paste('     The size of Plot_wrt is incorrect! It should be', toString(p)))
  } else if (any(Plot_wrt>1) || any(Plot_wrt<0) || any((Plot_wrt<1)*(Plot_wrt>0))){
    stop('The elements of Plot_wrt should be either 1 or 0.');
  } else if (sum(Plot_wrt)>2 || sum(Plot_wrt)<1){
    stop('A maximum (minimum) of 2 (1) elements of Plot_wrt (corresponding to the dimensions used in plotting) should be 1.')
  }
  Axes <- (Plot_wrt==1)
  NotAxes <- (Plot_wrt==0)
  SorL <- sum(Plot_wrt)

  ## range of drawing
  if (is.null(LB)){LB <- Model$Data$Xmin[Axes]}
  LB <- as.vector(LB)
  if (is.null(UB)){UB <- Model$Data$Xmax[Axes]}
  UB <- as.vector(UB)
  if (SorL == 1){
    if (length(LB)!=1 || length(UB)!=1){
      stop('LB and UB need to be a scalar for drawing a 1D plot')
    }
  } else {
    if (length(LB)!=2 || length(UB)!=2){
      stop('LB and UB need to be of size [1, 2] for drawing a 2D surface')
    }
  }
  if (any(LB>UB)){
    stop('The elements of LB should be smaller than the corresponding elements of UB')
  }
  ## settings for the fixed variables
  if ((p - SorL)>0){
    if (is.null(Values)){
      Values <- (Model$Data$Xmin[NotAxes] + Model$Data$Xmax[NotAxes])/2
      print('The fixed inputs, are set to the default values.')
    }
    if (length(Values)!=(p - SorL)){
      stop(paste('The size of Values should be [1, ', toString(p - SorL), '].'))
    }
  } else {
    if (!is.null(Values)){
      warning('The values assigned to "Values" for the variable fixed at plotting, are ignored.')
      Values <- NULL
    }
  }
  if (any(Values < Model$Data$Xmin[NotAxes]) || any(Values > Model$Data$Xmax[NotAxes])){
    warning('(Some of) The values assigned to the fixed variables are beyond the region where the model is fitted.')
  }
  ## Output ID, Resolution, axis labels, title, PI95
  if (q>1){
    if (is.null(Response_ID)){
      Response_ID <- 1
      print('There is more than 1 output in the model. The 1st one is chosen for drawing.')
    } else {
      if (length(Response_ID) != 1){
        stop(paste('Response_ID should be an integer between 1 and ', toString(q), '.'))
      } else if (Response_ID>q || Response_ID<1){
        stop(paste('Response_ID should be an integer between 1 and ', toString(q), '.'))
    }
    }
  } else {
    Response_ID <- 1
  }
  if (!is.null(res)){
    if (length(res)!=1){
      stop('"res" should be a scalar.')
    }
  }
  temp <- which(Axes)
  if (is.null(X1Label)){
    X1Label <- paste('X', toString(temp[1]))
  } else if (!is.character(X1Label)){
    stop('X1Label should be a string.')
  }
  if (is.null(X2Label)){
    X2Label <- paste('X', toString(temp[2]))
  } else if (!is.character(X2Label)){
    stop('X2Label should be a string.')
  }
  if (is.null(YLabel)){
    YLabel <- 'Response'
  } else if (!is.character(YLabel)){
    stop('YLabel should be a string.')
  }
  temp <- NULL; j<-0;
  for (i in which(NotAxes)){
    j <- j + 1
    temp <- paste(temp, '  X', toString(i), '=', toString(round(Values[j]*100)/100))
  }
  if (is.null(Title)){
    Title <- temp
  } else if (!is.character(Title)){
    stop('Title should be a string.')
  }
  if (is.null(PI95)){
    PI95 <- 1
  } else if (length(PI95)!=1){
    stop('PI95 should be a scaler. Set it to a non-zero value to turn the flag "on".')
  } else if (PI95!=0 && PI95!=1){
    PI95 <- 1
  }
  ## set the colors
  mycolors.trans <- rgb(c(0,255,0),
                       c(255,0,0),
                       c(0,0,255),alpha = 50,maxColorValue = 255)

  mycolors <- rgb(c(0,255,0),
                 c(255,0,0),
                 c(0,0,255),maxColorValue = 255)

  ## draw
  if (SorL == 1){
    x1 <- as.matrix(seq(LB, UB, length.out = res), res, 1)
    input <- matrix(0, res, p)
    input[, Axes] <- x1
    if (!is.null(Values)){
      input[, NotAxes]<- t(replicate(res, Values))
    }
    Y <- Predict(input, Model, PI95)
    y <- Y$YF[, Response_ID]
    plot(x=x1, y=y, type = "l", lty = 1, col = "blue", main = Title, xlab = X1Label, ylab = YLabel)
    if (PI95==1){
      yp <- y + 1.96*Y$MSE[, Response_ID];
      ym <- y - 1.96*Y$MSE[, Response_ID];
      lines(x1, yp, type = "l", lty = 2, col = "red")
      lines(x1, ym, type = "l", lty = 2, col = "red")
      legend(x= "topright", legend = c('Prediction', '95% PI'), col = c("blue", "red"),
             lty = c(1, 2))
    } else {
      legend(x = "topright", legend = 'Prediction')
    }
  } else {
    if (PI95==1){
      x1 <- seq(LB[1], UB[2], length.out = res)
      x2 <- seq(LB[1], UB[2], length.out = res)
      g <- expand.grid(x = x1, y = x2, gr = 1:3)
      input <- matrix(0, res^2, p)
      input[, Axes] <- as.matrix(g[1:res^2, -3])
      if (!is.null(Values)){
        input[, NotAxes]<- t(replicate(res^2, Values))
      }
      Y <- Predict(input, Model, PI95)
      y <- matrix(Y$YF[, Response_ID], res^2, 1)
      yp <- y + matrix(Y$MSE[, Response_ID], res^2, 1)
      ym <- y - matrix(Y$MSE[, Response_ID], res^2, 1)
      g$z[1:res^2] <- y
      g$z[(res^2+1):(2*res^2)] <- yp
      g$z[(2*res^2+1):(3*res^2)] <- ym

      wireframe(z ~ x * y, data = g, groups = g$gr, col.groups=mycolors.trans,
                xlab = list(X1Label, rot=30),
                ylab = list(X2Label, rot=-30), zlab = list(YLabel, rot=90),
                main = Title, scales = list(arrows = FALSE),
                key=list(text=list(c("Mean Prediction","95% Upper Bound","95% Lower Bound"),col=mycolors), space = "right",
                         lines=list(lty=c(1,1,1),col=mycolors), border = TRUE),
                par.settings = list(axis.line = list(col = "transparent")))
    } else {
      x1 <- seq(LB[1], UB[2], length.out = res)
      x2 <- seq(LB[1], UB[2], length.out = res)
      g <- expand.grid(x = x1, y = x2)
      input <- matrix(0, res^2, p)
      input[, Axes] <- as.matrix(g)
      if (!is.null(Values)){
        input[, NotAxes]<- t(replicate(res^2, Values))
      }
      Y <- Predict(input, Model, PI95)
      y <- matrix(Y$YF[, Response_ID], res^2, 1)
      g$z <- y

      wireframe(z~x*y, data = g, xlab = list(X1Label, rot=30),
                ylab = list(X2Label, rot=-30), zlab = list(YLabel, rot=90),
                main = Title, scales = list(arrows = FALSE),
                drape = TRUE, colorkey = TRUE)
    }
  }
  }
