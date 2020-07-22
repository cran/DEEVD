#' DEEVD
#'
#' Density Estimation by Extreme Value Distributions
#'
#' @description Weibull and Gumbel kernel related functions are presented. \code{\link{Weibull}} and \code{\link{Gumbel}} present estimated values and \code{\link{plot.Weibull}} and \code{\link{plot.Gumbel}} plot the densities. While \code{\link{mseweibull}} and \code{\link{msegumbel}} calculate the Mean Squared Error (MSE).
#'@author Javaria Ahmad Khan, Atif Akbar.
#'@references \itemize{
#'\item Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
#'\item Khan, J. A.; Akbar, A. Density Estimation by Gumbel Kernel. 2020. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#'}
"_PACKAGE"
#' Estimated Density Values by Gumbel kernel
#'
#' Estimated Kernel density values by using Gumbel Kernel.
#' @details The Gumbel kernel is developed by Khan and Akbar (2020). They provided evidence that performance of their proposed is better then Weibull kernel especially when data belongs to family of extreme distributions.
#' Gumbel Kernel is
#' \deqn{K_{Gumbel(x, \sqrt{h})}(j)=\frac{1}{\sqrt{h}}exp-\left( \frac{j-x}{\sqrt{h}} +exp\left( \frac{j-x}{\sqrt{h}}\right) \right)}

#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' Gumbel(y,200,h)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Khan, J. A.; Akbar, A. Density Estimation by Gumbel Kernel. 2020. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#' @seealso For esitimated values by Weibull kernel see \code{\link{Weibull}}. Further, for plot and MSE by Gumbel kernel see \code{\link{plot.Gumbel}} and \code{\link{msegumbel}}, respectively.
#' @export
Gumbel<-function(y,k,h){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)

  fhat <- rep(0, k)
  KGumbel <- matrix(rep(0, k * n), ncol = k)

  for(j in 1:k) {
    for(i in 1:n) {

      z<-(y[i]-(x[j]))/sqrt(h)
      KGumbel [i, j] <-(1/sqrt(h))*exp(-(z+exp(-(z))))
    }
    fhat[j] <- 1/n * (sum(KGumbel[, j]))
  }
  results <- list(x=x, y=fhat)
  class ( results ) <-c('list', 'Gumbel')
  results
}

#' Estimated Density Values by Weibull kernel
#'
#' Estimated Kernel density values by using Weibull Kernel.
#' @details The Weibull kernel is developed by Salha et al. (2014). They used it to nonparametric estimation of the probability density function (pdf) and the hazard rate function for independent and identically distributed (iid) data.
#' Weibull Kernel is
#' \deqn{ K_w\left( x, \frac{1}{h}\right)(t) =\frac{\Gamma(1+h)}{hx}\left[ \frac{t\Gamma(1+h)}{x}\right] ^{\frac{1}{h}-1} exp\left( -\left( \frac{t\Gamma(1+h)}{x}\right) ^\frac{1}{h}\right)}

#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' Weibull(y,200,h)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
#' @seealso For esitimated values by Gumbel kernel see \code{\link{Gumbel}}. Further, for plot and MSE by Weibull kernel see \code{\link{plot.Weibull}} and \code{\link{mseweibull}}, respectively.
#' @export
Weibull<-function(y,k,h){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)

  fhat <- rep(0, k)
  KWeibull <- matrix(rep(0, k * n), ncol = k)

  for(j in 1:k) {
    for(i in 1:n) {
      fn<-gamma(1+h)
      KWeibull[i, j] <-(fn/(h*x[j]))*((y[i]*fn)/x[j])^((1/h)-1)*exp(-((y[i]*fn)/x[j])^(1/h))
    }
    fhat[j] <- 1/n * (sum(KWeibull[, j]))
  }
  results <- list(x=x, y=fhat)
  class ( results ) <-c('list', 'Weibull')
  results
}

#' Density Plot by Gumbel kernel
#'
#' Plot density by using Gumbel Kernel.
#' @param x an object of class "Gumbel"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den <- Gumbel(y,200,h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Khan, J. A.; Akbar, A. Density Estimation by Gumbel Kernel. 2020. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#' @seealso For Weibull kernel see \code{\link{plot.Weibull}}. To calculate Gumbel estimated values see \code{\link{Gumbel}} and for
#' MSE by using Gumbel Kernel \code{\link{msegumbel}}.
#' @export
plot.Gumbel <- function(x,...) {
  plot(x$x, x$y,...)
}

#' Density Plot by Weibull kernel
#'
#' Plot density by using Weibull Kernel.
#' @param x an object of class "Weibull"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den <- Weibull(y,200,h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
#' @seealso For Gumbel kernel see \code{\link{plot.Gumbel}}. To calculate Weibull estimated values see \code{\link{Weibull}} and for
#' MSE by using Gumbel Kernel \code{\link{mseweibull}}.
#' @export
plot.Weibull <- function(x,...) {
  plot(x$x, x$y,...)
}

#' Calculate Mean Square Error( MSE) when Gumbel kernel is used.
#'
#' Calculate MSE by using Gumbel Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If Gumbel distribution is used scale=1 then use \code{"Gumbel"}.
#'     if use Weibull distribution with scale = 1 then use \code{"Weibull"}.
#'     if use Frechet distribution with scale=1 and shape=1 then use \code{"Frechet"}.
#' @import graphics
#' @import stats
#' @import evd
#' @author Javaria Ahmad Khan, Atif Akbar
#' @references Khan, J. A.; Akbar, A. Density Estimation by Gumbel Kernel. 2020. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#' @seealso For Weibull estimator MSE see \code{\link{mseweibull}}. For density estimation by using Gumbel Kernel \code{\link{plot.Gumbel}} and for estimated values
#' of density \code{\link{Gumbel}}.
#' @examples
#' y<-rweibull(350,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' msegumbel(y,200,h,"Weibull")
#' @return MSE
#' @export
msegumbel<-function(y,k,h,type){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  ftrue<-switch(type,

                Gumbel = dgumbel(x, loc=(mean(x)-0.5772), scale=1, log = FALSE),
                Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE),
                Frechet = dfrechet(x, loc=mean(x), scale=1, shape=1, log = FALSE)
  )

  fhat <- rep(0, k)
  KGumbel <- matrix(rep(0, k * n), ncol = k)

  for(j in 1:k) {
    for(i in 1:n) {

      z<-(y[i]-(x[j]))/sqrt(h)
      KGumbel [i, j] <-(1/sqrt(h))*exp(-(z+exp(-(z))))
    }# 1st loop end
    fhat[j] <- 1/n * (sum(KGumbel[, j]))

  }#2nd loop end
  return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
}#function end

#' Calculate Mean Square Error( MSE) when Weibull kernel is used.
#'
#' Calculate MSE by using Weibull Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If Gumbel distribution is used scale=1 then use \code{"Gumbel"}.
#'     if use Weibull distribution with scale = 1 then use \code{"Weibull"}.
#'     if use Frechet distribution with scale=1 and shape=1 then use \code{"Frechet"}.
#' @import graphics
#' @import stats
#' @import evd
#' @author Javaria Ahmad Khan, Atif Akbar
#' @references Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
#' @seealso For Gumbel estimator MSE see \code{\link{msegumbel}}. For density estimation by using Weibull Kernel \code{\link{plot.Weibull}} and for estimated values
#' of density \code{\link{Weibull}}.
#' @examples
#' y<-rweibull(350,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' mseweibull(y,200,h,"Weibull")
#' @return MSE
#' @export

mseweibull<-function(y,k,h,type){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  ftrue<-switch(type,
                Gumbel = dgumbel(x, loc=(mean(x)-0.5772), scale=1, log = FALSE),
                Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE),
                Frechet = dfrechet(x, loc=mean(x), scale=1, shape=1, log = FALSE)
  )
  fhat <- rep(0, k)
  KWeibull <- matrix(rep(0, k * n), ncol = k)

  for(j in 1:k) {
    for(i in 1:n) {

      fn<-gamma(1+h)
      KWeibull[i, j] <-(fn/(h*x[j]))*((y[i]*fn)/x[j])^((1/h)-1)*exp(-((y[i]*fn)/x[j])^(1/h))
    }# 1st loop end
    fhat[j] <- 1/n * (sum(KWeibull[, j]))

  }#2nd loop end
  return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
}#function end
