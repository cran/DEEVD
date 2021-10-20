#' DEEVD
#'
#' Density Estimation by Extreme Value Distributions
#'
#' @description Two extreme value distributions; Weibull and Gumbel kernel related functions are presented. \code{\link{Weibull}} and \code{\link{Gumbel}} present estimated values and \code{\link{plot.Weibull}} and \code{\link{plot.Gumbel}} plot the densities. While mean squared error can be calculated by using \code{\link{mse}}. Further, some related data sets are also presented.
#'@author Javaria Ahmad Khan, Atif Akbar.
#'@references \itemize{
#'\item Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
#'\item Khan, J. A., & Akbar, A. 2021. Density Estimation Using Gumbel Kernel Estimator. \emph{Open Journal of Statistics} \strong{11} (2), 319-328.
#'}
"_PACKAGE"
#' Estimate Density Values by Gumbel kernel
#'
#' The Gumbel kernel is developed by Khan and Akbar (2020). They provided evidence that performance of their proposed is better then Weibull kernel especially when data belongs to family of extreme distributions.
#' Gumbel Kernel is
#' \deqn{K_{Gumbel(x, \sqrt{h})}(j)=\frac{1}{\sqrt{h}}exp-\left( \frac{j-x}{\sqrt{h}} +exp\left( \frac{j-x}{\sqrt{h}}\right) \right)}
#' @details In this function, choice of bandwidth, number of grid points and scheme that how these grid points are generated are user based. If any parameter(s) is missing then function used default parameters.
#' But at least \code{x} or \code{k} should be specified otherwise \code{NA} will be produced. If \code{x} is missing then function will generate \code{k} grid points between minimum and maximum values of vector. Similarly, if
#' \code{k} is missing then function consider it same to length of main vector. In case if \code{h} is missing then function used normal scale rule bandwidth for non-normal data and described in Silverman (1986).

#' @param x scheme for generating grid points
#' @param y a numeric vector of positive values.
#' @param k gird points
#' @param h the bandwidth
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Khan, J. A., & Akbar, A. 2021. Density Estimation Using Gumbel Kernel Estimator. \emph{Open Journal of Statistics} \strong{11} (2), 319-328.
#' @seealso For other kernel see \code{\link{Weibull}}. To plot the density by using Gumbel kernel \code{\link{plot.Gumbel}} and to calculate MSE use \code{\link{mse}}.
#'
#' @examples
#' #Data can be simulated or real data
#' ## Number of grid points "k" should be at least equal to the data size.
#' ### If user define the generating scheme of gridpoints than number of gridpoints should
#' ####be equal or greater than "k"
#' ###### otherwise NA will be produduced.
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 100)
#' h <- 2
#' den <- Gumbel(x = xx, y = y, k = 200, h = h)
#'
#' ##If scheme for generating gridpoints is unknown
#' y <- rexp(50, 1)
#' h <- 3
#' den <- Gumbel(y = y, k = 90, h = h)
#'
#' ##If user do not mention the number of grid points
#' y <- rexp(23, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 90)
#'
#' \dontrun{
#' #any bandwidth can be used
#' require(KernSmooth)
#' h <- dpik(y)   #Direct Plug-In Bandwidth
#' den <- Gumbel(x = xx, y = y, h = h)
#' }
#'
#' #if bandwidth is missing
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 100)
#' den <- Gumbel(x = xx, y = y, k = 90)
#'
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}

#' @export
Gumbel <- function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  if(is.null(x))
    x = seq(min(y) + 0.05, max(y), length = k)
  if(is.null(k))
    k = length(y)
  if(is.null(h))
    h = 0.79 * IQR(y) * length(y) ^ (-1 / 5)

  KGumbel <- matrix(rep(0, k * n), ncol = k)

  for(j in 1:k) {
    for(i in 1:n) {

      z <- (y[i] - (x[j])) / sqrt(h)
      KGumbel [i, j] <- (1 / sqrt(h)) * exp( - (z + exp( - (z))))
    }}
  fhat <- colMeans(KGumbel)

  results <- list(x = x,
                  y = fhat)
  class ( results ) <- c('list',
                         'Gumbel')
  results
}

#' Estimated Density Values by Weibull kernel
#'
#' The Weibull kernel is developed by Salha et al. (2014). They used it to nonparametric estimation of the probability density function (pdf) and the hazard rate function for independent and identically distributed (iid) data.
#' Weibull Kernel is
#' \deqn{ K_w\left( x, \frac{1}{h}\right)(t) =\frac{\Gamma(1+h)}{hx}\left[ \frac{t\Gamma(1+h)}{x}\right] ^{\frac{1}{h}-1} exp\left( -\left( \frac{t\Gamma(1+h)}{x}\right) ^\frac{1}{h}\right)}

#' @details see the details in the \code{\link{Gumbel}}
#' @param x scheme for generating grid points
#' @param y a numeric vector of positive values
#' @param k number of gird points
#' @param h the bandwidth
#' @import stats
#' @examples
#' #Data can be simulated or real data
#' ## Number of grid points "k" should be at least equal to the data size.
#' ### If user define the generating scheme of gridpoints than number of gridpoints should
#' ####be equal or greater than "k"
#' ##### otherwise NA will be produduced.
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 100)
#' h <- 2
#' den <- Weibull(x = xx, y = y, k = 200, h = h)
#'
#' ##If scheme for generating gridpoints is unknown
#' y <- rexp(50, 1)
#' h <- 3
#' den <- Weibull(y = y, k = 90, h = h)
#'
#' ##If user do not mention the number of grid points
#' y <- rexp(23, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 90)
#'
#' \dontrun{
#' #any bandwidth can be used
#' require(KernSmooth)
#' h <- dpik(y)
#' den <- Weibull(x = xx, y = y, h = h)
#' }
#'
#' #if bandwidth is missing
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 100)
#' den <- Weibull(x = xx, y = y, k = 90)
#'
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
#' Silverman, B. W. 1986. \emph{Density Estimation}. Chapman & Hall/ CRC, London.
#' @seealso For Gumbel kernel see \code{\link{Gumbel}}. To plot its density see \code{\link{plot.Weibull}} and to calculate MSE \code{\link{mse}}.
#' @export

Weibull <- function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  if(is.null(x))
    x = seq(min(y) + 0.05, max(y), length = k)
  if(is.null(k))
    k = length(y)
  if(is.null(h))
    h = 0.79 * IQR(y) * length(y) ^ (-1 / 5)

  KWeibull <- matrix(rep(0, k * n), ncol = k)
  ###########Weibull Kernel##############
  for(j in 1:k) {
    for(i in 1:n) {
      fn <- gamma(1 + h)
      KWeibull[i, j] <- (fn / (h * x[j])) * ((y[i] * fn) / x[j]) ^ ((1 / h) - 1) * exp( - ((y[i] * fn) / x[j] ^ (1 / h)))
    }}
  fhat <- colMeans(KWeibull)

  results <- list(x = x,
                  y = fhat)
  class ( results ) <- c('list',
                         'Weibull')
  results
}

#' Density Plot by Gumbel kernel
#'
#' Plot kernel density by using Gumbel Kernel.
#' @param x an object of class "Gumbel"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' h <- 1.5
#' xx <- seq(min(y) + 0.05, max(y), length = 200)
#' den <- Gumbel(x = xx, y = y, k = 200, h = h)
#' plot(den, type = "l")
#'
#' ##other details can also be added
#' y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' grid <- seq(min(y) + 0.05, max(y), length = 200)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' gr <- Gumbel(x = grid, y = y, k = 200, h = h)
#' plot(gr, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' ## To add true density along with estimated
#' d1 <- density(y, bw = h)
#' lines(d1, type = "p", col = "red")
#' legend("topright", c("Real Density", "Density by Gumbel Kernel"),
#' col=c("red", "black"), lty=c(1,2))
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Khan, J. A., & Akbar, A. 2021. Density Estimation Using Gumbel Kernel Estimator. \emph{Open Journal of Statistics} \strong{11} (2), 319-328.
#' @seealso For Weibull kernel see \code{\link{plot.Weibull}}. To calculate Gumbel estimated values see \code{\link{Gumbel}} and for MSE \code{\link{mse}}.
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
#' y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' h <- 1.5
#' xx <- seq(min(y) + 0.05, max(y), length = 200)
#' den <- Weibull(x = xx, y = y, k = 200, h = h)
#' plot(den, type = "l")
#'
#' ##other details can also be added
#' y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' grid <- seq(min(y) + 0.05, max(y), length = 200)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' gr <- Weibull(x = grid, y = y, k = 200, h = h)
#' plot(gr, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' ## To add true density along with estimated
#' d1 <- density(y, bw = h)
#' lines(d1, type = "p", col = "green")
#' legend("topright", c("Real Density", "Density by Weibull Kernel"),
#' col=c("green", "black"), lty=c(1,2))
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
#' @seealso For Gumbel kernel see \code{\link{plot.Gumbel}}. To calculate Weibull estimated values see \code{\link{Weibull}} and for
#' MSE use \code{\link{mse}}.
#' @export
plot.Weibull <- function(x,...) {
  plot(x$x, x$y,...)
}

#' Calculate Mean Square Error( MSE) by using Extreme value distributions
#'
#'This function calculates the mean squared error (MSE) by using user specified kernel. But distribution of vector should be Exponential, Gamma, Gumbel, Frechet or Weibull. Any other choice of distribution will result \code{NaN}. This function is simillar to function \code{mse} in \pkg{DELTD}, but here more distributions are available for distribution vector.
#' @param kernel type of kernel which is to be used
#' @param type mention distribution of vector.If exponential distribution then use \code{"Exp"}.
#'     If use gamma distribution then type \code{"Gamma"}.If Gumbel distribution is used with scale=1, then use \code{"Gumbel"}. If Weibull distribution then use \code{"Weibull"}. If use Weibull distribution with scale = 1 then use \code{"Weibull"}.
#'     If use Frechet distribution with scale=1 and shape=1 then use \code{"Frechet"}.
#' @import stats
#' @import evd
#' @author Javaria Ahmad Khan, Atif Akbar
#' @references  \itemize{
#'\item Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
#'\item Khan, J. A., & Akbar, A. 2021. Density Estimation Using Gumbel Kernel Estimator. \emph{Open Journal of Statistics} \strong{11} (2), 319-328.
#'}
#' @examples
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 500)
#' h <- 2
#' gr <- Gumbel(x = xx, y = y, k = 200, h = h)
#' mse(kernel = gr, type = "Exp")
#' ## if distribution is other than mentioned \code{type} is used then NaN will be produced.
#' \dontrun{
#' mse(kernel = gr, type ="Beta")
#' }
#' @return Mean Squared Error (MSE)
#' @export
mse<-function(kernel,type){
  ftrue<-switch(type,
                Exp = dexp(kernel$x, (1 / mean(kernel$x))),
                Gamma = dgamma(kernel$x, (mean(kernel$x) / (var(kernel$x) / mean(kernel$x))), (var(kernel$x) / mean(kernel$x))),
                Gumbel = dgumbel(kernel$x, loc=(mean(kernel$x)-0.5772), scale=1, log = FALSE),
                Weibull = dweibull(kernel$x, ((sd(kernel$x)/mean(kernel$x))^-1.806), scale = 1, log = FALSE),
                Frechet = dfrechet(kernel$x, loc=mean(kernel$x), scale=1, shape=1, log = FALSE)
  )
  return(mean((ftrue-kernel$y)^2))#mse of fhat w.r.t. the true density
}


#' Flood discharge in per second from Rhone river
#'
#' A dataset containing the flood discharge in per second in cubic meter from Rhone river.
#'
#' @format A vector with 23 observations

#' @references  Gumbel, E. J. 1941. The return period of flood flows. \emph{The Annals of Mathematical Statistics}. \strong{12}, 163-190.
"Rhone"
#'

#' Flood discharge in per second from Mississippi river
#'
#' A dataset containing the flood discharge in per second in cubic meter from Mississippi  river.
#'
#' @format A vector with 23 observations

#' @references  Gumbel, E. J. 1941. The return period of flood flows. \emph{The Annals of Mathematical Statistics}. \strong{12}, 163-190.
"Mississippi"
#'
#'
#' Data of control patients of Suicide study
#'
#' A dataset which contains a length of treatment spells (in days) of control patients in suicide study.
#'
#' @format A vector with 86 observations

#' @references  Silverman, B. W. 1986. \emph{Density Estimation}. Chapman & Hall/ CRC, London.
"Suicide"
#'
