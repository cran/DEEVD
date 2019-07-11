#' Plot Density by Gumbel kernel.
#'
#' Plot Kernel density by using Gumbel Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @return Plot the density by using Gumbel Kernel
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' figGum(y,80,h)
#' @export
figGum<-function(y,k,h){
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
    d1<-density(y,bw=h)
    plot(x,fhat, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
    lines(d1,type="p",col="red")
    legend("topright", c("Real Density", "Density by  Gumbel  Kernel"),
           col=c("red", "black"), lty=c(1,2))
  }}

#' Plot Density by Weibull kernel.
#'
#' Plot Kernel density by using weibull Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @return Plot the density by using Weibull Kernel
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' figWeibull(y,80,h)
#' @export
#
figWeibull<-function(y,k,h){
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
    d1<-density(y,bw=h)
    plot(x,fhat, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
    lines(d1,type="p",col="red")
    legend("topright", c("Real Density", "Density by Weibull Kernel"),
           col=c("red", "black"), lty=c(1,2))
  }}


#' Calculate Mean Square Error( MSE) when Gumbel kernel is used.
#'
#' Calculate MSE by using Gumbel Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If Gumbeldistribution is used scale=1 then use "Gumbel".
#'     if use Weibull distribution with scale = 1 then use "Weibull".
#'     if use Frechet distribution with scale=1 and shape=1 then use "Frechet".
#' @import graphics
#' @import stats
#' @import evd
#' @examples
#' y<-rweibull(350,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' msegum(y,200,h,"Weibull")
#' @return MSE
#' @export
msegum<-function(y,k,h,type){
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
#' @param type mention distribution of vector.If Gumbeldistribution is used scale=1 then use "Gumbel".
#'     if use Weibull distribution with scale = 1 then use "Weibull".
#'     if use Frechet distribution with scale=1 and shape=1 then use "Frechet".
#' @import graphics
#' @import stats
#' @import evd
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
