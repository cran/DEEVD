% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEEVD.R
\name{plot.Weibull}
\alias{plot.Weibull}
\title{Density Plot by Weibull kernel}
\usage{
\method{plot}{Weibull}(x, ...)
}
\arguments{
\item{x}{an object of class "Weibull"}

\item{\dots}{Not presently used in this implementation}
}
\value{
nothing
}
\description{
Plot density by using Weibull Kernel.
}
\examples{
y <- rlnorm(100, meanlog = 0, sdlog = 1)
h <- 1.5
xx <- seq(min(y) + 0.05, max(y), length = 200)
den <- Weibull(x = xx, y = y, k = 200, h = h)
plot(den, type = "l")

##other details can also be added
y <- rlnorm(100, meanlog = 0, sdlog = 1)
grid <- seq(min(y) + 0.05, max(y), length = 200)
h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
gr <- Weibull(x = grid, y = y, k = 200, h = h)
plot(gr, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")

## To add true density along with estimated
d1 <- density(y, bw = h)
lines(d1, type = "p", col = "green")
legend("topright", c("Real Density", "Density by Weibull Kernel"),
col=c("green", "black"), lty=c(1,2))
}
\references{
Salha, R. B., El Shekh Ahmed, H. I., & Alhoubi, I. M. 2014. Hazard Rate Function Estimation Using Weibull Kernel. \emph{Open Journal of Statistics} \strong{4} (08), 650-661.
}
\seealso{
For Gumbel kernel see \code{\link{plot.Gumbel}}. To calculate Weibull estimated values see \code{\link{Weibull}} and for
MSE use \code{\link{mse}}.
}
\author{
Javaria Ahmad Khan, Atif Akbar.
}
