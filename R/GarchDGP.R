#' Normal errors with stationary GARCH(1,1) variance 
#' 
#' The processes \eqn{y} and \eqn{x} are generated by a predictive regression of the form \eqn{y_t=\beta x_{t-1}+\varepsilon_t} for \eqn{t=1,\cdots,T},
#' in which the regressors follow an AR(1) process - i.e. \eqn{x_t=\theta x_{t-1}+u_t}.  
#' The predictor's errors are distributed according to \eqn{u_t\sim N(0,1)}, whereas the disturbances 
#' of the predictive regression, \eqn{\varepsilon_t}, are distributed \eqn{\varepsilon_t\sim N(0,\sigma_t^2)}, 
#' where \eqn{\sigma_t^2=0.00037+0.0888\varepsilon_{t-1}^2+0.9024\sigma_{t-1}^2}. Examples of DGPs with Normal disturbances and stationary GARCH(1,1) variance can be found
#' in \insertCite{dufour2010exact;textual}{PredictiveDGP} and \insertCite{coudin2009finite;textual}{PredictiveDGP}. The initial value of the process \eqn{x} is generated by \eqn{x_0=\frac{w_0}{\sqrt{1-\theta^2}}}, 
#' where \eqn{w_t\sim N(0,1)}. Finally, the contemporaneous correlation between the disturbances \eqn{\varepsilon_t} and \eqn{u_t} is captured by \eqn{\rho\varepsilon_t+w_t\sqrt{1-\rho^2}}.
#' 
#' @param n the number of observations
#' @param beta the regressor coefficient of the predictive regression
#' @param theta the autocorrelation coefficient of the predictor
#' @param rho the contemporaneous correlation coefficient    
#' @keywords endogeneity persistency normal GARCH
#' @import pracma
#' @export 
#' @examples
#' GarchDGP(n=50, beta=0.5, theta=0.999, rho=0.9)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' 
#' 
#' 
GarchDGP <- function(n, beta, theta, rho){

	#Generate empty variables for the errors terms and their respective variances
	vari <- rep(0, times=n)
	eps <- rep(0, times=n)

	#Generate the initial error term of the predictive regression and its variance
	eps[1] <- rnorm(1, 0, 1)
	vari[1] <- 1

	#Generate the error terms with GARCH(1,1) variance for the predictive regression, and the standard normal errors for the regressor.
	for(l in 2:n){

		vari[l] <- 0.00037+0.0888*(eps[l-1]^2)+0.9024*vari[l-1]
		eps[l] <- rnorm(1,0,vari[l])

	}


	w <- rnorm(n, 0, 1)

	#Generate empty variables
	u <- rep(0, times=n)
	x <- rep(0, times=n)
	y <- rep(0, times=n)

	#Generate the initial predictor
	x[1] <- w[1]/(sqrt(1-theta^2))

	#Generate the error terms for the regressors. Note the possibility of feedback from the predictive regression errors to future x.
	for(i in 1:n){

		u[i] <- (rho*eps[i])+(w[i]*sqrt(1-rho^2))
	}

	#Generate the vectors X and Y
	for(k in 2:n){

		y[k] <- beta*x[k-1]+eps[k]
		x[k] <- theta*x[k-1]+u[k]

	}

	#Outcome variables
	z <- flipud(matrix(data=c(y, x), nrow=n, ncol=2))

	#Return output
	return(z)
}
