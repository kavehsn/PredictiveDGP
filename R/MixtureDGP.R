#' Mixture of Cauchy and Normal errors
#' 
#' The processes \eqn{y} and \eqn{x} are generated by a predictive regression of the form \eqn{y_t=\beta x_{t-1}+\varepsilon_t} for \eqn{t=1,\cdots,T},
#' in which the regressors follow an AR(1) process - i.e. \eqn{x_t=\theta x_{t-1}+u_t}.  
#' The predictor's errors are distributed according to \eqn{u_t\sim N(0,1)}, whereas the disturbances 
#' of the predictive regression, \eqn{\varepsilon_t}, are distributed \eqn{\varepsilon_t\sim s_t|\varepsilon_t^C|-(1-s_t)|\varepsilon_t^N|}, 
#' where \eqn{P(s_t=0)=P(s_t=1)=0.5} for all \eqn{t}. An example of a DGP with mixture perturbations can be found in \insertCite{dufour2010exact;textual}{PredictiveDGP}. The initial value of the process \eqn{x} is generated by \eqn{x_0=\frac{w_0}{\sqrt{1-\theta^2}}}, 
#' where \eqn{w_t\sim N(0,1)}. Finally, the contemporaneous correlation between the disturbances \eqn{\varepsilon_t} and \eqn{u_t} is captured by \eqn{\rho\varepsilon_t+w_t\sqrt{1-\rho^2}}.
#' 
#' 
#' 
#' @param n the number of observations
#' @param beta the regressor coefficient of the predictive regression
#' @param theta the autocorrelation coefficient of the predictor
#' @param rho the contemporaneous correlation coefficient    
#' @keywords endogeneity persistency normal mixture
#' @import pracma
#' @export 
#' @examples
#' MixtureDGP(n=50, beta=0.5, theta=0.999, rho=0.9)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt


MixtureDGP <- function(n, beta, theta, rho){
	
	#Setting the probability of a standard normal or Cauchy errors occurrence
	pp <- rbinom(n, 1, 0.5)
	
	#Generate the Cauchy error and standard normal error terms for the predictive regression
	eps_1 <- rnorm(n, 0, 1)
	eps_2 <- rt(n, df=1)

	#Generate the mixture error terms for the predictive regression, and the standard normal errors for the regressor
	eps <- pp*abs(eps_2)-(1-pp)*abs(eps_1)
	w <- rnorm(n, 0, 1)

	#Generate empty variables
	u <- rep(0, times=n)
	x <- rep(0, times=n)
	y <- rep(0,times=n)

	#Generate the initial predictor
	x[1] <- w[1]/(sqrt(1-theta^2))

	#Generate the error terms for the regressors. Note the possibility of feedback from the predictive regression errors to future x.
	for(i in 1:n){

		u[i] <- (rho*eps[i])+(w[i]*sqrt(1-rho^2))
	}
	
	#Generate the X and Y vectors
	for(k in 2:n){

		y[k] <- beta*x[k-1]+eps[k]
		x[k] <- theta*x[k-1]+u[k]

	}

	#Outcome variable
	z <- flipud(cbind(y, x))
	
	#Return output
	return(z)

}
