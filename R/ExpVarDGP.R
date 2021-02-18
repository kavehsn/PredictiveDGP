#' Normal errors with exponential variance 
#' 
#' 
#' The errors of the predictive regression are distributed according to \eqn{\varepsilon_t\sim N(0,\sigma_t^2)}, 
#' where \eqn{\sigma_t=\exp(0.5t)}. Examples of DGPs with Normal disturbances and exponential variance can be found in
#' \insertCite{dufour2010exact;textual}{Rdpack} and \insertCite{coudin2009finite;textual}{Rdpack}.
#' 
#' @param n the number of observations.
#' @param beta the regressor coefficient of the predictive regression.
#' @param theta the autocorrelation coefficient of the predictor.
#' @param rho the contemporaneous correlation coefficient.    
#' @keywords Endogeneity, persistency, Normal
#' @import pracma
#' @export 
#' @examples
#' ExpVarDGP(50,0.5,0.999,0.9)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt

ExpVarDGP<-function(n,beta,theta,rho){

	stdev<-rep(0,times=n)
	eps<-rep(0,times=n)

	for(l in 1:n){

		stdev[l]<-exp(0.5*l)
		eps[l]<-rnorm(1,0,stdev[l]^2)

	}


	w<-rnorm(n,0,1)

	u<-rep(0,times=n)
	x<-rep(0,times=n)
	y<-rep(0,times=n)


	x[1]<-w[1]/(sqrt(1-theta^2))


	for(i in 1:n){

		u[i]<-(rho*eps[i])+(w[i]*sqrt(1-rho^2))
	}

	for(k in 2:n){

		y[k]<-beta*x[k-1]+eps[k]
		x[k]<-theta*x[k-1]+u[k]

	}

	z<-flipud(cbind(y,x))

	return(z)
}
