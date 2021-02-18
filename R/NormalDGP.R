#' Standard Normal errors
#' 
#' 
#' The original framework initially studied by \insertCite{mankiw1986we;textual}{Rdpack}. 
#' The errors of the predictive regression are distributed according to \eqn{\varepsilon_t\sim N(0,1)}.
#' 
#' @param n the number of observations.
#' @param beta the regressor coefficient of the predictive regression.
#' @param theta the autocorrelation coefficient of the predictor.
#' @param rho the contemporaneous correlation coefficient.    
#' @keywords Endogeneity, persistency, Normal
#' @import pracma
#' @export 
#' @examples
#' NormalDGP(50,0.5,0.999,0.9)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#
NormalDGP<-function(n,beta,theta,rho){

	eps<-rnorm(n,0,1)
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
