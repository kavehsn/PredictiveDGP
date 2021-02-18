#' Normal errors with stationary GARCH(1,1) variance 
#' 
#' 
#' The errors of the predictive regression are distributed according to \eqn{\varepsilon_t\sim N(0,\sigma_t^2)}, 
#' where \eqn{\sigma_t^2=0.00037+0.0888\varepsilon_{t-1}^2+0.9024\sigma_{t-1}^2}. Examples of DGPs with Normal disturbances and stationary GARCH(1,1) variance can be found
#' in \insertCite{dufour2010exact;textual}{Rdpack} and \insertCite{coudin2009finite;textual}{Rdpack}.
#' 
#' @param n the number of observations.
#' @param beta the regressor coefficient of the predictive regression.
#' @param theta the autocorrelation coefficient of the predictor.
#' @param rho the contemporaneous correlation coefficient.    
#' @keywords Endogeneity, persistency, Normal
#' @import pracma
#' @export 
#' @examples
#' GarchDGP(50,0.5,0.999,0.9)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' 
#' 
#' 
GarchDGP<-function(n,beta,theta,rho){

	vari<-rep(0,times=n)
	eps<-rep(0,times=n)

	eps[1]<-rnorm(1,0,1)
	vari[1]<-1

	for(l in 2:n){

		vari[l]<-0.00037+0.0888*(eps[l-1]^2)+0.9024*vari[l-1]
		eps[l]<-rnorm(1,0,vari[l])

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

	z<-flipud(matrix(data=c(y,x),nrow=n,ncol=2))

	return(z)
}