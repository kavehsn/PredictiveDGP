#' Normal errors with break in variance 
#' 
#' 
#' The errors of the predictive regression are distributed according to \eqn{\varepsilon_t\sim N(0,1)} for \eqn{t\neq \lceil T/2\rceil} and
#' \eqn{\varepsilon_t\sim 1000 N(0,1)} for \eqn{t= lceil T/2\rceil} respectively. An example of standard Normal disturbances with break in variance can be found in 
#' \insertCite{dufour2010exact;textual}{Rdpack}.\rho\varepsilon_t+w_t\sqrt{1-\rho^2}}.
#' 
#' @param n the number of observations.
#' @param beta the regressor coefficient of the predictive regression.
#' @param theta the autocorrelation coefficient of the predictor.
#' @param rho the contemporaneous correlation coefficient.    
#' @keywords Endogeneity, persistency, Normal
#' @import pracma
#' @export 
#' @examples
#' BIVDGP(50,0.5,0.999,0.9)
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' 
#' 
BIVDGP<-function(n,beta,theta,rho){

	eps<-rnorm(n,0,1)
	brk<-round(n/2)
	eps[brk]<-1000*rnorm(1,0,1)
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