#' CBSfunc.R
#'
#' Calculate either the Area Under the Curve (AUC) of a CBS function, or calculate the y coordinates of CBS function given x
#' @param xpos Vector of real numbers of length 1+3n (n = 1, 2, 3, ...), corresponding to Bezier points' x-coordinates of a CBS function
#' @param ypos Vector of real numbers of length 1+3n (n = 1, 2, 3, ...), corresponding to Bezier points' y-coordinates of a CBS function
#' @param x Vector of real numbers, corresponding to x-coordinates of a CBS function. Default value is Null.
#' @param javaswitch Boolean indicator. If true, uses java function to calculate y given x. If false, uses R function to do so. By default, it's true as java function is much faster.
#' @return If x is provided, return y coordinates corresponding to x.
#' If x is not provided, return AUC.
#' @examples
#' CBSfunc(c(0,0.3,0.6,1),c(0.5, 0.2, 0.7, 0.9))
#' CBSfunc(c(0,0.3,0.6,1),c(0.5, 0.2, 0.7, 0.9),seq(0,1,0.1))
#' CBSfunc(c(0,0.3,0.6,1),c(0.5, 0.2, 0.7, 0.9),seq(0,1,0.1),FALSE)

CBSfunc <- function(xpos,ypos,x = NULL,javaswitch = TRUE){
  if (is.null(x)){ #x is not provided. hence calculating AUC
    return(CBSAUC(xpos,ypos))
  } else { # x is provided. hence calculating yhat
    # if javaswitch is true, try calling the java function. If not, use R implementation of numerical approximation, which is 13~14 times slower
    if (javaswitch) { # call the java function
      y <- .jcall("CBScalc", returnSig = "[D","getyhat",xpos,ypos,.jarray(x))
    } else {
      # R implementation of numerical approximation
      y <- rep(NA, length(x))
      for(i in seq(4, length(xpos), 3)) {
        if (i==4){idx = x <=xpos[i]} else {idx = xpos[i-3] < x & x <= xpos[i]}
        # calculate cubic equation coefficients
        a = -xpos[i-3]+3*xpos[i-2]-3*xpos[i-1]+xpos[i]
        b = 3*xpos[i-3]-6*xpos[i-2]+3*xpos[i-1]
        c = -3*xpos[i-3]+3*xpos[i-2]
        d = xpos[i-3]-x[idx]

        t = 0.5 # initial point
        ft = a*t^3 + b*t^2 + c*t+d # f(t)
        delta = 0.25

        while(max(abs(ft))>0.0000001){ # using bisection, because other methods kept finding root outside of 0 1 (e.g., Newton-Raphson, Halley's Method)
          t = t-sign(ft)*delta
          delta = delta/2
          ft = a*t^3 + b*t^2 + c*t+d
        }
        # calculate using deCasteljau's algorithm for numerical stability
        y[idx] = (1-t)*((1-t)*((1-t)*ypos[i-3]+t*ypos[i-2])+t*((1-t)*ypos[i-2]+t*ypos[i-1])) + t*((1-t)*((1-t)*ypos[i-2]+t*ypos[i-1])+t*((1-t)*ypos[i-1]+t*ypos[i]))
      }
    }
    return(y)
  }
}

CBSAUC <- function(xpos,ypos){
  AUC = 0;
  for(i in seq(1,length(xpos)-1,3)){
    AUC = AUC + partialAUC(xpos[i],xpos[i+1],xpos[i+2],xpos[i+3],ypos[i],ypos[i+1],ypos[i+2],ypos[i+3])
  }
  return(AUC)
}

partialAUC <- function(x1,x2,x3,x4,y1,y2,y3,y4){return((6*x2*y1-6*x1*y2-10*x1*y1-3*x1*y3+3*x3*y1-x1*y4-3*x2*y3+3*x3*y2+x4*y1-3*x2*y4+3*x4*y2-6*x3*y4+6*x4*y3+10*x4*y4)/20)}
