#' CBS_ITC.R
#'
#' Fit either a 1-piece or 2-piece CBS latent utility function to binary intertemporal choice data.
#'
#' The data should be consisted of n trials (ideally n > 100) in which a participant made choices between two reward options.
#' At least one of the two options should be delayed such that option 1 is receiving \code{Amt1} in delay1 and option 2 is receiving
#' \code{Amt2} in delay2. If one of the options are immediate, than delay may be 0. If the participant chose option 1, \code{choice}
#' should be 1 and 0 otherwise. Note that before running this function, the delay should be normalized to [0 1]. For example if
#' delay1 ranges from 0 to max_delay1, and delay2 ranges from 20 to max_delay2, one could normalize all delays by \code{normD = max{max_delay1,max_delay2}}.
#' Hence, \code{Var1} should be delay1/normD and \code{Var2} should be delay2/normD. This normalization is not done internally because
#' the choice of the normalization constant affects how the CBS function can fit the data.
#'
#' @param choice Vector of 0s and 1s. 1 if the choice was option 1, 0 if the choice was option 2
#' @param Amt1 Vector of real numbers. Reward amount of choice 1.
#' @param Var1 Vector of real numbers. Delay until the reward of choice 1. Should be normalized within 0 and 1
#' @param Amt2 Vector of real numbers. Reward amount of choice 2.
#' @param Var2 Vector of real numbers. Delay until the reward of choice 2. Should be normalized within 0 and 1
#' @param numpiece Either 1 or 2. Number of CBS pieces to use.
#' @return A list containing the following:
#' \itemize{
#'     \item \code{type}: either 'CBS1' or 'CBS2' depending on the number of pieces
#'     \item \code{LL}: log likelihood of the model
#'     \item \code{numparam}: number of total parameters in the model
#'     \item \code{scale}: scaling factor of the logit model
#'     \item \code{xpos}: x coordinates of the fitted CBS function
#'     \item \code{ypos}: y coordinates of the fitted CBS function
#'     \item \code{AUC}: area under the curve of the fitted CBS function
#'     \item \code{Origmodel} : NlcOptim object of the fitted CBS function. Not intended to be useful, but contains diagnostic information such as convergence condition, iterations, etc.
#' }

# Because NlcOptim package does not have a multistart function (and is slower than MATLAB's fmincon),
# we're using fewer starting points. (But NlcOptim seems at least miles faster than other non-linear constraint optimization packages)

CBS_ITC <- function(choice,Amt1,Var1,Amt2,Var2,numpiece){
  if(any(Var1>1) | any(Var1<0) | any(Var2>1) | any(Var2<0) ){
    normD = max(Var1,Var2); Var1 = Var1/normD; Var2 = Var2/normD # normalizing delay to [0 1] for easier parameter search
  }
  numparam <- 6*numpiece; numfit <- 10*numpiece

  lb <- c(-36,rep(0, numparam-1)) # lower bounds
  ub <- c(36,rep(1, numparam-1)) # upper bounds

  if (numpiece == 1){ # active parameters (6): logbeta, x2, x3, y2, y3, y4
    A = rbind(c(0,0,0,0,-1,1),c(0,0,0,-1,0,1)); B = matrix(numeric(2), nrow=2, ncol=1); # linear constraints: y4-y3<0, y4-y2<0
    confun = NULL # no non-linear constraints
    x0 = c(0,1/3,2/3,2/3,1/3,0.01)
  } else if (numpiece == 2){ # active parameters (12): logbeta, x2,x3,x4,x5,x6, y2,y3,y4,y5,y6,y7
    # linear constraints:
    A = rbind(c(0,1,0,-1,0,0, numeric(6)),  # x2-x4<0
              c(0,0,1,-1,0,0, numeric(6)),  # x3-x4<0
              c(0,0,0,1,-1,0, numeric(6)),  # x4-x5<0
              c(0,0,0,1,0,-1, numeric(6)),  # x4-x6<0
              c(numeric(6), -1,0,1,0,0,0),  # y4-y2<0
              c(numeric(6), 0,-1,1,0,0,0),  # y4-y3<0
              c(numeric(6), 0,0,-1,1,0,0),  # y5-y4<0
              c(numeric(6), 0,0,-1,0,1,0),  # y6-y4<0
              c(numeric(6), 0,0,0,-1,0,1),  # y7-y5<0
              c(numeric(6), 0,0,0,0,-1,1))  # y7-y6<0
    B = matrix(numeric(10), nrow=10, ncol=1);
    confun = ITCtwopiece_nonlincon
    x0 = c(0,1/6,2/6,3/6,4/6,5/6, 5/6,4/6,3/6,2/6,1/6,0.01)
  }
  # optimizer input and options
  funcinput <- list(X = x0, objfun = function(x) ITCnegLL(x,Amt1,Var1,Amt2,Var2,choice), confun = confun, A = A, B = B,
                    Aeq = NULL, Beq = NULL, lb = lb, ub = ub, tolX = 1e-04,tolFun = 1e-04, tolCon = 1e-04, maxnFun = 1e+07, maxIter = 4000)

  # fitting
  mdl <- ITCfitter(funcinput) # initial fit
  for(i in 2:numfit){
    funcinput$X = ITCrandstartpoint(numpiece) # new starting position
    newmdl <- ITCfitter(funcinput) # new fit
    if(newmdl$fn < mdl$fn){ mdl <- newmdl } # if new fit is better, swap
  }

  # organizing output
  LL <- -mdl$fn*length(choice)
  xpos <- c(0,mdl$par[2:(numparam/2)],1)
  ypos <- c(1,mdl$par[(numparam/2 +1):numparam])
  out <- list("type"=paste("CBS",numpiece,sep=""), "LL" = LL, "numparam" = numparam, "scale" = exp(mdl$par[1]),"xpos"=xpos, "ypos"=ypos, "AUC" = CBSfunc(xpos,ypos),"Origmodel"=mdl)
  return(out)
}

ITCfitter <- function(inputlist){
  mdl <- NULL
  try(mdl <- do.call(solnl, inputlist),silent=TRUE) # the NlcOptim package fails now and then for unknown reasons
  while(is.null(mdl)){ # if it's not fit
    inputlist$X = ITCrandstartpoint(length(inputlist$X)) # new starting position
    try(mdl <- do.call(solnl, inputlist),silent=TRUE)
  }
  return(mdl)
}

ITCnegLL <- function(x,A1,V1,A2,V2,Ch){ # objective function: negative log likelihood
  cutoff <- length(x)/2
  yhat1 <- CBSfunc(c(0,x[2:cutoff],1), c(1,tail(x,-cutoff)), V1)
  yhat2 <- CBSfunc(c(0,x[2:cutoff],1), c(1,tail(x,-cutoff)), V2)
  DV <- A1*yhat1 - A2*yhat2 # diff between utilities
  DV[Ch==0] = -DV[Ch==0] # utility toward choice
  reg = -exp(x[1])*DV # scaling by noise parameter
  logp = -log(1+exp(reg)) # directly calculating logp
  logp[reg>709] = -reg[reg>709]; # log(realmax) is about 709.7827. (e.g., try log(exp(709)) vs. log(exp(710)))
  return(-mean(logp)) # making per-trial LL
}

ITCtwopiece_nonlincon <- function(x){
  minhandle = 0.1
  x3 = x[3]; x4 = x[4]; x5 = x[5]; y3 = x[8]; y4 = x[9]; y5 = x[10]
  # non-linear inequalities: 0.1^2 -(x4-x3)^2 -(y4-y3)^2 < 0, 0.1^2-(x5-x4)^2-(y5-y4)^2 < 0
  c_ineq = c(minhandle^2 -(x4-x3)^2 -(y4-y3)^2, minhandle^2 -(x5-x4)^2 -(y5-y4)^2)
  # non-linear equalities: (x4-x3)/(y4-y3) = (x5-x4)/(y5-y4)
  ceq = (x4-x3)*(y5-y4)-(x5-x4)*(y4-y3)
  return(list(c=c_ineq,ceq=ceq))
}

ITCrandstartpoint <- function(numpiece){
  # function for finding random starting points that abide the constraints. Mostly because solnl can't seem to start from points that do not abide by constraints
  if(numpiece == 1){
    randpoints = runif(3)
    return(c(0,runif(2),randpoints[randpoints!=min(randpoints)],min(randpoints)))
  } else {
    randx3 = runif(1,min=0.1,max=0.6); randy3 = runif(1,min=0.4,max=0.9) # point no. 3
    randx5 = runif(1,min=randx3+0.2,max=0.9); randy5 = runif(1,min=0.1,max=randy3-0.2) # point no.5
    randx4 = 0.5*randx3+0.5*randx5; randy4 = 0.5*randy3+0.5*randy5 # point 4
    randx2 = runif(1,min=0.05,max=randx4); randy2 = runif(1,min=randy4,max=0.95) # point no.2
    randy7 = runif(1,min=0,max=randy5) # point no. 7
    randx6 = runif(1,min=randx4,max = 0.95); randy6 = runif(1,min=randy7+0.05,max=randy4)
    return(c(0,randx2,randx3,randx4,randx5,randx6,randy2,randy3,randy4,randy5,randy6,randy7))
  }
}
