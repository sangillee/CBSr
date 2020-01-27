# CBS_RC.R
# Arthur Lee. Modified Jan.25.2020
# Function for fitting CBS 1-piece or 2-piece function to Risky choice data.
# choice would be 1 if they chose option 1 or 0 if they chose option 2
# Var should be scaled probability
# Because NlcOptim package does not have a multistart function (and is slower than MATLAB),
# we're using fewer starting points. (But NlcOptim seems at least miles faster than other non-linear constraint optimization packages)

CBS_RC <- function(choice,Amt1,Var1,Amt2,Var2,numpiece){
  if(any(Var1>1) | any(Var1<0) | any(Var2>1) | any(Var2<0) ){stop("prob not within [0 1]")}
  numparam <- 6*numpiece; numfit <- 10*numpiece

  lb <- c(-36,rep(0, 6*numpiece-2)) # lower bounds
  ub <- c(36,rep(1, 6*numpiece-2)) # upper bounds

  if (numpiece == 1){ # active parameters (5): logbeta, x2, x3, y2, y3
    A = NULL; B = NULL; # no linear constraints
    confun = NULL # no non-linear constraints
    x0 = c(0,1/3,2/3,1/3,2/3)
  } else if (numpiece == 2){ # active parameters (11): logbeta, x2,x3,x4,x5,x6, y2,y3,y4,y5,y6
    # linear constraints:
    A = rbind(c(0,1,0,-1,0,0, numeric(5)),  # x2-x4<0
              c(0,0,1,-1,0,0, numeric(5)),  # x3-x4<0
              c(0,0,0,1,-1,0, numeric(5)),  # x4-x5<0
              c(0,0,0,1,0,-1, numeric(5)),  # x4-x6<0
              c(numeric(6), 1,0,-1,0,0),  # y2-y4<0
              c(numeric(6), 0,1,-1,0,0),  # y3-y4<0
              c(numeric(6), 0,0,1,-1,0),  # y4-y5<0
              c(numeric(6), 0,0,1,0,-1))  # y4-y6<0


    B = matrix(numeric(8), nrow=8, ncol=1);
    confun = RCtwopiece_nonlincon
    x0 = c(0,1/6,2/6,3/6,4/6,5/6, 1/6,2/6,3/6,4/6,5/6)
  }
  # optimizer input and options
  funcinput <- list(X = x0, objfun = function(x) RCnegLL(x,Amt1,Var1,Amt2,Var2,choice), confun = confun, A = A, B = B,
                    Aeq = NULL, Beq = NULL, lb = lb, ub = ub, tolX = 1e-04,tolFun = 1e-04, tolCon = 1e-04, maxnFun = 1e+07, maxIter = 4000)

  # fitting
  mdl <- RCfitter(funcinput) # initial fit
  for(i in 2:numfit){
    funcinput$X = RCrandstartpoint(numpiece) # new starting position
    newmdl <- RCfitter(funcinput) # new fit
    if(newmdl$fn < mdl$fn){ mdl <- newmdl } # if new fit is better, swap
  }

  # organizing output
  LL <- -mdl$fn*length(choice)
  xpos <- c(0,mdl$par[2:((numparam+1)/2)],1)
  ypos <- c(0,mdl$par[((numparam+1)/2 +1):numparam],1)
  out <- list("type"=paste("CBS",numpiece,sep=""), "LL" = LL, "numparam" = numparam, "scale" = exp(mdl$par[1]),"xpos"=xpos, "ypos"=ypos, "AUC" = CBSfunc(xpos,ypos),"Origmodel"=mdl)
  return(out)
}

RCfitter <- function(inputlist){
  mdl <- NULL
  try(mdl <- do.call(solnl, inputlist),silent=TRUE) # the NlcOptim package fails now and then for unknown reasons
  while(is.null(mdl)){ # if it's not fit
    inputlist$X = RCrandstartpoint(length(inputlist$X)) # new starting position
    try(mdl <- do.call(solnl, inputlist),silent=TRUE)
  }
  return(mdl)
}

RCnegLL <- function(x,A1,V1,A2,V2,Ch){ # objective function: negative log likelihood
  cutoff <- (length(x)+1)/2
  yhat1 <- CBSfunc(c(0,x[2:cutoff],1), c(0,tail(x,-cutoff),1), V1)
  yhat2 <- CBSfunc(c(0,x[2:cutoff],1), c(0,tail(x,-cutoff),1), V2)
  DV <- A1*yhat1 - A2*yhat2 # diff between utilities
  DV[Ch==0] = -DV[Ch==0] # utility toward choice
  reg = -exp(x[1])*DV # scaling by noise parameter
  logp = -log(1+exp(reg)) # directly calculating logp
  logp[reg>709] = -reg[reg>709]; # log(realmax) is about 709.7827. (e.g., try log(exp(709)) vs. log(exp(710)))
  return(-mean(logp)) # making per-trial LL
}

RCtwopiece_nonlincon <- function(x){
  minhandle = 0.1
  x3 = x[3]; x4 = x[4]; x5 = x[5]; y3 = x[8]; y4 = x[9]; y5 = x[10]
  # non-linear inequalities: 0.1^2 -(x4-x3)^2 -(y4-y3)^2 < 0, 0.1^2-(x5-x4)^2-(y5-y4)^2 < 0
  c_ineq = c(minhandle^2 -(x4-x3)^2 -(y4-y3)^2, minhandle^2 -(x5-x4)^2 -(y5-y4)^2)
  # non-linear equalities: (x4-x3)/(y4-y3) = (x5-x4)/(y5-y4)
  ceq = (x4-x3)*(y5-y4)-(x5-x4)*(y4-y3)
  return(list(c=c_ineq,ceq=ceq))
}

RCrandstartpoint <- function(numpiece){
  # function for finding random starting points that abide the constraints. Mostly because solnl can't seem to start from points that do not abide by constraints
  if(numpiece == 1){
    randpoints = runif(3)
    return(c(0,runif(2),randpoints[randpoints!=min(randpoints)],min(randpoints)))
  } else {
    randx3 = runif(1,min=0.1,max=0.6); randy3 = runif(1,min=0.1,max=0.6) # point no. 3
    randx5 = runif(1,min=randx3+0.2,max=0.9); randy5 = runif(1,min=randy3+0.2,max=0.9) # point no.5
    randx4 = 0.5*randx3+0.5*randx5; randy4 = 0.5*randy3+0.5*randy5 # point 4
    randx2 = runif(1,min=0.05,max=randx4); randy2 = runif(1,min=0.05,max=randy4) # point no.2
    randx6 = runif(1,min=randx4,max = 0.95); randy6 = runif(1,min=randy4,max=0.95)
    return(c(0,randx2,randx3,randx4,randx5,randx6,randy2,randy3,randy4,randy5,randy6))
  }
}
