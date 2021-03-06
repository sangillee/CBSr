% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CBS_ITC.R
\name{CBS_ITC}
\alias{CBS_ITC}
\title{CBS_ITC}
\usage{
CBS_ITC(choice, Amt1, Delay1, Amt2, Delay2, numpiece, numfit = NULL)
}
\arguments{
\item{choice}{Vector of 0s and 1s. 1 if the choice was option 1, 0 if the choice was option 2.}

\item{Amt1}{Vector of positive real numbers. Reward amount of choice 1.}

\item{Delay1}{Vector of positive real numbers. Delay until the reward of choice 1.}

\item{Amt2}{Vector of positive real numbers. Reward amount of choice 2.}

\item{Delay2}{Vector of positive real numbers. Delay until the reward of choice 2.}

\item{numpiece}{Either 1 or 2. Number of CBS pieces to use.}

\item{numfit}{Number of model fits to perform from different starting points. If not provided, numfit = 10*numpiece}
}
\value{
A list containing the following:
\itemize{
    \item \code{type}: either 'CBS1' or 'CBS2' depending on the number of pieces
    \item \code{LL}: log likelihood of the model
    \item \code{numparam}: number of total parameters in the model
    \item \code{scale}: scaling factor of the logit model
    \item \code{xpos}: x coordinates of the fitted CBS function
    \item \code{ypos}: y coordinates of the fitted CBS function
    \item \code{AUC}: area under the curve of the fitted CBS function. Normalized to be between 0 and 1.
    \item \code{normD} : The domain of CBS function runs from 0 to \code{normD}. Specifically, this is the constant used to normalize all delays between 0 and 1, since CBS is fitted in a unit square first and then scaled up.
}
}
\description{
Fit either a 1-piece or 2-piece CBS latent utility function to binary intertemporal choice data.
}
\details{
The input data has n choices (ideally n > 100) between two reward options.
Option 1 is receiving \code{Amt1} in \code{Delay1} and Option 2 is receiving \code{Amt2} in \code{Delay2} (e.g., $40 in 20 days vs. $20 in 3 days).
One of the two options may be immediate (i.e., delay = 0; e.g., $40 in 20 days vs. $20 today).
\code{choice} should be 1 if option 1 is chosen, 0 if option 2 is chosen.
}
\examples{
# Fit example ITC data with 2-piece CBS function.
# Load example data (included with package).
# Each row is a choice between option 1 (Amt at Delay) vs option 2 (20 now).
Amount1 = ITCdat$Amt1
Delay1 = ITCdat$Delay1
Amount2 = 20
Delay2 = 0
Choice = ITCdat$Choice

# Fit the model
out = CBS_ITC(Choice,Amount1,Delay1,Amount2,Delay2,2)

# Plot the choices (x = Delay, y = relative amount : 20 / delayed amount)
plot(Delay1[Choice==1],20/Amount1[Choice==1],type = 'p',col="blue",xlim=c(0, 180), ylim=c(0, 1))
points(Delay1[Choice==0],20/Amount1[Choice==0],type = 'p',col="red")

# Plot the fitted CBS
x = 0:out$normD
lines(x,CBSfunc(out$xpos,out$ypos,x),col="black")
}
