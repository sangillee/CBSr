# CBSr 1.0.4

R Package for using Cubic Bezier Spline as a function approximator especially in modeling latent utility functions.  
While Cubic Bezier Splines (CBS) are heavily used in the graphics software industry, it can also be used as a flexible  
function approximation tool given the right constraints. The CBS package provides a method to calculate the y value  
from a x value given an appropriately constrained CBS curve. It then uses this method to approximate latent utility  
functions in intertemporal choice and risky choice data.  

CBSr requires rJava package to run (https://CRAN.R-project.org/package=rJava).  
rJava package requires Java development kit (JDK) which needs to be installed on your computer.  
On unix systems, R needs to be reconfigured after JDK installation by using command 'R CMD javareconf'.  
See rJava webpage under installation section for more details (http://www.rforge.net/rJava/).

Please use the package on CRAN as it is the latest stable build.  
Package on GITHUB is a development version that is not released yet.  

PsyArXiv : Lee, Glaze, Bradlow, Kable, (2019) <https://doi.org/10.31234/osf.io/2ugwr>.  
CRAN: https://CRAN.R-project.org/package=CBSr  
GITHUB: https://github.com/sangillee/CBSr
