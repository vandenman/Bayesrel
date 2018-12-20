# bayesrel
R-Package "bayesrel" provides both Bayesian and frequentist internal consistency estimates

The package provides the most common internal consistency estimates, being: 
    Coefficient Alpha, Guttman's Lambda-2/-4/-6, greatest lower bound and Mcdonald's Omega. 
    
    The Bayesian estimates are provided with credible intervals. The method for the Bayesian estimates 
    is sampling from the posterior inverse wishart for the covariance matrix based measures, and 
    gibbs sampling from the joint conditional distributions of a single factor model in the case of Omega.
    The Bayesian posterior mcmc objects are provided, which allow for further analysis.
    
    The frequentist estimates are provided with bootstrapped confidence intervals. The user can choose between non-parametric or paramteric bootstrap. For Alpha the interval can also be analytic. 
    The frequentist Omega can be calculated using a PFA or a CFA. 
    
    A graphical predictive posterior check can be done for the single factor model. The often used if-item-dropped statistics can be calculated, too. 
    The package also allows for the calculation of the probability of an estimator being bigger than a given threshold.
    There is also the functionality to plot the posterior against the prior for a chosen estimator with cutoffs. In addition to that a plot for the multiple posteriors for the if-item-dropped cases can be done.
    
   
