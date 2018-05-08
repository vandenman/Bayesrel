# bayesrel
R-Package "bayesrel" provides both Bayesian and Frequentist reliability estimates

The package provides the most common internal consistency estimates, being: 
    coefficient alpha, guttman's lambda-2 and lambda-6, greatest lower bound and mcdonald's omega-h.
    the frequentist estimates are provided with bootstrapped confidence intervals, 
    the bayesian estimates are provided with credible intervals. The method for the bayesian estimates 
    is sampling from the posterior inverse wishart for the covariance matrix based measures, and 
    gibbs sampling from the joint conditional distributions of a factor model in the case of omega 
    additionally, the bayesian posterior mcmc objects are provided, which allow for 
    further analysis
    BEWARE: this is a work in progress. the basic features all work, yet the help and documentation is widely incomplete.
