# logrank
Comparing survival curves of two groups using the log rank test<br/>
Comparison of two survival curves can be done using a statistical
hypothesis test called the log rank test. It is used to test the null
hypothesis that there is no difference between the population survival
curves (i.e. the probability of an event occurring at any time point is
the same for each population). This function uses the Kaplan-Meier
procedure to estimate the survival function (KMPLOT), so if it misses, 
logrank will try to download it from FEX.

Syntax: 	logrank(x1,x2,alpha,censflag)
     
    Inputs:
          X1 and X2 (mandatory)- Nx2 data matrix:
                    (X:,1) = survival time of the i-th subject
                    (X:,2) = censored flag 
                            (0 if not censored; 1 if censored)
          note that if X is a vector, all the flags of the second column
          will be set to 0 (all data are not censored).
          ALPHA (optional) - significance level (default 0.05) 
          CENSFLAG (optional) - Censored Plot flag (default 0). If 0
          censored data will be plotted spreaded on the horizontal
          segment; if 1 they will be plotted at the given time of
          censoring.
    Outputs:
          Kaplan-Meier plot
          Log-rank statistics

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2008). LogRank: Comparing survival curves of two groups using the log rank test
http://www.mathworks.com/matlabcentral/fileexchange/22317
