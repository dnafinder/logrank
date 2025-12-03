[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/logrank&file=logrank.m)

# logrank

## 📘 Overview
logrank is a MATLAB function for comparing survival curves of two groups using the classical log-rank test. It is designed for right-censored time-to-event data and operates on Kaplan–Meier survival estimates produced by the companion function kmplot.

The log-rank test evaluates the null hypothesis that there is no difference between two population survival curves; in other words, the hazard (instantaneous event rate) is assumed to be the same in both groups at all times under H₀. Deviations from this assumption are summarized by the Mantel–Haenszel statistic and expressed via a hazard ratio, confidence interval, and p-value.

Internally, logrank:

- calls kmplot to construct survival tables and Kaplan–Meier estimates
- reconstructs the number at risk in each group at each distinct event time
- computes the log-rank statistic
- derives a Mantel–Haenszel hazard ratio and its confidence interval
- reports the z-statistic (with Yates correction) and a two-sided p-value
- produces a combined Kaplan–Meier plot for visual comparison

## ✨ Features
- Accepts:
  - vectors of survival times (all uncensored), or
  - N-by-2 matrices [time, censorFlag] for each group
- Utilizes kmplot for:
  - Kaplan–Meier survival estimation
  - basic hazard rate estimates per group
  - handling of censored values and their graphical display
- Computes:
  - log-rank test statistic (UL)
  - standard error based on the pooled risk structure
  - z-statistic with Yates’ continuity correction
  - two-sided p-value
- Provides a Mantel–Haenszel hazard ratio (HR) and (1 – α)·100% confidence interval, where α is user-specified
- Generates a combined Kaplan–Meier plot showing:
  - survival curve for group 1 (blue)
  - survival curve for group 2 (red)
  - censored observations (black “+” symbols)
- Displays results in clear tabular layout using array2table and cell2table

## 📥 Installation
1. Download or clone the repository:
   https://github.com/dnafinder/logrank

2. Also download the companion function kmplot:
   https://github.com/dnafinder/kmplot

3. Add the folders to your MATLAB path:
      addpath('path_to_logrank_folder')
      addpath('path_to_kmplot_folder')

4. Verify that MATLAB can see both functions:
      which logrank
      which kmplot

If MATLAB returns the paths to logrank.m and kmplot.m, the installation is successful.

## ⚙️ Requirements
- MATLAB (any recent version)
- Statistics and Machine Learning Toolbox (required for:
  - tabulate
  - basic statistical computations)

logrank also requires kmplot to be present on the MATLAB path in order to build Kaplan–Meier survival tables.

## 📈 Usage

Basic comparison of two survival curves (both as N-by-2 matrices):

    % x1(:,1) = survival times for group 1
    % x1(:,2) = censor flag (0 = uncensored, 1 = censored)
    % x2(:,1) = survival times for group 2
    % x2(:,2) = censor flag (0 = uncensored, 1 = censored)
    logrank(x1, x2);

Using vectors of event times (all uncensored in both groups):

    t1 = [2; 3; 5; 7; 9; 10];
    t2 = [1; 4; 6; 8; 12; 15];
    logrank(t1, t2);

Specify a different significance level alpha (affects both test and HR confidence interval):

    logrank(x1, x2, 0.01);   % alpha = 0.01 → 99% CI for HR

Control plotting of censored observations (passed through to kmplot):

    % cflag = 0: spread censored marks along the horizontal segment
    % cflag = 1: place censored marks at the actual censoring time
    logrank(x1, x2, 0.05, 0);
    logrank(x1, x2, 0.05, 1);

## 🔢 Inputs

logrank(x1, x2)  
logrank(x1, x2, alpha)  
logrank(x1, x2, alpha, cflag)

- x1, x2  
  - Type:
    - numeric vector (N×1), or
    - numeric matrix (N×2)
  - Description:
    - If a vector:
      - interpreted as survival times for the corresponding group
      - all observations are internally treated as uncensored
      - converted to [x(:) zeros(N,1)]
    - If N-by-2:
      - x(:,1) = survival times
      - x(:,2) = censoring indicator (0 = uncensored, 1 = censored)

- alpha (optional)  
  - Type: scalar numeric  
  - Default: 0.05  
  - Description:
    - significance level for the log-rank test
    - used for the two-sided hypothesis test and for constructing the (1 – α)·100% confidence interval of the hazard ratio

- cflag (optional)  
  - Type: scalar numeric (0 or 1)  
  - Default: 0  
  - Description:
    - controls how censored observations are plotted, relying on kmplot:
      - 0 → censored marks are spread along the horizontal segment to avoid overlap
      - 1 → censored marks are plotted exactly at their censoring time

## 📤 Outputs
logrank does not return outputs directly; instead, it:

- Produces a combined Kaplan–Meier plot for the two groups
- Prints a set of descriptive tables in the command window including:
  - Hazard rate estimates for each curve (from kmplot)
  - Mantel–Haenszel hazard ratio and its confidence interval
  - Log-rank statistic (UL), its standard error (SE), z-statistic, alpha, two-sided p-value, and a textual decision (Accept Ho / Reject Ho)

## 🧠 Interpretation
- The Kaplan–Meier curves for the two groups show the estimated survival functions S₁(t) and S₂(t).
- The log-rank test compares observed vs expected deaths under the null hypothesis H₀: S₁(t) = S₂(t) for all t (or equivalently, equal hazard functions).
- The Mantel–Haenszel hazard ratio (HR) summarizes the relative hazard across time:
  - HR > 1 suggests a higher event rate in the numerator group (Curve_1)
  - HR < 1 suggests a higher event rate in the denominator group (Curve_2)
- The confidence interval for HR reflects uncertainty at the chosen alpha:
  - If it does not include 1, this indicates a statistically significant difference in hazards at the (1 – α)·100% level.
- The z-statistic incorporates Yates’ continuity correction, and the associated p-value is two-sided:
  - Small p-values (p < alpha) lead to rejection of H₀ → “Reject Ho”
  - Larger p-values lead to “Accept Ho” (i.e., not enough evidence to claim a difference)

## 📌 Notes
- Censoring is assumed to be non-informative and independent of the event process.
- The test is most sensitive to proportional differences in hazards between the two groups; strong time-varying effects may require alternative or complementary methods.
- When providing data as vectors, all observations are treated as uncensored. To properly model censored data, always use the N-by-2 format with explicit censoring flags.
- kmplot must be available on the MATLAB path; otherwise, logrank will stop with an error explaining that kmplot is required.

## 🧪 Example

The following example is the same as in the function help section. Suppose you have two groups stored in x and y, each as an N-by-2 matrix [time, censorFlag]:

    load('logrank.mat')
    logrank(x, y)

A typical Kaplan–Meier plot and log-rank output may resemble (figure and values indicative):

    % Kaplan–Meier plot: group 1 (blue), group 2 (red), censored (+)

    LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS
    ------------------------------------------------------------------------------------------
                       Curve_1     Curve_2 
                       ________    ________
     
        Hazard_rate    0.054452    0.043966
     
        Mantel_Haenszel_Hazard_ratio    Confidence_Interval
        ____________________________    ___________________
     
        2.3016                          1.1452    4.6257   
     
          UL      Standard_error      z       alpha    two_tailed_p_value      Comment  
        ______    ______________    ______    _____    __________________    ___________
     
        6.5723    2.8079            2.1626    0.05     0.030574              'Reject Ho'

![](https://github.com/dnafinder/logrank/blob/master/logrank.jpg)

## 🧾 Citation
If you use logrank in research, analysis, or publications, please cite:

Cardillo G. (2008). LogRank: Comparing survival curves of two groups using the log rank test.  
Available at: https://github.com/dnafinder/logrank

## 👤 Author
Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

## 📄 License
The code is provided as-is, without any explicit warranty.  
Please refer to the repository for licensing details if a LICENSE file is present.  
logrank is distributed under the terms specified in the LICENSE file:
https://github.com/dnafinder/logrank
