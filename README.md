# ApaOxIRA
# Apatite Oxygen Isotopes interpretation risk analysis
## Aims

ApaOxIRA enable to choose an optimal sampling scheme that minimize this risk or to abandon the study before any interpretation if the risk in too important.

ApaOxIRA estimate the risk of considering biological apatite δ18O data (essentially in thanatocenoses) being exempt of a diagenesis effect. When used on users data, ApaOxIRA computes statistics for normality and proceed to risk analysis by using a Monte-Carlo simulation for the given sample

## Version 
The current version is mk03 2023  

## Authors 

### JP Flandrois 
jean-pierre.Flandrois@univ-lyon1.fr

Author and maintainer
### C Lecuyer
christophe.lecuyer@univ-lyon1.fr

Original idea and design, equations

## User Licence
The program is made available under the [CeCILL2.1](http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt) licence.

## Rapid description

### Hypothesis

In condition of sedimentary deposits in short time and limited space the distribution of biological apatite δ18O data is expected to be _issued from apopulation following a Normal Law_ (H0). 

If the samples have been submitted to a diagenesis process the distribution of biological apatite δ18O data do not fit to a Normal Law (blurred signal).

### Estimation the risk level linked to sample size
* It is possible to compute (under some hypothesis defined in the companion "yaml" file) a table of the over-interpretation risk related to the sample number. Standard deviation of the estimates is also given (bootstrap).

#### Usage

* The table gives the risk of over-interpretation (due to the observation of the normal distribution "by chance" when the number of samples increases.
* Before any isotopic analysis, this enables to choose a number of sample that lower the risk under a desired level (5% ideally)
* If performed lately this enables to increase the sampling size to fulfill the desired level of risk

### Analysing the  experimental data

When used on users data, ApaOxIRA computes _statistics for normality_ and proceed to a _risk analysis of over-interpretation_ by using a Monte-Carlo simulation for the given sample

#### Interpretation

* If the observed apatite δ18O are _not_ following a Normal distribution, there is (following the H0 definition) either:
    * a diagenesis process
    * the spatio-temporal event that has given the deposit does not correspond to a suffisently limited time or space to allows the Normal Law hypothesis.

* If the biological apatite δ18O data are issued from a Normal distribution:
    * the thanatocenose and the absence of diagenesis is expected
        * but *unfortunately* the Normal distribution _may have been observed only by chance_, and this is heavilly dependant of the number of samples. 
        * ApaOxSis enables to compute the proportion of such "Normal distribution by chance" for a sample of given number within a original population submitted to diagenesis. 
        * It estimate the risk level to _disagree with H0 even if a Normal distribution is observed_ (risk of over-interpretation). This is obtained by the Monte-Carlo process described lower with the sample size as hypothesis.



### Conclusion

ApaOxIRA enables: 

* hopefully to validate the hypothesis 
* or to choose an optimal sampling scheme that minimize this risk (increasing the sampling or limiting the space distribution or both)
* or to abandon the study before any interpretation if the risk in too important for the current number of sample.

## Description of usage:

`` ApaOxIRA [-h] [-r | -x | -s ] i o``

### positional arguments:
*  i                  Parameter file (yaml file, WARNING: structure is mandatory) **or** user data
*  o                  The _directory_ that will gather the results
#
### optional arguments:
*  -h, --help        : Show this help message and exit
*  -r, --Water       : Computes the table linking the number of samples to the risk of over-interpretation
*  -x, --Apatite     : Same as -r but repeated 30 times
*  -s, --Temperature : User data are given, output is made of several files concerning statistics and the estimated level of risk of over-interpretation
#
### exemples

#### Table of over-interpretation Risk with sample size

``python3.11 ApaOxIRA.py -r parameters.yaml resultsFiles``
``python3.11 ApaOxIRA.py -x parameters.yaml resultsFiles``

* parameters.yaml a file containing the parameters of the model and some concerning the computer
* here "resultsFile", is a table linking the number of samples to the risk of over-interpretation.

#### Analysis of experimental data

``python3.11 ApaOxIRA.py -s sampleData.txt resultsName ``

* "sampleData.txt" = one column table of apatite δ18O values
* the outputs are files prefixed by the "resultName".

### Language
ApaOxIS is written in python3 and that the optimal version of python is > 3.10: python3.11 is 10% more rapid than python3.10 here!

## How it works :

The risk of not identifying a diagenesis process by using the Normality criteria (option -r) is estimated by an iterative Monte-Carlo simulation. 

Monte Carlo (MC) simulation is based on the creation of an array of A/W values taken from a uniform distribution between given threshold values of δ18O. This represents the diversity of the deposit situation and physical state of the biological apatite. The array is composed of subarrays with length corresponding to the number n of samples in a sample set. Typically 100,000 sets of length n are simulated. 

The second phase uses equation (4) to compute the δ18O in biological apatite δ18OAf at the equilibrium on each item of the whole array given the temperature T, initial δ18O in biological apatite (δ18OAi) and the δ18O in water δ18OWi. These two steps may be repeated once, using the array issued from the previous computation (the array of δ18OAf values) as initial value of δ18OA with different values of A/W, δ18OWi and T.) 

The resulting δ18OAf array is then optionnaly used to simulate the analytical uncertainty by using Normal law random sampling whose mean being the individal δ18OAf and given standard deviation (default 0.05). 

Finaly each subarray representing a sample set of lenth n of δ18OAf values is submitted to the normality tests. Shapiro-Wilk and the Anderson-Darling are combined. Both have demonstrated a high potential in normality detection, but they are using different approaches and have different sensitivity to skew and tails. Normality is rejected if anyone test reject it at the given risk α. The number of sets in the array were normality cannot be rejected is computed. The result is expressed as a percentage of the number of sample sets. As Normality is unexpected from the construction but may arise by chance (and the shorter the sample, the higher the chance ), this is the MC estimated risk of the sample strategy with sets of length n. Lastly, a bootstrap (BT) estimator of the mean and standard deviation is computed. 

The whole simulation is done iteratively with increasing values for the length of the sets with given limits for the exploration and the parameters described in a file using the yaml format.   
                         
Statistical (option -s) analysis of experimental data: 
    * Basic statistics (mean, variance, skew and kurtosis) are computed for a table of data (δ18O, one column only). 
    * The normality tests (Shapiro-Wilk, Anderson-Darling) are used and a joined normality test is computed that enable us to detect or reject normality. 
    * The interpretation criteria being that of the MC simulation parameters, the results may be compared to the MC simulation. MC simulation performed for the given number of samples and basic parameters give us the _risk_ of getting a normal distributed population by chance after diagenesis. 
    * Additionaly an histogram, the adjusted Normal law and the quantile-quantile plot is also returned. 
