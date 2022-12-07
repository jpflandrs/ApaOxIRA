# ApaOxIRA
#Apatite Oxygen Isotopes interpretation risk analysis
## Aims
Statistics for normality and risk analysis by using a Monte-Carlo simulation
of the alteration of the original oxygen isotope composition of biogenic apatites
## Version 2023  mk03

## Authors 

### JP Flandrois 
jean-pierre.Flandrois@univ-lyon1.fr

Author and maintainer
### C Lecuyer
christophe.lecuyer@univ-lyon1.fr

Original idea and design, equations

## How it works :

The risk of not identifying a diagenesis process by using the Normality criteria (-r) is estimated by an iterative Monte-Carlo simulation. Monte Carlo simulation is based on the creation of an array of A/W values taken from a uniform distribution between given threshold values of δ18O. This represents the diversity of the deposit situation and physical state of the biological apatite. The array is composed of subarrays with length corresponding to the number n of samples in a sample set. Typically 100,000 sets of length n are simulated. The second phase uses equation (4) to compute the δ18O in biological apatite δ18OAf at the equilibrium on each item of the whole array given the temperature T, initial δ18O in biological apatite (δ18OAi) and the δ18O in water δ18OWi. (These two steps may be repeated once, using the array issued from the previous computation (the array of δ18OAf values) as initial value of δ18OA with different values of A/W, δ18OWi and T.) The resulting δ18OAf array is then optionnaly used to simulate the analytical uncertainty by using Normal law random sampling whose mean being the individal δ18OAf and given standard deviation (default 0.05). Finaly each subarray representing a sample set of lenth n of δ18OAf values is submitted to the normality tests. Shapiro-Wilk and the Anderson-Darling are combined. Both have demonstrated a high potential in normality detection, but they are using different approaches and have different sensitivity to skew and tails. Normality is rejected if anyone test reject it at the given risk α. The number of sets in the array were normality cannot be rejected is computed. The result is expressed as a percentage of the number of sample sets. As Normality is unexpected from the construction but may arise by chance (and the shorter the sample, the higher the chance ), this is the MC estimated risk of the sample strategy with sets of length n. Lastly, a bootstrap (BT) estimator of the mean and standard deviation is computed. The whole simulation is done iteratively with increasing values for the length of the sets with given limits for the exploration and the parameters described in a file using the yaml format.   
                         
Statistical (-s) analysis of experimental data: Basic statistics (mean, variance, skew and kurtosis) are computed for a list of data (δ18O, one column only). The normality tests (Shapiro-Wilk, Anderson-Darling) are used and a joined normality test is computed that enable us to detect or reject normality. The interpretation criteria being that of the MC simulation parameters, the results may be compared to the MC simulation. MC simulation performed for the given number of samples and basic parameters give us the _risk_ of getting a normal distributed population by chance after diagenesis. Additionaly an histogram, the adjusted Normal law and the quantile-quantile plot is also returned. Statistical analysis if the data centered on normality is provided as well as graphical outputs (histogram, qq-plot, fitting to the normal distribution)
