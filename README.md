# Data

The data needed to reproduce analyses are available here: 

# Code

Scripts are labelled in the order they were used.

Script 0 merges plant-hummingbird interaction data from different countries in one unique dataset, and merge that dataset with species traits.

Script 1 prepares the data and the Bayesian model (output as a text file) for analyses.

Script 2 fits the Bayesian model, for each country. It was run on a HPC plateform.

Script 3 used the output of the model to analyse how the proportion of zeros (missing interactions) varies with the sampling pressure.

Script 4 extracts the mechanisms (trait complementarity and exploitation barrier) from the model, to do Figure 2.

Script 5 was used to predict the interaction network for each site.

Script 6 studies how trait distributions (bill and corolla lengths) vary across sites and how it translates to proportion of forbidden links.

Script 7 performs robustness simulations

Script 8 analyses the output of robustness simulations to estimate the effect of forbidden links on robustness.
