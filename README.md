### PROJECT STRUCTURE
.
└── EPHI_paper

    ├── EPHI_paper.Rproj
    
    ├── data_zenodo
    
    ├── scripts
    
        ├── additional scripts
    
        └── numbered scripts
        
    ├── README.md

### DOWNLOADING THE DATA

To allow reproducibility, processed data and results of Bayesian models are available in the subfolder named "data_zenodo" at this URL: https://doi.org/10.5281/zenodo.14418365. This subfolder contains the data used as input for the Bayesian model ran in Script 2, the output of that script as well as the output of robustness simulations (Script 7), to allow reproducing analyses without the need of a High Performance Computing platform. Using these processed data, one can thus start running the series of scripts from Script 3.

Raw data will be available here: https://doi.org/10.5281/zenodo.14185547. They are currently under embargo until the data paper describing them will be published, but can be requested to francois[dot]duchenne[dot]bio[at]gmail.com. Using these raw data, one can start running the series of scripts from the beginning (Script 0).

### R CODES

They are stored in the "scripts" subfolder of the zenodo repository (https://doi.org/10.5281/zenodo.14418365). They are also available at: https://github.com/f-duchenne/EPHI_paper.

To reproducing the analysis, the user will need the "R2jags" R package. Before installing this package the user needs to install the JAGS program, following instructions available here: http://mcmc-jags.sourceforge.net

Scripts are numbered in the order they were used.

To reproduce analyses, open the "EPHI_paper.Rproj" with Rstudio (R studio can be downloaded for free here: "https://posit.co/downloads/"). Opening the R project under R studio will allow the "here" R package (https://cran.r-project.org/web/packages/here/vignettes/here.html) to locate automatically the data folder on your computer. The user should be able to run all the scripts without the need to change any paths, excepting for the scripts that run on HPC (scripts 2 and 7). 

Script 0 merges plant-hummingbird interaction data from different countries in one unique dataset, and merge that dataset with species traits.

Script 1 prepares the data and the Bayesian model (output as a text file) for analyses.

Script 2 fits the Bayesian model, for each country. It was run on a HPC plateform.

Script 3 uses the output of the Bayesian model to analyse how the proportion of zeros (missing interactions) varies with the sampling pressure.

Script 4 extracts the mechanisms (trait complementarity and exploitation barrier) from the model, to do Figure 2.

Script 5 predicts the interaction network for each site.

Script 6 studies how trait distributions (bill and corolla lengths) vary across sites and how it translates to proportion of forbidden links.

Script 7 performs robustness simulations

Script 8 analyses the output of robustness simulations to estimate the effect of forbidden links on robustness.

The subfolder scripts/additional scripts contains a R function that is called by the main scripts (function_to_predict_from_bayesian.r) and older scripts that has been used during the construction of the study but are not used any more in the final version, but could still be useful.
