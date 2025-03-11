# PP-paper-examples
This repository contains all the code to run the examples for the Poisson process 14C summarisation paper:

Heaton TJ, Bard E & Al-assam S. (to appear) "A New Approach to Radiocarbon Summarisation: Rigorous Identification of Variations/Changepoints in the Occurrence Rate of Radiocarbon Samples using a Poisson Process." _Journal of Archaeological Science_ (XXXX)

All code uses the `carbondate` R library that is provided to accompany the manuscript. This library can be found on CRAN or on Github at [https://tjheaton.github.io/carbondate/](https://tjheaton.github.io/carbondate/)

The scripts use the current development version 1.0.1.9000 available using `devtools::install_github("TJHeaton/carbondate")` but all the code (except the SI figures showing the individual realisations and conditional mean in simulation study 1) will also run using version 1.0.1 available through CRAN.    


## Repository Format 
The paper provides a series of analyses. The code to reproduce each of them (and to generate the figures used in the manuscript) can all be found here. Additionally, we provide the Dale Guthrie (2006) radiocabron (14C) dates. The output plots generated are also included.  

### R scripts
All scripts can be found in the `R/` directory. The only scripts which need to be directly runs/sourced are entitled e.g. `001_`, `002_`, ... all others scripts are called within these master scripts. These master scripts are all independent of one another (they are self-contained)

### Data
The original Dale Guthrie (2006) data is found in the `data/PleistoceneDates/` directory

### Output
The output plots are stored in the `output` directory (which contains relevant sub-directories dependent upon the examples run). On running the scripts, the figures in these directories will be overwritten if the flag `write_plots_to_file <- TRUE` (otherwise the plots will appear in the current plotting device).     

## Specific Analyses

The specific analyses one might wish to run are described below. All scripts can be run independently of one another. If you want to write the plots to file (this will create the figures in the manuscript) then set the flag `write_plots_to_file <- TRUE` at the start of each script: 

### Analysis of Late-pleistocene Megafauna (Dale Guthrie, 2006)
This code can be run by sourcing:

- *001_Analyse_Guthrie_Pleistocene_c14_Dates.R* 

This will fit the Poisson process model to the c14 dates corresponding to humans, alces, bison and mammoth. For each species, it will plot the posterior mean of the Poisson process sample occurrence rate, a histogram of the posterior number of changepoints, and histograms of the posterior locations of those changepoints (conditional on their number). 

In the code, we have specified the prior on the number of changepoints to have a mean of 6. The output plots will be saved in `output/PriorMeanChangepoints6/`.  If you change the prior mean on the number of changepoints, you will need to create a new sub-directory with the relevant name, e.g., `output/PriorMeanChangepoints8/` if you change it to 8. 

### Simulation Study on Artificial Examples
The manuscript provides two simulation studies. To run these, source:

- *002_Simulation_Study_1.R* --- to run the simulation study example 1, a single uniform phase (a Poisson process with two changepoints) 
- *003_Simulation_Study_2.R* --- to run the simulation study example 2, a Poisson process with four changepoints

Again plots of the posterior mean of the sample occurrence rate, histograms of the posterior number of changepoints, and histograms of their locations are generated in `output/SimulationStudy1` and `output/SimulationStudy2` respectively. For simulation study 1, we also create the plots in the SI showing individual posterior realisations of the rate (and the mean rate conditional on a specified number of changepoints)  


### Evidence of SPD failure
The demonstrations/illustrations of how SPDs fail to provide reliable or accurate estimates are reproduced by

- *004_SPDFailure.R* --- creates plots showing how SPDs are overly variable (Figure 2) and how an SPD of a single determination makes no sense (Figure 3)
- *005_SPD_Bootstrapping_Failure.R* --- creates plot illustrating the failure of bootstrapping and confidence intervals for SPDs (Figure 4) 












