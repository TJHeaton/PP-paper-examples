# PP-paper-examples
This repository contains all the code to run the examples for the Poisson process 14C summarisation paper:

Heaton TJ, Bard E & Al-assam S. (to appear) "A New Approach to Radiocarbon Summarisation: Rigorous Identification of Variations/Changepoints in the Occurrence Rate of Radiocarbon Samples using a Poisson Process." _Journal of Archaeological Science_ (XXXX)

All code uses the `carbondate` R library that is provided to accompany the manuscript. This library can be found on CRAN or on Github at [https://tjheaton.github.io/carbondate/](https://tjheaton.github.io/carbondate/). The specific version of `carbondate` used is 1.1.0.

## Installation of carbondate
The easiest way to install the latest `carbondate` release is via CRAN, by typing the following into your R console:
``` r
install.packages("carbondate")
```
You can alternatively install the development version of carbondate from
[GitHub](https://github.com/) with:
``` r
devtools::install_github("TJHeaton/carbondate")
```
Once you have installed the library with either of the above methods, you need to load it using:
```{r}
library(carbondate)
```

### Other required packages 
To implement the manuscript code, you will also need to have installed the `ggplot2` library using:
``` r
install.packages("ggplot2") # This code used version 3.4.2
```

## Repository Format 
The paper provides a series of analyses. The code to reproduce each of them (and to generate the figures used in the manuscript) can all be found here. Additionally, we provide the Dale Guthrie (2006) radiocarbon (14C) dates. The output plots generated are also included. You will need to ensure you have downloaded all the folders and subfolders.  

### R scripts
All scripts can be found in the `R/` directory. The only scripts which need to be directly runs/sourced are entitled e.g. `001_`, `002_`, ... all others scripts are called within these master scripts. These master scripts are all independent of one another (they are self-contained).

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

In the code, we have specified the prior on the number of changepoints to have a mean of 6. The output plots will be saved in `output/PleistoceneMegafauna/`.  If you change the prior mean on the number of changepoints, the names of the plots in this directory will change automatically, e.g., the plot of the changepoint locations will shift to `FitPP_Alces_Locations_Changepoints_Prior_8_Internal_Changes` if you specify the prior on the number of changepoints to have a mean of 8. 

### Simulation Study on Artificial Examples
The manuscript provides three simulation studies. To run these, source:

- *002_Simulation_Study_1_Single_Uniform.R* --- to run the simulation study example 1, a single uniform phase (a Poisson process with two changepoints).
- *003_Simulation_Study_2_Four_Changepoints.R* --- to run the simulation study example 2, a Poisson process with four changepoints.
- *004_Simulation_Study_3__Exponential_Growth.R* --- to run the simulation study example 3 corresponding to a Poisson process with exponential growth in the occurrence rate

Again plots of the posterior mean of the sample occurrence rate, histograms of the posterior number of changepoints, and histograms of their locations are generated. For simulation study 1, we also create the plots in the SI showing individual posterior realisations of the rate (and the mean rate conditional on a specified number of changepoints).  

Should you set `write_plots_to_file <- TRUE` the plots will be written to `output/SimulationStudy1`, `output/SimulationStudy2` and `output/SimulationStudy3` respectively. Otherwise, they will appear in the plotting window.

### Evidence of SPD failure
The demonstrations/illustrations of how SPDs fail to provide reliable or accurate estimates are reproduced by

- *005_SPDFailure.R* --- creates plots showing how SPDs are overly variable (Figure 2) and how an SPD of a single determination makes no sense (Figure 3)
- *006_SPD_Bootstrapping_Failure.R* --- creates plot illustrating the failure of bootstrapping and confidence intervals for SPDs (Figure 4) 

All plots from these examples will be stored in `output/SPDFailure` should you set `write_plots_to_file <- TRUE`. Otherwise, they will appear in the plotting window.

### Additional Analysis of Late-Pleistocene Equus (Dale Guthrie, 2006)
We have also provided code to analyse the occurrence of c14 samples corresponding to horse (Equus) in Alaska and the Yukon (Dale Guthrie, 2006). This can be run using:

- *007_Analyse_Equus.R* 

This analysis is analogous to the analysis of the other megafauna - it is only not shown in the manuscript for space saving reasons. The code fits the Poisson process model to the Equus samples, plots the posterior mean of the Poisson process sample occurrence rate, a histogram of the posterior number of changepoints, and histograms of the posterior locations of those changepoints (conditional on their number). We have specified the prior on the number of changepoints to have a mean of 6 but this can be altered. 

The output plots will be saved in `output/PleistoceneMegafauna/`. If you change the prior mean on the number of changepoints, the names of the plots in this directory will change automatically, e.g., the plot of the changepoint locations will shift to `FitPP_Equus_Locations_Changepoints_Prior_8_Internal_Changes` if you specify the prior on the number of changepoints to have a mean of 8. 









