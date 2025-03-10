# PP-paper-examples
This repository contains all the code to run the examples for the Poisson process 14C summarisation paper:

Heaton TJ, Bard E & Al-assam S. (to appear) "A New Approach to Radiocarbon Summarisation: Rigorous Identification of Variations/Changepoints in the Occurrence Rate of Radiocarbon Samples using a Poisson Process." _Journal of Archaelogical Science_ (XXXX)

all code uses the R library that is provided alongside the manuscript.  


## Repository Format 
The paper provides a series of analyses. The code to reproduce each of them (and to generate the figures used in the manuscript) can all be found here. We also provide the Dale Guthrie (2006) 14C dates and the outputs.  

### R scripts
All scripts can be found in the `R/` directory. The only scripts which need to be directly runs/sourced are entitled e.g. `001_`, `002_`, ... all others scripts are called within these master scripts.

### Data
The original Dale Guthrie (2006) data is found in the `data/PleistoceneDates/` directory

### Output
The output plots are stored in the `output` directory (which contains relevant subdirectories dependent upon the exmaples run). On running the scripts, the figures in these directories will be overwritten.     

## Specific Analyses

The specific analyses one might wish to run are: 

### Analysis of Late-pleistocene Megafauna (Dale Guthrie, 2006)
This code can be run by sourcing:

*001_Analyse_Guthrie_Pleistocene_c14_Dates.R* 

This will fit the Poisson process model to the c14 dates corresponding to humans, alces, bison and mammoth. For each speciies, it will plot the posterior mean of the Poisson process sample occurrence rate, a histogram of the posterior number of changepoints, and histograms of the posterior locations of those changepoints (conditional on their number). All these plots will be written to the relevant files in `output` - they will not appear in the plotting window within Rstudio.  

In the code, we have specified the prior on the number of changepoints to have a mean of 6. The output plots will be saved in `output/PriorMeanChangepoints6/`.  If you change the prior mean on the number of changepoints, you will need to create a new subdirectory with the relevant name, e.g., `output/PriorMeanChangepoints8/` if you change it to 8. 

### Simulation Study on Artificial Examples
The manuscript provides two simulation studies. To run these, source:

- *002_Simulation_Study_1.R* --- to run the simulation study example 1, a single uniform phase (a Poisson process with two changepoints) 
- *003_Simulation_Study_2.R* --- to run the simulation study example 2, a Poisson process with four changepoints

Again plots of the posterior mean of the sample occurrence rate, histograms of the posterior number of changepoints, and histograms of their locations are generated in `output/SimulationStudy1` and `output/SimulationStudy2` respectively. For simulation study 1, we also create the plots in the SI showing individual posterior realisations of the rate (and the mean rate conditional on a specified number of changepoints)  


### Evidence of SPOD failure
The demonstration of how SPDs fail to provide valid estimates are reproduced by

- 004_SPDFailure.R --- 
- 005_SPD_Bootstrapping_Failure.R --- to run the simulation study example 2, a Poisson process with four changepoints




I recommend substantially updating the README file to include: (1) a citation to the manuscript, and (2) a note that it accompanies a manuscript with the title and authors of the manuscript; (3) a brief description of the files included in the repository and their precise relationship to the figures and tables in the manuscript; (4) the names and version numbers of the key pieces of software used; (5) brief instructions to the user about how to get started working with the project, eg. which file to open first to reproduce the results presented in the paper.










