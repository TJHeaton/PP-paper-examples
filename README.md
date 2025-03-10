# PP-paper-examples
This repository contains all the code to run the examples for the Poisson process 14C summarisation paper:

Heaton TJ, Bard E & Al-assam S. (to appear) "A New Approach to Radiocarbon Summarisation: Rigorous Identification of Variations/Changepoints in the Occurrence Rate of Radiocarbon Samples using a Poisson Process." _Journal of Archaelogical Science_ (XXXX)

## Repository Format 
The paper provides a series of analyses, the code to reproduce them (and to generate figures used in manuscripts) can be found here. We also provide the Dale Guthrie (2006) 14C dates and the outputs.  

### R scripts
The scripts can be found in the R/ directory. The only scripts which need to be directly runs/sourced are entitled e.g. `001_`, `002_`, ... all others scripts are called within these master scripts.

### Data
The original Dale Guthrie (2006) data is found in the `data/PleistoceneDates/` directory

### Output
The output plots are stored in the `output` directory (which contains relevant subdirectories dependent upon the exmaples run). On running the scripts, the figures in these directories will be overwritten.     

## Specific Analyses

The specific analyses one might wish to run are: 

### Analysis of Late-pleistocene Megafauna (Dale Guthrie, 2006)
This code can be run by sourcing:
*001_Analyse_Guthrie_Pleistocene_c14_Dates.R* 
This will fit the Poisson process model to the humans, alces, bison and mammoth c14 dates and plot the outputs. In the code, I have specified the prior on the number of changepoints to  have a mean of 6. The output plots will be saved in `output/PriorMeanChangepoints6/`.  If you change the prior mean on the number of changepoints you will need to create a new subdirectory with the relevant name, e.g., `output/PriorMeanChangepoints8/` if you change it to 8. 

### Simulation Study



### Evidence of SPOD failure





I recommend substantially updating the README file to include: (1) a citation to the manuscript, and (2) a note that it accompanies a manuscript with the title and authors of the manuscript; (3) a brief description of the files included in the repository and their precise relationship to the figures and tables in the manuscript; (4) the names and version numbers of the key pieces of software used; (5) brief instructions to the user about how to get started working with the project, eg. which file to open first to reproduce the results presented in the paper.










