# This file analyses the Guthrie (2006) data of Pleistocene Megafauna:
# Specifically:
# Fits Poisson Process model to each species and creates plots/figures shown in real-life case study
# The plots from the PP are merged into Figs 9 and 10

# Load carbondate library
library(carbondate)

# Load other libraries
library(ggplot2)

# Decide if write plots to file
write_plots_to_file <- FALSE

# Decide if want plotting output as pdf (TRUE) or png (FALSE)
pdf_output <- FALSE

# Store user par parameters (and revert back to these at the end)
oldpar <- par(no.readonly = TRUE)

# 14C ages at which to cutoff analysis
cutoffages <- c(6000, 25000)

# Choose a common prior on number of changepoints
prior_n_internal_changepoints_lambda <- 6


##############################################
### Read in the Guthrie data

# Read in 14C information on human occupation
data <- read.csv("data/PleistoceneDates/Humanc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
keep <- which(data$X14C < cutoffages[2] & data$X14C > cutoffages[1])
data <- data[keep,]
removena <- which(is.na(data$X14C) | is.na(data$X1..Sigma))
if(length(removena) != 0) { # Remove any values where 14C or sigma is NA
  data <- data[-removena,]
}
n <- nrow(data)
x <- data$X14C
xsig <- data$X1..Sigma
Human14C <- data.frame(C14age = x, C14sig = xsig, Species = "Human")

# Read in 14C information on bison
data <- read.csv("data/PleistoceneDates/Bisonc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
keep <- which(data$X14C < cutoffages[2] & data$X14C > cutoffages[1])
data <- data[keep,]
removena <- which(is.na(data$X14C) | is.na(data$X1..Sigma))
if(length(removena) != 0) { # Remove any values where 14C or sigma is NA
  data <- data[-removena,]
}
n <- nrow(data)
x <- data$X14C
xsig <- data$X1..Sigma
Bison14C <- data.frame(C14age = x, C14sig = xsig, Species = "Bison")

# Read in 14C information on mammoth
data <- data <- read.csv("data/PleistoceneDates/Mammothc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
keep <- which(data$X14C < cutoffages[2] & data$X14C > cutoffages[1])
data <- data[keep,]
removena <- which(is.na(data$X14C) | is.na(data$X1..Sigma))
if(length(removena) != 0) { # Remove any values where 14C or sigma is NA
  data <- data[-removena,]
}
n <- nrow(data)
x <- data$X14C
xsig <- data$X1..Sigma
Mammoth14C <- data.frame(C14age = x, C14sig = xsig, Species = "Mammoth")

# Read in 14C information on alces/moose
data <- read.csv("data/PleistoceneDates/Alcesc14Dates.csv", header = TRUE, na.strings = c("", "greater than"))
keep <- which(data$X14C < cutoffages[2] & data$X14C > cutoffages[1])
data <- data[keep,]
removena <- which(is.na(data$X14C) | is.na(data$X1..Sigma))
if(length(removena) != 0) { # Remove any values where 14C or sigma is NA
  data <- data[-removena,]
}
n <- nrow(data)
x <- data$X14C
xsig <- data$X1..Sigma
Alces14C <- data.frame(C14age = x, C14sig = xsig, Species = "Alces")

# Combine all data
Pleist14C <- rbind(Human14C, Bison14C, Mammoth14C, Alces14C)
Pleist14C$Species <- factor(Pleist14C$Species,
                            levels = c("Human", "Alces", "Bison", "Mammoth"))



##################################################
# Fit a PP model to each of the datasets and create later plots
##################################################

# Find the plausible calendar age range
min_calendar_age_range <- intcal20$calendar_age_BP[which.min(abs(intcal20$c14_age - cutoffages[1]))] - 250
max_calendar_age_range <- intcal20$calendar_age_BP[which.min(abs(intcal20$c14_age - cutoffages[2]))] + 250
calendar_age_range <- c(min_calendar_age_range, max_calendar_age_range)

# Read in cold/warm periods
Heinrich <- read.csv("data/PleistoceneDates/HeinrichDates.csv", header = FALSE)
names(Heinrich) <- c("Event", "calage")
YD <- data.frame(Event = c("end YD", "start YD"), calage = c(11703, 12896) - 50)
Heinrich <- rbind(YD, Heinrich)


# Run code to fit PP model and create plotted output

### Each example will create three plots:
# Plot i) The posterior mean of the sample occurrence rate
# Plot ii) The posterior estimates of the locations of the changepoints (conditional on number)
# Plot iii) A histogram of the posterior number of changepoints
### These are merged into Figs 9 and 10

set.seed(17)
source("R/FitPPHumans.R")
source("R/FitPPAlces.R")
source("R/FitPPBison.R")
source("R/FitPPMammoth.R")


# Return to user specified par parameters (before running code)
par(oldpar)







