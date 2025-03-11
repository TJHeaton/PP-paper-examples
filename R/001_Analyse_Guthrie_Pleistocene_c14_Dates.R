# This file analyses the Dale Guthrie (2006) data of Pleistocene Megafauna:
# Specifically:
# 1 - Creates Figure 1 - A violin plot of Guthrie's 14C determinations
# 2 - Fits Poisson Process model to each species and creates plots/figures shown in real-life case study

# Install development version 1.0.1.9000 of carbondate library
devtools::install_github("TJHeaton/carbondate")
library(carbondate)

# Violin plot of the radiocarbon ages of the different species in Yukon
library(ggplot2)
library(patchwork)

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




##################################################################
### Create initial violin plots and plot human 14C determinations against IntCal20 curve

# Create violin plot
ggplot2::theme_update(plot.tag = element_text(face = "bold"))

p <- ggplot(Pleist14C, aes(x=Species, y=C14age, fill = Species))
p <- p + geom_violin(position="dodge", alpha=0.5)
# p <- p + coord_flip()
p <- p + geom_jitter(shape=16, position=position_jitter(0.02))
p <- p + labs(title = expression(paste("Density of ", ""^14, "C dates for Animals in Yukon and Alaska")), x="Species", y = expression(paste("Radiocarbon Age (", ""^14, "C yr", " BP)")))
p <- p + theme(aspect.ratio = 1)
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + theme(legend.position = "none")
p <- p + labs(tag = "A")
p

ggsave("output/Pleistocene14CDates.pdf", width = 5.98, height = 5.03, device = "pdf") # Height was 6.03

# Also create a png output
png("output/Pleistocene14CDates.png", width = 5.98, height = 5.03, units = "in", res = 480)
p
dev.off()


# Create panel B showing 14C determinations of humans against IntCal20
source("R/PlotHumansvsIntCal20.R")





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



# Run code to fit PP model and create plotted output for later figures
set.seed(17)
source("R/FitPPHumans.R")
source("R/FitPPAlces.R")
source("R/FitPPBison.R")
source("R/FitPPMammoth.R")


# Return to user specified par parameters (before running code)
par(oldpar)







