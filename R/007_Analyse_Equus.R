# Analysis of Equus data (using in-built library)
set.seed(19)

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


# Fit PP model to the equus data
equus_PP_fit_output <- PPcalibrate(
  rc_determinations = equus$c14_age,
  rc_sigmas = equus$c14_sig,
  calibration_curve = intcal20,
  calendar_age_range = calendar_age_range,
  calendar_grid_resolution = 1,
  prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
  n_iter = 1e5,
  n_thin = 10,
  show_progress = TRUE)

## Plot the posterior mean occurrence rate
out_file_name <- paste("output/PleistoceneMegafauna/FitPP_Equus_PriorMean_",
                       prior_n_internal_changepoints_lambda,
                       "_Changes", sep = "")

# Decide if write plots to a file
if(write_plots_to_file) {
  if(pdf_output) {
    pdf(paste(out_file_name, ".pdf", sep = ""),
        width = 16,
        height = 4)
  } else {
    png(paste(out_file_name, ".png", sep = ""),
        width = 16,
        height = 4,
        units = "in", res = 480)
  }
}

# Plot the posterior mean for rate of equus
equus_PP_mean_fit <- PlotPosteriorMeanRate(equus_PP_fit_output,
                                           show_individual_means = FALSE)
# We do not show the individual posterior mean calendar ages as it is a bit confusing

# Create new overlay (with same calendar age scale) on which we plot the cold/warm periods
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- rev(range(equus_PP_mean_fit$calendar_age_BP))
ylim <- c(0, 3 * max(equus_PP_mean_fit$rate_mean))
plot(
  NULL,
  NULL,
  type = "n",
  ylim = ylim,
  xlim = xlim,
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  xaxs = "i",
  yaxs = "i")
# Overlay the true underlying changepoints
for(i in seq(1, dim(Heinrich)[1], by = 2)) {
  tempcal <- Heinrich$calage[c(i, i+1)]
  tempx <- c(tempcal, rev(tempcal))
  tempy <- rep(par("usr")[3:4], c(2,2))
  polygon(tempx, tempy, col = rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50"))
  name <- as.character(Heinrich$Event[i])
  text(mean(tempcal), max(tempy), labels = substr(name, nchar(name)-1, nchar(name)), pos = 1)
}

text(x = 21000, y = 0.25 *max(ylim), labels = "Equus", col = "magenta", cex = 2 )

xtick <- seq(6000, 30000, by= 500)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.025)

if(write_plots_to_file) {
  dev.off()
}


##############################################
#Plot changepoint locations
out_file_name <- paste("output/PleistoceneMegafauna/FitPP_Equus_PriorMean_",
                       prior_n_internal_changepoints_lambda,
                       "_Locations_Changepoint", sep = "")

# Decide if write plots to a file
if(write_plots_to_file) {
  if(pdf_output) {
    pdf(paste(out_file_name, ".pdf", sep = ""),
        width = 16,
        height = 4)
  } else {
    png(paste(out_file_name, ".png", sep = ""),
        width = 16,
        height = 4,
        units = "in", res = 480)
  }
}

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
PlotPosteriorChangePoints(equus_PP_fit_output, n_changes = c(4,5, 6))
par(new = TRUE)
xlim <- rev(range(equus_PP_mean_fit$calendar_age_BP))
ylim <- c(0, 3 * max(equus_PP_mean_fit$rate_mean))
plot(
  NULL,
  NULL,
  type = "n",
  ylim = ylim,
  xlim = xlim,
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  xaxs = "i",
  yaxs = "i")

# Overlay the hotter/cooler time periods of Heinrich Events
for(i in seq(1, dim(Heinrich)[1], by = 2)) {
  tempcal <- Heinrich$calage[c(i, i+1)]
  tempx <- c(tempcal, rev(tempcal))
  tempy <- rep(par("usr")[3:4], c(2,2))
  polygon(tempx, tempy, col = rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50"))
  name <- as.character(Heinrich$Event[i])
  text(mean(tempcal), max(tempy), labels = substr(name, nchar(name)-1, nchar(name)), pos = 1)
}

text(x = 21000, y = 0.25 *max(ylim), labels = "Equus", col = "magenta", cex = 2 )

xtick<-seq(6000, 30000, by= 500)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.025)

if(write_plots_to_file) {
  dev.off()
}

#######################################################################
### Plot the posterior number of internal changes

out_file_name <- paste("output/PleistoceneMegafauna/FitPP_Equus_PriorMean_",
                       prior_n_internal_changepoints_lambda,
                       "_Number_Changepoint", sep = "")

# Decide if write plots to a file
if(write_plots_to_file) {
  if(pdf_output) {
    pdf(paste(out_file_name, ".pdf", sep = ""),
        width = 5,
        height = 5)
  } else {
    png(paste(out_file_name, ".png", sep = ""),
        width = 5,
        height = 5,
        units = "in", res = 480)
  }
  # Change plotting size for pdf/png files only
  par(cex.axis = 1.6,
      cex.lab = 1.6,
      mar = c(5, 4.5, 0.85, 0.55) + 0.1)
}




PlotNumberOfInternalChanges(equus_PP_fit_output)

if(write_plots_to_file) {
  dev.off()
}

