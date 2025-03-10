# Fit PP process to HUMANS and plot using Guthrie (2006)
# NOTE: You must have sourced 001_Analyse_Guthrie_Pleistocene_c14_Dates.R beforehand to get data


# Fit PP model to the human data
human_PP_fit_output <- PPcalibrate(
  rc_determinations = Human14C$C14age,
  rc_sigmas = Human14C$C14sig,
  calibration_curve = intcal20,
  calendar_age_range = calendar_age_range,
  calendar_grid_resolution = 1,
  prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
  n_iter = 1e5,
  n_thin = 10,
  show_progress = TRUE)

## Plot the posterior mean occurrence rate
out_file_name <- paste("output/PriorMeanChangepoints",
                       prior_n_internal_changepoints_lambda,
                       "/FitPP_Human_PriorMean_",
                       prior_n_internal_changepoints_lambda,
                       "_Changes", sep = "")

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

# Plot the posterior mean for rate of human occupation
human_PP_mean_fit <- PlotPosteriorMeanRate(human_PP_fit_output,
                                           show_individual_means = FALSE)
# We do not show the individual posterior mean calendar ages as it is a bit confusing

# Create new overlay (with same calendar age scale) on which we plot the cold/warm periods
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- rev(range(human_PP_mean_fit$calendar_age_BP))
ylim <- c(0, 3 * max(human_PP_mean_fit$rate_mean))
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

text(x = 21000, y = 0.25 *max(ylim), labels = "Human", col = "magenta", cex = 2 )

xtick <- seq(6000, 30000, by= 500)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.025)
dev.off()


##############################################
#Plot changepoint locations
out_file_name <- paste("output/PriorMeanChangepoints",
                       prior_n_internal_changepoints_lambda,
                       "/FitPP_Human_PriorMean_",
                       prior_n_internal_changepoints_lambda,
                       "_Locations_Changepoint", sep = "")
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

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
PlotPosteriorChangePoints(human_PP_fit_output, n_changes = c(4,5, 6))
par(new = TRUE)
xlim <- rev(range(human_PP_mean_fit$calendar_age_BP))
ylim <- c(0, 3 * max(human_PP_mean_fit$rate_mean))
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

text(x = 21000, y = 0.25 *max(ylim), labels = "Human", col = "magenta", cex = 2 )

xtick<-seq(6000, 30000, by= 500)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.025)
dev.off()


#######################################################################
### Plot the posterior number of internal changes


out_file_name <- paste("output/PriorMeanChangepoints",
                       prior_n_internal_changepoints_lambda,
                       "/FitPP_Human_PriorMean_",
                       prior_n_internal_changepoints_lambda,
                       "_Number_Changepoint", sep = "")

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

par(cex.axis = 1.6,
    cex.lab = 1.6,
    mar = c(5, 4.5, 0.85, 0.55) + 0.1)


PlotNumberOfInternalChanges(human_PP_fit_output)

dev.off()


