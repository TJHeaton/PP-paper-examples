# This code implements simulation study 1
# An example where there is a single uniform phase model (i.e. two change points)
library(carbondate)

# Select plotting parameters
plot_scale_fac <- 1.4
plot_height <- 9
plot_width <- 8
SPD_colour <- grey(0.1, alpha = 0.1)
add_true_rate <- TRUE
add_SPD <- TRUE
lab_adj <- 0.02
panel_label_cex <- plot_scale_fac * 4/3

# Decide if want plotting output as pdf (TRUE) or png (FALSE)
pdf_output <- FALSE

# Create some simulated calendar dates from min_cal to max_cal
set.seed(14)
min_cal <- 2050
max_cal <- 2100

n_samp <- 50
theta_true <- sample(min_cal:max_cal, size = n_samp, replace = TRUE)

# Fix the maximum calendar age range for analysis
calendar_age_range_BP <- c(1850, 2350)

#####################################################
# Sample/simulate some 14C determinations for these underlying calendar ages

# Choose +/- obs_radiocarbon_sigma as measurement uncertianty
obs_radiocarbon_sigmas <- rep(25, n_samp)

# Create simulated 14C determinations
obs_radiocarbon_ages <- rnorm(n_samp,
                          mean = approx(intcal20$calendar_age_BP,
                                        intcal20$c14_age,
                                        theta_true)$y,
                          sd = sqrt(obs_radiocarbon_sigmas^2 +
                                      (approx(intcal20$calendar_age_BP,
                                              intcal20$c14_sig,
                                              theta_true)$y)^2) )



# Fit Poisson process to this data and plot
PP_fit_output_simulation_1 <- PPcalibrate(
  rc_determinations = obs_radiocarbon_ages,
  rc_sigmas = obs_radiocarbon_sigmas,
  calibration_curve = intcal20,
  calendar_age_range = calendar_age_range_BP,
  calendar_grid_resolution = 1,
  n_iter = 1e5,
  n_thin = 10,
  show_progress = TRUE)


######### Now create plots to show results

# Now plot the Poisson process rate
if(pdf_output) {
  pdf("Plots/SimulationStudy1/PosteriorMean.pdf",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac)
} else {
  png("Plots/SimulationStudy1/PosteriorMean.png",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac,
      units = "in", res = 480)
}

## Plot 1: The posterior mean occurrence rate
simulated_PP_mean_fit <- PlotPosteriorMeanRate(PP_fit_output_simulation_1)
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- rev(range(simulated_PP_mean_fit$calendar_age_BP))
ylim <- c(0, 3 * max(simulated_PP_mean_fit$rate_mean))
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

# Add SPD
if(add_SPD) {
  SPD_fit <- FindSummedProbabilityDistribution(
    calendar_age_range_BP= range(simulated_PP_mean_fit$calendar_age_BP),
    rc_determinations= obs_radiocarbon_ages,
    rc_sigmas = obs_radiocarbon_sigmas,
    calibration_curve=intcal20,
    plot_output = FALSE, plot_pretty = FALSE)
  polygon(
    c(SPD_fit$calendar_age_BP, rev(SPD_fit$calendar_age_BP)),
    50 * c(SPD_fit$probability, rep(0, length(SPD_fit$probability))),
    border = NA,
    col = SPD_colour)
}

# Overlay the true underlying changepoints
abline(v = min_cal, col = rgb(219,62,177, alpha = 255*0.25, maxColorValue = 255), lty = 1, lwd = 4)
abline(v = max_cal, col = rgb(219,62,177, alpha = 255*0.25, maxColorValue = 255), lty = 1, lwd = 4)

mtext(LETTERS[1], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)

dev.off()



#######################################################################
### Plot 2: The posterior number of internal changes
if(pdf_output) {
  pdf("Plots/SimulationStudy1/Number_Changepoints.pdf",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac)
} else {
  png("Plots/SimulationStudy1/Number_Changepoints.png",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac,
      units = "in", res = 480)
}

PlotNumberOfInternalChanges(PP_fit_output_simulation_1)

mtext(LETTERS[2], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)
dev.off()




#######################################################################
## Plot 3: The locations of the changes
if(pdf_output) {
  pdf("Plots/SimulationStudy1/Changepoint_Locations.pdf",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac)
} else {
  png("Plots/SimulationStudy1/Changepoint_Locations.png",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac,
      units = "in", res = 480)
}


par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)

PlotPosteriorChangePoints(PP_fit_output_simulation_1)
xlim <- rev(range(simulated_PP_mean_fit$calendar_age_BP))
ylim <- c(0, 3 * max(simulated_PP_mean_fit$rate_mean))
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
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
abline(v = min_cal, col = rgb(219,62,177, alpha = 255*0.25, maxColorValue = 255), lty = 1, lwd = 4)
abline(v = max_cal, col = rgb(219,62,177, alpha = 255*0.25, maxColorValue = 255), lty = 1, lwd = 4)

# Add panel label
mtext(LETTERS[3], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)

dev.off()

