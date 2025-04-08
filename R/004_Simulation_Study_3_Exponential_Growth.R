# This file performs Simulation Study 3
# Creates a Poisson process with an exponentially increasing occurrence rate and then fits a PP model to it
# This simulation is analogous to Example 1 from Crema and Shoda (PloS, 2021)

# Creates three plots used for Fig 8:
# Plot i) The posterior mean of the sample occurrence rate
# Plot ii) A histogram of the posterior number of changepoints
# Plot iii) The posterior estimates of the locations of the changepoints (conditional on number)


# Load carbondate library
library(carbondate)

# Set the mean number of internal changepoints in the prior
# (the default in the library function will select 3)
prior_n_internal_changepoints_lambda <- 3


# Plotting parameters
plot_scale_fac <- 1.4
plot_height <- 9
plot_width <- 8
pdf_output <- FALSE
abline_col <- rgb(219,62,177, alpha = 255*0.25, maxColorValue = 255)
abline_lty <- 1
abline_lwd <- 2
SPD_colour <- grey(0.1, alpha = 0.1)
add_true_rate <- TRUE
add_SPD <- TRUE
lab_adj <- 0.02
panel_label_cex <- plot_scale_fac * 4/3

# Decide if write plots to file
write_plots_to_file <- FALSE

# Store user par parameters (and revert back to these at the end)
oldpar <- par(no.readonly = TRUE)

# Fix seed for reproducibility
set.seed(16)


# Select exponential growth model to match Crema and Shoda (2021)

# Start/end dates and growth rate
min_cal <- 4000
max_cal <- 6000
growth_rate <- 0.003
n_expected_samples <- 500

# Create calendar age grid on which to sample events
t <- seq(min_cal, max_cal, by = 1)
t <- t[-length(t)]

# Set upper/lower calendar range for Poisson process
calendar_age_range_BP <- range(t) + c(-200, 200)

# Create exponential lambda with chosen growth rate
lambda <- exp(growth_rate * (max_cal - t))
# Adjust so that expect n_expected_samples in the time period
lambda <- (n_expected_samples/sum(lambda)) * lambda
true_rate <- lambda

# Sample calendar ages (according to Poisson process)
n_samp <- rpois(1, sum(true_rate))
theta_true <- sample(t, n_samp, replace = TRUE, prob = true_rate)

# Sample some corresponding 14C determinations
rc_sigmas <- rep(25, n_samp)
rc_determinations <- rnorm(n_samp,
                           mean = approx(intcal20$calendar_age_BP,
                                         intcal20$c14_age,
                                         theta_true)$y,
                           sd = sqrt(rc_sigmas^2 +
                                       (approx(intcal20$calendar_age_BP,
                                               intcal20$c14_sig,
                                               theta_true)$y)^2) )

# Find normalising constant so that SPD can be compared to true_rate
norm_constant_true_rate <- sum(true_rate)

############################################
### Fit a PP to the data

# Fit Poisson process to this data and plot
PP_fit_output_simulation_3 <- PPcalibrate(
  rc_determinations = rc_determinations,
  rc_sigmas = rc_sigmas,
  prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
  calibration_curve = intcal20,
  calendar_age_range = calendar_age_range_BP,
  calendar_grid_resolution = 1,
  n_iter = 1e5,
  n_thin = 10,
  show_progress = TRUE)

#####################################################
######### Create the plots of the posterior

## Plot i) The posterior mean occurrence rate

# Decide if write plots to a file
if(write_plots_to_file) {
  if(pdf_output) {
    pdf(paste("output/SimulationStudy3/PosteriorMean_Prior_",
              prior_n_internal_changepoints_lambda,
              "_Internal_Changes.pdf", sep = ""),
        width = 3 * plot_width/plot_scale_fac,
        height = plot_height/plot_scale_fac)
  } else {
    png(paste("output/SimulationStudy3/PosteriorMean_Prior_",
              prior_n_internal_changepoints_lambda,
              "_Internal_Changes.png", sep = ""),
        width = 3 * plot_width/plot_scale_fac,
        height = plot_height/plot_scale_fac,
        units = "in", res = 480)
  }
}

simulated_PP_mean_fit <- PlotPosteriorMeanRate(PP_fit_output_simulation_3,
                                               show_individual_means = FALSE)
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


# Add SPD (normalised so comparable to truerate)
if(add_SPD) {
  SPD_initial_fit <- FindSummedProbabilityDistribution(
    calendar_age_range_BP= range(simulated_PP_mean_fit$calendar_age_BP),
    rc_determinations= rc_determinations,
    rc_sigmas = rc_sigmas,
    calibration_curve=intcal20,
    plot_output = FALSE, plot_pretty = FALSE)
  polygon(
    c(SPD_initial_fit$calendar_age_BP, rev(SPD_initial_fit$calendar_age_BP)),
    norm_constant_true_rate * c(SPD_initial_fit$probability, rep(0, length(SPD_initial_fit$probability))),
    border = NA,
    col = SPD_colour)
}

if(add_true_rate) {
  lines(c(min(t), t, max(t)), c(0,true_rate,0), lty = 1, col = "red")
}

mtext(LETTERS[1], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)

if(write_plots_to_file) {
  dev.off()
}


###########################################################################
### Plot ii) The posterior number of internal changes

# Decide if write plots to a file
if(write_plots_to_file) {
  if(pdf_output) {
    pdf(paste("output/SimulationStudy3/Number_Changepoints_Prior_",
              prior_n_internal_changepoints_lambda,
              "_Internal_Changes.pdf", sep = ""),
        width = 1.5 * plot_width/plot_scale_fac,
        height = plot_height/plot_scale_fac)
  } else {
    png(paste("output/SimulationStudy3/Number_Changepoints_Prior_",
              prior_n_internal_changepoints_lambda,
              "_Internal_Changes.png", sep = ""),
        width = 1.5 * plot_width/plot_scale_fac,
        height = plot_height/plot_scale_fac,
        units = "in", res = 480)
  }
}


PlotNumberOfInternalChanges(PP_fit_output_simulation_3)

mtext(LETTERS[2], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)

if(write_plots_to_file) {
  dev.off()
}



#######################################################################
## Plot iii) The locations of the changes

# Decide if write plots to a file
if(write_plots_to_file) {
  if(pdf_output) {
    pdf(paste("output/SimulationStudy3/Changepoint_Locations_Prior_",
              prior_n_internal_changepoints_lambda,
              "_Internal_Changes.pdf", sep = ""),
        width = 1.5 * plot_width/plot_scale_fac,
        height = plot_height/plot_scale_fac)
  } else {
    png(paste("output/SimulationStudy3/Changepoint_Locations_Prior_",
              prior_n_internal_changepoints_lambda,
              "_Internal_Changes.png", sep = ""),
        width = 1.5 * plot_width/plot_scale_fac,
        height = plot_height/plot_scale_fac,
        units = "in", res = 480)
  }
}

# Make sure that the x-axis is the same as for the plot of posterior mean
# Requires xaxs and yas = "i"
par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)

PlotPosteriorChangePoints(PP_fit_output_simulation_3, n_changes = c(5, 6))
# Overlay the endpoints of the exponential
abline(v = min_cal, col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = max_cal, col = abline_col, lty = abline_lty, lwd = abline_lwd)


mtext(LETTERS[3], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)

if(write_plots_to_file) {
  dev.off()
}


# Return to user specified par parameters (before running code)
par(oldpar)

