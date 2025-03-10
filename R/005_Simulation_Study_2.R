# This file performs simulation study 2
# Creates the rate illustration plot (with four internal changes)
# and then fits a PP model to it

library(carbondate)

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




# Create a simulate lambda (piecewise constant) and plot it with the data
set.seed(16)

# Artificial start and end dates
min_cal <- 1950
max_cal <- 3100
t <- seq(min_cal, max_cal, by = 1)
t <- t[-length(t)]
calendar_age_range_BP <- range(t) + c(-200, 200)


# Create a simulated lambda (piecewise constant)
lambda <- c(rep(0.15, 350), rep(0.7, 400), rep(0.2, length(t) - (400 + 350)))

# Choose beta(mutiplier)
beta <- 0.4
true_rate <- beta * lambda

# Sample calendar ages (according to Poisson process)
n_samp <- rpois(1, sum(true_rate))
theta_true <- sample(t, n_samp, replace = TRUE, prob = true_rate)



# Sample some 14C determinations
rc_sigmas <- rep(25, n_samp)
rc_determinations <- rnorm(n_samp,
                           mean = approx(intcal20$calendar_age_BP,
                                         intcal20$c14_age,
                                         theta_true)$y,
                           sd = sqrt(rc_sigmas^2 +
                                       (approx(intcal20$calendar_age_BP,
                                               intcal20$c14_sig,
                                               theta_true)$y)^2) )


# Normalise so SPD can be compared to true_rate
norm_constant_true_rate <- sum(true_rate)


#########################################################
### Create plot to show the underlying (true) PP rate
if(pdf_output) {
  pdf("Plots/IllustratePP.pdf",
      width = 8,
      height = 3)
} else {
  png("Plots/IllustratePP.png",
      width = 8,
      height = 3,
      units = "in", res = 480)
}

# Layout as two plots
par(
  mfrow = c(1,2),
  mgp = c(3, 0.7, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(3, 4.5, 2, 2) + 0.1,
  las = 1)

# Create panel A showing rate and calendar ages
plot(c(min(t), t, max(t)),
     c(0,lambda,0),
     lty = 1, lwd = 2,
     col = "red", type = "l",
     ylim = c(0, max(lambda) + 0.1),
     xlim = c(max_cal, min_cal) + c(100, -100), yaxs = "i",
     xlab = "",
     ylab = expression(paste("Occurrence Rate ", lambda, "(", theta, ")" )),
     main = "Occurence of Samples")

# Add ylabel and smaller ticks
mtext("Calendar Age (cal yr BP)",
      side = 1,
      line = 2)
xtick<-seq(1500, 4000, by = 20)
ytick<-seq(0, 1, by = 0.04)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)
rug(theta_true, side = 1,
    col = "purple",
    lwd = 2, ticksize = 0.05,
    quiet = TRUE)
abline(v = t[350], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[350 + 400], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[1], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[length(t)], col = abline_col, lty = abline_lty, lwd = abline_lwd)

# Add panel label
mtext(LETTERS[1], side = 3, adj = lab_adj, cex = 0.8 * panel_label_cex/plot_scale_fac, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)



# Create panel B
# Plot 14C determinations against IntCal20
plot(intcal20$calendar_age_BP,
     intcal20$c14_age,
     col = "blue",
     ylim = range(rc_determinations) + c(-4,4) * max(rc_sigmas),
     xlim = c(max_cal, min_cal) + c(100, -100),
     xlab = "",
     ylab = expression(paste("Radiocarbon age (", ""^14, "C ", "yr BP)")),
     type = "l", main = expression(paste(""^14,"C Changepoint Analysis")))

# Add label and smaller ticks
mtext("Calendar Age (cal yr BP)",
      side = 1,
      line = 2)
legend("topright",
       lty = 1, col = "red",
       expression(paste(lambda, "(", theta, ")" )))
xtick<-seq(1500, 4000, by=20)
ytick<-seq(1500, 4000, by=20)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)

lines(intcal20$calendar_age_BP,
      intcal20$c14_age + 2 * intcal20$c14_sig,
      lty = 2,
      col = "blue" )
lines(intcal20$calendar_age_BP,
      intcal20$c14_age - 2 * intcal20$c14_sig,
      lty = 2,
      col = "blue" )
polygon(c(rev(intcal20$calendar_age_BP), intcal20$calendar_age_BP),
        c(rev(intcal20$c14_age - 2 * intcal20$c14_sig),
          intcal20$c14_age + 2 * intcal20$c14_sig),
        col=rgb(0,0,1,.3), border=NA)
rug(rc_determinations, side = 2, quiet = TRUE)

# Plot the true rate along the bottom
par(new = TRUE)
plot(c(min(t), t, max(t)), c(0,true_rate,0), lty = 1, col = "red", type = "l",
     ylim = c(0, 10 * max(true_rate)), xlim = c(max_cal, min_cal) + c(100, -100),
     axes = FALSE, xlab = NA, ylab = NA, yaxs = "i")

mtext(LETTERS[2], side = 3, adj = lab_adj, cex = 0.8 * panel_label_cex/plot_scale_fac, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)

dev.off()




############################################
### Now fit a PP to the data

# Fit Poisson process to this data and plot
PP_fit_output_simulation_2 <- PPcalibrate(
  rc_determinations = rc_determinations,
  rc_sigmas = rc_sigmas,
  calibration_curve = intcal20,
  calendar_age_range = calendar_age_range_BP,
  calendar_grid_resolution = 1,
  n_iter = 1e5,
  n_thin = 10,
  show_progress = TRUE)



#####################################################
######### Now create the plots of the posterior

# Now plot the Poisson process rate
if(pdf_output) {
  pdf("Plots/SimulationStudy2/PosteriorMean.pdf",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac)
} else {
  png("Plots/SimulationStudy2/PosteriorMean.png",
      width = 3 * plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac,
      units = "in", res = 480)
}

## Plot 1: The posterior mean occurrence rate
simulated_PP_mean_fit <- PlotPosteriorMeanRate(PP_fit_output_simulation_2)
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
# Overlay the true underlying changepoints
abline(v = t[350], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[350 + 400], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[1], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[length(t)], col = abline_col, lty = abline_lty, lwd = abline_lwd)

# Add SPD (nomrlaised so comparable to truerate )
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

dev.off()



###########################################################################
### Plot 2: The posterior number of internal changes
if(pdf_output) {
  pdf("Plots/SimulationStudy2/Number_Changepoints.pdf",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac)
} else {
  png("Plots/SimulationStudy2/Number_Changepoints.png",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac,
      units = "in", res = 480)
}

PlotNumberOfInternalChanges(PP_fit_output_simulation_2)

mtext(LETTERS[2], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)
dev.off()




#######################################################################
## Plot 3: The locations of the changes
if(pdf_output) {
  pdf("Plots/SimulationStudy2/Changepoint_Locations.pdf",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac)
} else {
  png("Plots/SimulationStudy2/Changepoint_Locations.png",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac,
      units = "in", res = 480)
}



PlotPosteriorChangePoints(PP_fit_output_simulation_2, n_changes = c(3, 4, 5))
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
abline(v = t[350], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[350 + 400], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[1], col = abline_col, lty = abline_lty, lwd = abline_lwd)
abline(v = t[length(t)], col = abline_col, lty = abline_lty, lwd = abline_lwd)

mtext(LETTERS[3], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)

dev.off()



############################################################
### Plot 4: The posterior rates in the different periods
if(pdf_output) {
  pdf("Plots/SimulationStudy2/PosteriorHeights.pdf",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac)
} else {
  png("Plots/SimulationStudy2/PosteriorHeights.png",
      width = plot_width/plot_scale_fac,
      height = plot_height/plot_scale_fac,
      units = "in", res = 480)
}

PlotPosteriorHeights(PP_fit_output_simulation_2, n_changes = c(3, 4, 5))

# Overlay the true underlying rates
abline(v = unique(c(0,true_rate)), col = abline_col, lty = abline_lty, lwd = abline_lwd)

mtext(LETTERS[4], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 5/3 * (-1), outer = FALSE)

dev.off()










