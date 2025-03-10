# Illustrate the failings of SPDs
library(carbondate)

# Plotting parameters
plot_height <- 4
top_margin <- 2.5
lab_adj <- 0.02
panel_label_cex <- 4/3

# Decide if want plotting output as pdf (TRUE) or png (FALSE)
pdf_output <- FALSE


##############################################################
#### Plots/Illustration 1
#### Set of two plots
#### Plot A - Show SPDs overly variable and interpretation hard
#### Plot B - Polya Urn (not used in paper but just to show it works)
##############################################################

# Plot 1
out_file_name <- "Plots/SPDOverlyVariable"

if(pdf_output) {
  pdf(paste(out_file_name, ".pdf", sep = ""),
      width = 12,
      height = plot_height)
} else {
  png(paste(out_file_name, ".png", sep = ""),
      width = 12,
      height = plot_height,
      units = "in", res = 480)
}

oldpar <- par(no.readonly = TRUE)

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, top_margin, 2) + 0.1, las = 1)

# Fit SPD to inbuilt two normals dataset
two_normals_spd <- FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(2500, 6000),
  rc_determinations= two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20,
  plot_output = TRUE, plot_pretty = FALSE)

par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, top_margin, 2) + 0.1,
    las = 1)
xlim <- rev(range(two_normals_spd$calendar_age_BP))
ylim <- c(0, 3 * max(two_normals_spd$probability))
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
# Show true underlying calendar age density
weights_true <- c(0.45, 0.55)
cluster_means_true_calBP <- c(3500, 5000)
cluster_precisions_true <- 1 / c(200, 100)^2

# Find and plot true exact density
truedens <- function(t, w, truemean, trueprec) {
  dens <- 0
  for(i in 1:length(w)) {
    dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
  }
  dens
}
curve(truedens(
  x,
  w = weights_true,
  truemean = cluster_means_true_calBP,
  trueprec = cluster_precisions_true),
  from = 2500, to = 7000, n = 401,
  lwd = 2,
  col = "red", add = TRUE)


# Reset plotting parameters
par(oldpar)

dev.off()


#####
# Plot B - show how a DPMM works to recreate correct calendar age density
out_file_name <- "Plots/DPMMFit"

if(pdf_output) {
  pdf(paste(out_file_name, ".pdf", sep = ""),
      width = 12,
      height = plot_height)
} else {
  png(paste(out_file_name, ".png", sep = ""),
      width = 12,
      height = plot_height,
      units = "in", res = 480)
}

oldpar <- par(no.readonly = TRUE)

polya_urn_output <- PolyaUrnBivarDirichlet(
  rc_determinations = two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20,
  n_iter = 1e4,
  show_progress = FALSE)

two_normals_DPMM <- PlotPredictiveCalendarAgeDensity(
  output_data = polya_urn_output,
  show_SPD = TRUE)

par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- rev(range(two_normals_DPMM[[1]]$calendar_age_BP))
ylim <- c(0, 3 * max(two_normals_DPMM[[1]]$density_mean))
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
# Show true underlying calendar age density
weights_true <- c(0.45, 0.55)
cluster_means_true_calBP <- c(3500, 5000)
cluster_precisions_true <- 1 / c(200, 100)^2

# Find and plot true exact density
truedens <- function(t, w, truemean, trueprec) {
  dens <- 0
  for(i in 1:length(w)) {
    dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
  }
  dens
}
curve(truedens(
  x,
  w = weights_true,
  truemean = cluster_means_true_calBP,
  trueprec = cluster_precisions_true),
  from = 2500, to = 7000, n = 401,
  lwd = 2,
  col = "red", add = TRUE)

# Reset plotting parameters
par(oldpar)

dev.off()




##############################################################
#### Plots/Illustration 2
#### Single plot to show single determination SPD doesn't make sense
#### there needs to be some modelling in the calendar age domain
##############################################################


out_file_name <- "Plots/SPDSingleDetermination"

if(pdf_output) {
  pdf(paste(out_file_name, ".pdf", sep = ""),
      width = 12,
      height = plot_height)
} else {
  png(paste(out_file_name, ".png", sep = ""),
      width = 12,
      height = plot_height,
      units = "in", res = 480)
}

# Illustrate how single determination calibrates to multiple peaks
## Single determination around 2100 cal BP
test_radiocarbon_date <- approx(intcal20$calendar_age_BP, intcal20$c14_age,
                                xout = 2010)$y + 40
test_radiocarbon_sigma <- 30

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, top_margin, 2) + 0.1, las = 1)


calibration_result <- CalibrateSingleDetermination(
  rc_determination = test_radiocarbon_date,
  rc_sigma = test_radiocarbon_sigma,
  F14C_inputs = FALSE,
  calibration_curve = intcal20,
  plot_output = TRUE, plot_pretty = FALSE)

dev.off()



##############################################################
#### Plots/Illustration 3
#### Set of four plots to show failure of confidence intervals
##############################################################


# Sample 50 from a uniform distribution

####### Create simulated truth
# Create some simulated calendar dates from mincal to maxcal
set.seed(14)
mincal <- 2050
maxcal <- 2100

nsamp <- 50
theta_true <- sample(mincal:maxcal, size = nsamp, replace = TRUE)

calendar_age_range_BP=c(1900, 2200) # It will fit 400 cal yrs ither side

# Sample some 14C determinations
radiocarbon_sigmas_SPD_confint<- rep(25, nsamp)
radiocarbon_ages_SPD_confint <- rnorm(nsamp,
                                      mean = approx(intcal20$calendar_age_BP,
                                                    intcal20$c14_age,
                                                    theta_true)$y,
                                      sd = sqrt(radiocarbon_sigmas_SPD_confint^2 +
                                                  (approx(intcal20$calendar_age_BP,
                                                          intcal20$c14_sig,
                                                          theta_true)$y)^2) )



## Plot A - Sampling from a uniform phase
# pdf("Plots/BootstrapFailure1.pdf", width = 8, height = plot_height)
png("Plots/BootstrapFailure1.png", width = 8, height = plot_height, units = "in", res = 480)
par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)

x_lim_plot <- c(2600,1500)

plot(intcal20$calendar_age_BP, intcal20$c14_age, col = "blue",
     ylab = expression(paste("Radiocarbon age (", ""^14, "C yr BP)")),
     xlab = "Calendar Age (cal yr BP)",
     xlim = x_lim_plot, ylim = c(1400, 2700),
     type = "l", main = "")

# Add ylabel and smaller ticks
#mtext("Calendar Age (cal yr BP)", side = 1, line = 2)
xtick<-seq(1400, 4000, by=20)
ytick<-seq(1400, 4000, by=20)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)


# Plot calibration curve
intcal20$ub <- intcal20$c14_age + 1.96 * intcal20$c14_sig
intcal20$lb <- intcal20$c14_age - 1.96 * intcal20$c14_sig
lines(intcal20$calendar_age_BP, intcal20$ub, lty = 3, col = "blue" )
lines(intcal20$calendar_age_BP, intcal20$lb, lty = 3, col = "blue")
polygon(c(rev(intcal20$calendar_age_BP), intcal20$calendar_age_BP),
        c(rev(intcal20$lb), intcal20$ub), col=rgb(0,0,1,.3), border=NA)
rug(radiocarbon_ages_SPD_confint, side = 2, ticksize = 0.03, lwd = 1, col = "red")

legend_labels <- c(
  substitute(paste(""^14, "C determination ")),
  "IntCal20",
  expression(paste("2", sigma, " interval")))
lty <- c(1, 1, 2)
lwd <- c(1, 1, 1)
pch <- c(NA, NA, NA)
col <- c(grDevices::rgb(1, 0, 0, .5), "blue", "blue")
legend("topright", legend = legend_labels, lty = lty, lwd=lwd, pch = pch, col = col)

# Create polygon of true density
t <- seq(mincal, maxcal, by = 1)
truerate <- rep(1, length(t))/length(t) # So integrates to 1

# Plot the true rate along the bottom
par(new = TRUE)
plot(c(min(t), t, max(t)), c(0,truerate,0), lty = 1, col = "red", type = "l",
     ylim = c(0, 1.2 * max(truerate)), xlim = x_lim_plot,
     axes = FALSE, xlab = NA, ylab = NA, yaxs = "i")
dens_polygon <- cbind(c(min(t), t, max(t)),
                      c(0,truerate,0))

polygon(dens_polygon, col = grDevices::rgb(1, 0, 0, .5))


mtext(LETTERS[1], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)
dev.off()



# Plot B - showing the SPD of this density (doesn't look like the original)
# pdf("Plots/BootstrapFailure2.pdf", width = 8, height = plot_height)
png("Plots/BootstrapFailure2.png", width = 8, height = plot_height, units = "in", res = 480)

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)

SPD_initial_fit <- FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(1900, 2200),
  rc_determinations= radiocarbon_ages_SPD_confint,
  rc_sigmas = radiocarbon_sigmas_SPD_confint,
  calibration_curve=intcal20,
  plot_output = TRUE, plot_pretty = FALSE)


par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)
xlim <- rev(range(SPD_initial_fit$calendar_age_BP))
ylim <- c(0, 3 * max(SPD_initial_fit$probability))
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
polygon(dens_polygon,
        border = NA,
        col = grDevices::rgb(1, 0, 0, .25))
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)


# Now resample 50 observations from SPD
set.seed(32)
resampled_calendar_ages <- sample(
  x = SPD_initial_fit$calendar_age_BP,
  size = nsamp,
  replace = TRUE,
  prob = SPD_initial_fit$probability)

resampled_radiocarbon_ages <- rnorm(nsamp,
                                    mean = approx(intcal20$calendar_age_BP,
                                                  intcal20$c14_age,
                                                  resampled_calendar_ages)$y,
                                    sd = sqrt(radiocarbon_sigmas_SPD_confint^2 +
                                                (approx(intcal20$calendar_age_BP,
                                                        intcal20$c14_sig,
                                                        resampled_calendar_ages)$y)^2) )
rug(resampled_calendar_ages, col = "purple")

mtext(LETTERS[2], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)

dev.off()

# Plot C - SPD on resampled/bootstrap sample
#pdf("Plots/BootstrapFailure3.pdf", width = 8, height = plot_height)
png("Plots/BootstrapFailure3.png", width = 8, height = plot_height, units = "in", res = 480)

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)

SPD_bootstrap_fit <- FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(1900, 2200),
  rc_determinations= resampled_radiocarbon_ages,
  rc_sigmas = radiocarbon_sigmas_SPD_confint,
  calibration_curve=intcal20,
  plot_output = TRUE, plot_pretty = FALSE)

par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)
xlim <- rev(range(SPD_bootstrap_fit$calendar_age_BP))
ylim <- c(0, 3 * max(SPD_bootstrap_fit$probability))
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
polygon(dens_polygon,
        border = NA,
        col = grDevices::rgb(1, 0, 0, .25))

bootstrap_SPD_fit_polygon <- cbind(c(min(SPD_bootstrap_fit$calendar_age_BP),
                                     SPD_bootstrap_fit$calendar_age_BP,
                                     max(SPD_bootstrap_fit$calendar_age_BP)),
                                   c(0,
                                     SPD_bootstrap_fit$probability,
                                     0)
)

polygon(bootstrap_SPD_fit_polygon,
        border = NA,
        col = grDevices::rgb(1, 0, 1, .25))

axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)

mtext(LETTERS[3], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)

dev.off()






#########################################################################
## Create code which calculates a bootstrap-based confidence interval
# set.seed(14)
# mincal <- 2050
# maxcal <- 2100
#
# nsamp <- 50
# theta_true <- sample(mincal:maxcal, size = nsamp, replace = TRUE)
#
#
# # Create polygon of true density
# t <- seq(mincal, maxcal, by = 1)
# truerate <- rep(1, length(t))/length(t) # So integrates to 1
#
# # Plot the true rate along the bottom
# dens_polygon <- cbind(c(min(t), t, max(t)),
#                       c(0,truerate,0))
#
#
#
# # Sample some 14C determinations
# radiocarbon_sigmas_SPD_confint<- rep(25, nsamp)
# radiocarbon_ages_SPD_confint <- rnorm(nsamp,
#                                       mean = approx(intcal20$calendar_age_BP,
#                                                     intcal20$c14_age,
#                                                     theta_true)$y,
#                                       sd = sqrt(radiocarbon_sigmas_SPD_confint^2 +
#                                                   (approx(intcal20$calendar_age_BP,
#                                                           intcal20$c14_sig,
#                                                           theta_true)$y)^2) )
#
#
# calendar_age_range_BP <- c(1900, 2200)


# SPD_initial_fit <- FindSummedProbabilityDistribution(
#   calendar_age_range_BP = calendar_age_range_BP,
#   rc_determinations = radiocarbon_ages_SPD_confint,
#   rc_sigmas = radiocarbon_sigmas_SPD_confint,
#   calibration_curve=intcal20,
#   plot_output = TRUE, plot_pretty = FALSE)


FindBootstrapSPDInterval <- function(SPD_initial_fit,
                                     radiocarbon_obs_sigmas,
                                     calendar_age_range_BP,
                                     n_boot = 100,
                                     CI_probs = c(0.025, 0.975)) {

  n_samp <- length(radiocarbon_obs_sigmas)
  n_resamples_needed <- n_samp * n_boot

  # Create matrix of resamples (each row corresponds to a bootstrap resample i)
  resampled_calendar_ages <- sample(
    x = SPD_initial_fit$calendar_age_BP,
    size = n_resamples_needed,
    replace = TRUE,
    prob = SPD_initial_fit$probability)

  # Mean of determinations resampled from SPD
  resampled_radiocarbon_obs_mean <- approx(intcal20$calendar_age_BP,
                                           intcal20$c14_age,
                                           resampled_calendar_ages)$y

  # Sd of determinations resampled from SPD
  # Combination of lab-reported/observed and calibration curve
  resampled_radiocarbon_obs_sigmas <- sqrt(rep(radiocarbon_obs_sigmas, n_boot)^2 +
                                             (approx(intcal20$calendar_age_BP,
                                                     intcal20$c14_sig,
                                                     resampled_calendar_ages)$y)^2)


  resampled_radiocarbon_ages <- rnorm(n_resamples_needed,
                                      mean = resampled_radiocarbon_obs_mean,
                                      sd = resampled_radiocarbon_obs_sigmas)

  # Convert to matrix
  bootstrap_resample_radiocarbon_ages <- matrix(resampled_calendar_ages, byrow = TRUE, nrow = n_boot)

  # Now create SPD for each row
  n_calendar_grid <- length(SPD_initial_fit$calendar_age_BP)
  bootstrap_SPD <- matrix(NA, ncol = n_calendar_grid, nrow = n_boot)

  for(i in 1:n_boot) {
  bootstrap_SPD[i,] <- FindSummedProbabilityDistribution(
    calendar_age_range_BP = calendar_age_range_BP,
    rc_determinations = bootstrap_resample_radiocarbon_ages[i,],
    rc_sigmas = radiocarbon_obs_sigmas,
    calibration_curve= intcal20,
    plot_output = FALSE, plot_pretty = FALSE)$probability
  }

  # Now find CI
  bootstrap_SPD_lower_CI <- apply(bootstrap_SPD, 2, quantile, probs = CI_probs[1])
  bootstrap_SPD_upper_CI <- apply(bootstrap_SPD, 2, quantile, probs = CI_probs[2])

  retlist <- list(calendar_age_BP = SPD_initial_fit$calendar_age_BP,
                  bootstrap_SPD_lower_CI = bootstrap_SPD_lower_CI,
                  bootstrap_SPD_upper_CI = bootstrap_SPD_upper_CI)
  return(retlist)
}

bootstrap_SPD_CI <- FindBootstrapSPDInterval(SPD_initial_fit = SPD_initial_fit,
                         radiocarbon_obs_sigmas = radiocarbon_sigmas_SPD_confint,
                         calendar_age_range_BP = calendar_age_range_BP,
                         n_boot = 5000)




#############################
# Create plot showing SPD bootstrap intervals

png("Plots/BootstrapFailure4.png", width = 8, height = plot_height, units = "in", res = 480)

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)

SPD_initial_fit <- FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(1900, 2200),
  rc_determinations= radiocarbon_ages_SPD_confint,
  rc_sigmas = radiocarbon_sigmas_SPD_confint,
  calibration_curve=intcal20,
  plot_output = TRUE, plot_pretty = FALSE)




par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)
xlim <- rev(range(SPD_initial_fit$calendar_age_BP))
ylim <- c(0, 3 * max(SPD_initial_fit$probability))
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
polygon(dens_polygon,
        border = NA,
        col = grDevices::rgb(1, 0, 0, .25))


# Plot the bootstrap SPD 95% confidence intervals
lines(bootstrap_SPD_CI$calendar_age_BP,
      bootstrap_SPD_CI$bootstrap_SPD_lower_CI,
      col = "black",
      lty = 2)

lines(bootstrap_SPD_CI$calendar_age_BP,
      bootstrap_SPD_CI$bootstrap_SPD_upper_CI,
      col = "black",
      lty = 2)

mtext(LETTERS[4], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)


dev.off()


