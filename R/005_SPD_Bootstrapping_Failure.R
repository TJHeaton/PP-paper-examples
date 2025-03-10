##############################################################
# File to illustrate bootstrap failure of SPDs
#### Set of four plots to show failure of bootstrap confidence intervals for SPDs
##############################################################

# Illustrate the failings of SPDs
library(carbondate)

# Plotting parameters
plot_height <- 4
top_margin <- 2.5
lab_adj <- 0.02
panel_label_cex <- 4/3
denscale <- 3

# Decide if want plotting output as pdf (TRUE) or png (FALSE)
pdf_output <- FALSE


####################################################################
####### Create simulated truth
# Create n_samp simulated calendar dates from uniform distribution U[min_cal, max_cal]
set.seed(14)
min_cal <- 2050
max_cal <- 2100

n_samp <- 50
theta_true <- sample(min_cal:max_cal, size = n_samp, replace = TRUE)

calendar_age_range_BP <- c(1900, 2200)
# It will fit an SPD considering 400 cal yrs either side (this is encoded into function)

# Sample some 14C determinations corresponding to these calendar ages
obs_radiocarbon_sigmas <- rep(25, n_samp)
obs_radiocarbon_ages <- rnorm(n_samp,
                                      mean = approx(intcal20$calendar_age_BP,
                                                    intcal20$c14_age,
                                                    theta_true)$y,
                                      sd = sqrt(obs_radiocarbon_sigmas^2 +
                                                  (approx(intcal20$calendar_age_BP,
                                                          intcal20$c14_sig,
                                                          theta_true)$y)^2) )


# Create initial SPD (helps find initial plotting parameters too so can make consistent plots)
SPD_initial_fit <- FindSummedProbabilityDistribution(
  calendar_age_range_BP = calendar_age_range_BP,
  rc_determinations = obs_radiocarbon_ages,
  rc_sigmas = obs_radiocarbon_sigmas,
  calibration_curve=intcal20,
  denscale = denscale, # Specify to ensure consistent plotting across panels
  plot_output = TRUE, plot_pretty = FALSE)

# Choose xlim and ylim to make plots as consistent as possible across panels
x_lim_plot <-  rev(range(SPD_initial_fit$calendar_age_BP)) # should be calendar_age_range_BP extended by 400 cal yrs

# Choose ylim to make consistent with automated SPD plots
cal_age_ind_min <- which.min(abs(intcal20$calendar_age_BP - min(x_lim_plot)))
cal_age_ind_max <- which.min(abs(intcal20$calendar_age_BP - max(x_lim_plot)))
calendar_age_indices <- cal_age_ind_min:cal_age_ind_max
calibration_curve_ub <- intcal20$c14_age + 2 * intcal20$c14_sig
calibration_curve_lb <- intcal20$c14_age - 2 * intcal20$c14_sig

# ylim for the 14C ages
y_lim_plot_c14_ages <- c(
  min(calibration_curve_lb[calendar_age_indices]),
  max(calibration_curve_ub[calendar_age_indices])
)
y_lim_plot_c14_ages <- y_lim_plot_c14_ages + 0.05 * c(-3, 1) * diff(y_lim_plot_c14_ages)

# The same but in terms of the calendar age density plots
y_lim_plot_density <- c(0, 3 * max(SPD_initial_fit$probability))

# Bit of tidying up
rm(calibration_curve_lb,
   calibration_curve_ub,
   calendar_age_indices,
   cal_age_ind_max,
   cal_age_ind_min)

## Plot A - Illustrate sampling 14C determinations  from a uniform phase
out_file_name <- "output/SPDFailure/SPDBootstrapFailure/BootstrapFailure/BootstrapFailure1"

if(pdf_output) {
  pdf(paste(out_file_name, ".pdf", sep = ""),
      width = 8,
      height = plot_height)
} else {
  png(paste(out_file_name, ".png", sep = ""),
      width = 8,
      height = plot_height,
      units = "in", res = 480)
}

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)


plot(intcal20$calendar_age_BP, intcal20$c14_age, col = "blue",
     ylab = expression(paste("Radiocarbon age (", ""^14, "C yr BP)")),
     xlab = "Calendar Age (cal yr BP)",
     xlim = x_lim_plot, ylim = y_lim_plot_c14_ages,
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
rug(obs_radiocarbon_ages, side = 2, ticksize = 0.03, lwd = 1, col = "red")

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
t <- seq(min_cal, max_cal, by = 1)
truerate <- rep(1, length(t))/length(t) # So integrates to 1

# Plot the true rate along the bottom
par(new = TRUE)
plot(c(min(t), t, max(t)), c(0,truerate,0), lty = 1, col = "red", type = "l",
     ylim = y_lim_plot_density,
     xlim = x_lim_plot,
     axes = FALSE,
     xlab = NA,
     ylab = NA,
     yaxs = "i")
dens_polygon <- cbind(c(min(t), t, max(t)),
                      c(0,truerate,0))

polygon(dens_polygon, col = grDevices::rgb(1, 0, 0, .5))


mtext(LETTERS[1], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)
dev.off()



#####################################################################################
# Plot B - showing the SPD of this density (doesn't look like the original)
out_file_name <- "output/SPDFailure/SPDBootstrapFailure/BootstrapFailure/BootstrapFailure2"

if(pdf_output) {
  pdf(paste(out_file_name, ".pdf", sep = ""),
      width = 8,
      height = plot_height)
} else {
  png(paste(out_file_name, ".png", sep = ""),
      width = 8,
      height = plot_height,
      units = "in", res = 480)
}

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)

SPD_initial_fit <- FindSummedProbabilityDistribution(
  calendar_age_range_BP = calendar_age_range_BP,
  rc_determinations = obs_radiocarbon_ages,
  rc_sigmas = obs_radiocarbon_sigmas,
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
  size = n_samp,
  replace = TRUE,
  prob = SPD_initial_fit$probability)

resampled_radiocarbon_ages <- rnorm(n_samp,
                                    mean = approx(intcal20$calendar_age_BP,
                                                  intcal20$c14_age,
                                                  resampled_calendar_ages)$y,
                                    sd = sqrt(obs_radiocarbon_sigmas^2 +
                                                (approx(intcal20$calendar_age_BP,
                                                        intcal20$c14_sig,
                                                        resampled_calendar_ages)$y)^2) )
rug(resampled_calendar_ages, col = "purple")


mtext(LETTERS[2], side = 3, adj = lab_adj, cex = panel_label_cex, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)



# Add rug to colour obs c14ages to make consistent with panel A
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)
plot(
  NULL,
  NULL,
  type = "n",
  ylim = y_lim_plot_c14_ages,
  xlim = x_lim_plot,
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  xaxs = "i",
  yaxs = "i")
rug(obs_radiocarbon_ages, side = 2, ticksize = 0.03, lwd = 1, col = "red")

dev.off()


##############################################################
# Plot C - SPD on resampled/bootstrap sample

# First make sure the density axis is the same for initial and bootstrap SPD
SPD_bootstrap_fit <- FindSummedProbabilityDistribution(
  calendar_age_range_BP = calendar_age_range_BP,
  rc_determinations = resampled_radiocarbon_ages,
  rc_sigmas = obs_radiocarbon_sigmas,
  calibration_curve = intcal20,
  plot_output = TRUE, plot_pretty = FALSE)
plot_denscale_bootstrap <- denscale * max(SPD_initial_fit$probability)/max(SPD_bootstrap_fit$probability)


out_file_name <- "output/SPDFailure/SPDBootstrapFailure/BootstrapFailure/BootstrapFailure3"

if(pdf_output) {
  pdf(paste(out_file_name, ".pdf", sep = ""),
      width = 8,
      height = plot_height)
} else {
  png(paste(out_file_name, ".png", sep = ""),
      width = 8,
      height = plot_height,
      units = "in", res = 480)
}

par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)

SPD_bootstrap_fit <- FindSummedProbabilityDistribution(
  calendar_age_range_BP = calendar_age_range_BP,
  rc_determinations = resampled_radiocarbon_ages,
  rc_sigmas = obs_radiocarbon_sigmas,
  calibration_curve = intcal20,
  denscale = plot_denscale_bootstrap,
  plot_output = TRUE, plot_pretty = FALSE)

par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)


plot(
  NULL,
  NULL,
  type = "n",
  ylim = y_lim_plot_density,
  xlim = x_lim_plot,
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

# Add rug to colour bootstrap resample c14ages
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)
plot(
  NULL,
  NULL,
  type = "n",
  ylim = y_lim_plot_c14_ages,
  xlim = x_lim_plot,
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  xaxs = "i",
  yaxs = "i")
rug(resampled_radiocarbon_ages, side = 2, col = "purple")

dev.off()




#############################
# Plot D - Create plot showing SPD bootstrap intervals


## First write a function to calculate supposed bootstrap CI for SPD
Not_WorkingFindBootstrapSPDInterval <- function(SPD_initial_fit,
                                     obs_radiocarbon_sigmas,
                                     calendar_age_range_BP,
                                     n_boot = 100,
                                     CI_probs = c(0.025, 0.975)) {

  n_samp <- length(obs_radiocarbon_sigmas)
  n_resamples_needed <- n_samp * n_boot

  # Create resamplesd calendar ages
  # We need n_resamples_needed which we will then convert to a matrix
  # where each row corresponds to a bootstrap resample i
  resampled_calendar_ages <- sample(
    x = SPD_initial_fit$calendar_age_BP,
    size = n_resamples_needed,
    replace = TRUE,
    prob = SPD_initial_fit$probability)

  # Mean of 14C determinations resampled from SPD
  resampled_radiocarbon_mean <- approx(intcal20$calendar_age_BP,
                                           intcal20$c14_age,
                                           resampled_calendar_ages)$y

  # Sd of determinations resampled from SPD
  # Combination of lab-reported/observed and calibration curve
  resampled_radiocarbon_sigmas <- sqrt(rep(obs_radiocarbon_sigmas, n_boot)^2 +
                                             (approx(intcal20$calendar_age_BP,
                                                     intcal20$c14_sig,
                                                     resampled_calendar_ages)$y)^2)


  resampled_radiocarbon_ages <- rnorm(n_resamples_needed,
                                      mean = resampled_radiocarbon_mean,
                                      sd = resampled_radiocarbon_sigmas)

  # Convert resampled 14C ages to matrix (each row is a single resample)
  bootstrap_resample_radiocarbon_ages <- matrix(resampled_calendar_ages,
                                                byrow = TRUE,
                                                nrow = n_boot)

  # Now create SPD for each row
  n_calendar_grid <- length(SPD_initial_fit$calendar_age_BP)
  bootstrap_SPD <- matrix(NA, ncol = n_calendar_grid, nrow = n_boot)

  for(i in 1:n_boot) {
    bootstrap_SPD[i,] <- FindSummedProbabilityDistribution(
      calendar_age_range_BP = calendar_age_range_BP,
      rc_determinations = bootstrap_resample_radiocarbon_ages[i,],
      rc_sigmas = obs_radiocarbon_sigmas,
      calibration_curve = intcal20,
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




# Function to calculate supposed bootstrap CI for SPD
FindBootstrapSPDInterval_MethodB <- function(SPD_initial_fit,
                                     obs_radiocarbon_sigmas,
                                     calendar_age_range_BP,
                                     n_boot = 100,
                                     CI_probs = c(0.025, 0.975)) {

  n_samp <- length(obs_radiocarbon_sigmas)

  # Create matrix that will contain bootstrap resampled SPD for each set of resamples
  n_calendar_grid <- length(SPD_initial_fit$calendar_age_BP)
  bootstrap_SPD <- matrix(NA, ncol = n_calendar_grid, nrow = n_boot)

  for(i in 1:n_boot) {
    resampled_calendar_ages <- sample(
      x = SPD_initial_fit$calendar_age_BP,
      size = n_samp,
      replace = TRUE,
      prob = SPD_initial_fit$probability)

    resampled_radiocarbon_ages <- rnorm(n_samp,
                                        mean = approx(intcal20$calendar_age_BP,
                                                      intcal20$c14_age,
                                                      resampled_calendar_ages)$y,
                                        sd = sqrt(obs_radiocarbon_sigmas^2 +
                                                    (approx(intcal20$calendar_age_BP,
                                                            intcal20$c14_sig,
                                                            resampled_calendar_ages)$y)^2) )

    Temp_bootstrap_SPD <- FindSummedProbabilityDistribution(
      calendar_age_range_BP = calendar_age_range_BP,
      rc_determinations = resampled_radiocarbon_ages,
      rc_sigmas = obs_radiocarbon_sigmas,
      calibration_curve = intcal20,
      plot_output = FALSE, plot_pretty = FALSE)

    if(!identical(range(Temp_bootstrap_SPD$calendar_age_BP), c(1500, 2600))) stop("Broken SPD")

    bootstrap_SPD[i,] <- Temp_bootstrap_SPD$probability
#   lines(Temp_bootstrap_SPD$calendar_age_BP, Temp_bootstrap_SPD$probability)
  }

  # Now find CI
  bootstrap_SPD_lower_CI <- apply(bootstrap_SPD, 2, quantile, probs = CI_probs[1])
  bootstrap_SPD_upper_CI <- apply(bootstrap_SPD, 2, quantile, probs = CI_probs[2])

  retlist <- list(calendar_age_BP = SPD_initial_fit$calendar_age_BP,
                  bootstrap_SPD_lower_CI = bootstrap_SPD_lower_CI,
                  bootstrap_SPD_upper_CI = bootstrap_SPD_upper_CI)
  return(retlist)
}



set.seed(18)
bootstrap_SPD_CI <- FindBootstrapSPDInterval_MethodB(SPD_initial_fit = SPD_initial_fit,
                                             obs_radiocarbon_sigmas = obs_radiocarbon_sigmas,
                                             calendar_age_range_BP = calendar_age_range_BP,
                                             n_boot = 10000)


# Make actual plot
out_file_name <- "output/SPDFailure/SPDBootstrapFailure/BootstrapFailure/BootstrapFailure4"

if(pdf_output) {
  pdf(paste(out_file_name, ".pdf", sep = ""),
      width = 8,
      height = plot_height)
} else {
  png(paste(out_file_name, ".png", sep = ""),
      width = 8,
      height = plot_height,
      units = "in", res = 480)
}



# Make actual plot
par(mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)

SPD_initial_fit <- FindSummedProbabilityDistribution(
  calendar_age_range_BP = calendar_age_range_BP,
  rc_determinations = obs_radiocarbon_ages,
  rc_sigmas = obs_radiocarbon_sigmas,
  calibration_curve = intcal20,
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

# Add rug to colour obs c14ages to make consistent with panel A
par(new = TRUE,
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 0.1, 2) + 0.1,
    las = 1)
plot(
  NULL,
  NULL,
  type = "n",
  ylim = y_lim_plot_c14_ages,
  xlim = x_lim_plot,
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  xaxs = "i",
  yaxs = "i")
rug(obs_radiocarbon_ages, side = 2, ticksize = 0.03, lwd = 1, col = "red")

dev.off()
