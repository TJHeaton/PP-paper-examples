# Illustrate the failings of SPDs

# Performs two analyses:
# Analysis 1: Simulates samples from a calendar age model that is a mixture of two normals
# - Show that SPDs are overly variable, overly spread and cannot be easily interpreted.
# - Demonstrates that the DPMM (Heaton, 2022) is able to reconsruct the underlying calndar age model
# This forms Fig 2

# Analysis 2: Calibrates SPD from a single 14C determination
# - Show that SPDs paradigm of multiple periods of activity is evidently false
# This forms Figure 3

# Load carbondate library
library(carbondate)


# Plotting parameters
plot_height <- 4
top_margin <- 2.5
lab_adj <- 0.02
panel_label_cex <- 4/3

# Decide if want plotting output as pdf (TRUE) or png (FALSE)
pdf_output <- FALSE

# Decide if write plots to file
write_plots_to_file <- FALSE

# Store user par parameters (and revert back to these at the end)
oldpar <- par(no.readonly = TRUE)

##############################################################
#### Analysis 1 - Recosntructing a mixture of two normals
#### Set of two plots - the first is used as Fig 2
#### Plot i) - Show SPD is overly variable and interpretation hard
#### Plot ii)  - Show Polya Urn DPMM reconstruction works (not used in paper)
##############################################################

# Plot i) - Show that SPDs overly variable and interpretation hard
out_file_name <- "output/SPDFailure/SPDOverlyVariable"

# Decide if write plots to a file
if(write_plots_to_file) {
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
}

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

if(write_plots_to_file) {
  dev.off()
}


#####
# Plot ii) - show how a DPMM works to recreate correct calendar age density
out_file_name <- "output/SPDFailure/DPMMFit"

# Decide if write plots to a file
if(write_plots_to_file) {
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

if(write_plots_to_file) {
  dev.off()
}




##############################################################
#### Analysis 2
#### Single plot to show single determination SPD doesn't make sense
#### there needs to be some modelling in the calendar age domain
##############################################################


out_file_name <- "output/SPDFailure/SPDSingleDetermination"

# Decide if write plots to a file
if(write_plots_to_file) {
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

if(write_plots_to_file) {
  dev.off()
}


# Return to user specified par parameters (before running code)
par(oldpar)

