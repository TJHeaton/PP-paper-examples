# Plot the humans 14C data and reconstruction against the IntCal20 curve
curve_interval_sigma <- 2

# Plotting calendar age limits for clarity
plot_human_mincal <- 9000
plot_human_maxcal <- 16000

## Select file name for output
out_file_name <- "Plots/IllustrateHumans14CDeterminatons"

if(pdf_output) {
  pdf(paste(out_file_name, ".pdf", sep = ""),
      width = 7,
      height = 5.03)
} else {
  png(paste(out_file_name, ".png", sep = ""),
      width = 7,
      height = 5.03,
      units = "in", res = 480)
}


par(mar = c(3, 4.5, 0.3, 2) + 0.1, las = 1)
plot(intcal20$calendar_age_BP, intcal20$c14_age, col = "blue",
     ylim = range(Human14C$C14age) + c(-1.5,1.5) * max(Human14C$C14sig),
     xlim = c(plot_human_maxcal, plot_human_mincal) + c(100, -100),
     xlab = "",
     ylab = expression(paste("Radiocarbon age (", ""^14, "C yr BP)")),
     type = "l", main = "")

# Add ylabel and smaller ticks
mtext("Calendar Age (cal yr BP)", side = 1, line = 2)
xtick<-seq(8000, 17000, by=200)
ytick<-seq(5000, 16000, by=200)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)


# Add some question marks along the bottom
n_rep <- 20
text(rep("?", n_rep),
     x = seq(9000, 16000, length = n_rep),
     y = 7550, font = 2)

# Plot calibration curve
intcal20$ub <- intcal20$c14_age + curve_interval_sigma * intcal20$c14_sig
intcal20$lb <- intcal20$c14_age - curve_interval_sigma * intcal20$c14_sig
lines(intcal20$calendar_age_BP, intcal20$ub, lty = 3, col = "blue" )
lines(intcal20$calendar_age_BP, intcal20$lb, lty = 3, col = "blue")
polygon(c(rev(intcal20$calendar_age_BP),
          intcal20$calendar_age_BP),
        c(rev(intcal20$lb), intcal20$ub),
        col=rgb(0,0,1,.3),
        border=NA)
rug(Human14C$C14age, side = 2, ticksize = 0.03, lwd = 1, col = "red")


legend_labels <- c(
  substitute(paste(""^14, "C determination ")),
  "IntCal20",
  expression(paste("2", sigma, " interval")))
lty <- c(1, 1, 2)
lwd <- c(1, 1, 1)
pch <- c(NA, NA, NA)
col <- c(grDevices::rgb(1, 0, 0, .5), "blue", "blue")
legend("topright", legend = legend_labels, lty = lty, lwd=lwd, pch = pch, col = col)

# Add panel label
mtext(LETTERS[2], side = 3, adj = 0.01, cex = 1.2, font = 2,
      line = 3.4/3 * (-1), outer = FALSE)

dev.off()

