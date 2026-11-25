# Conceptual mean curves, not simulated observations or fitted results.
# Run from any directory: Rscript scripts/mm_model_concept.R [output.pdf]
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
bundle <- dirname(dirname(script))
output <- if (length(args)) args[1L] else
  file.path(bundle, "figures", "mm_model_concept_rev1.pdf")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

draw_mm_concept <- function() {
  ink <- "#243746"
  teal <- "#007F82"
  orange <- "#B85A28"
  guide <- "#87939B"
  pale <- "#E3F0EF"
  s <- seq(0, 8.4, length.out = 1001L)
  mm <- function(s, scale = 1) scale * s / (1 + s)
  stopifnot(mm(1) == 0.5, mm(1, 1.25) == 1.25 / 2,
            mm(1, 0.75) == 0.75 / 2)

  par(mfrow = c(1, 2), mar = c(3.8, 4.7, 3.1, 0.7),
      mgp = c(2.5, 0.7, 0), tcl = -0.24,
      col = ink, col.axis = ink, col.lab = ink,
      family = "Helvetica", las = 1, xaxs = "i", yaxs = "i")
  start_panel <- function(title, subtitle, y_label = FALSE) {
    plot(NA_real_, NA_real_, type = "n", xlim = c(0, 8.4),
         ylim = c(0, 1.34), axes = FALSE, xlab = "", ylab = "")
    axis(1, at = c(0, 1, 4, 8),
         labels = expression(0, K[m], 4 * K[m], 8 * K[m]),
         lwd = 0, lwd.ticks = 0.65, cex.axis = 0.94)
    axis(2, at = c(0, 0.5, 1),
         labels = expression(0, V[plain(max)] / 2, V[plain(max)]),
         lwd = 0, lwd.ticks = 0.65, cex.axis = 0.94)
    segments(0, 0, 8.4, 0, col = ink, lwd = 0.7)
    segments(0, 0, 0, 1.34, col = ink, lwd = 0.7)
    mtext("Substrate concentration, S", side = 1, line = 2.25, cex = 0.94)
    if (y_label) mtext("Mean reaction velocity", side = 2, line = 3.5,
                      las = 0, cex = 0.94)
    mtext(title, side = 3, line = 1.6, adj = 0, cex = 1, font = 2)
    mtext(subtitle, side = 3, line = 0.35, adj = 0, cex = 0.86,
          col = "#596873")
  }

  start_panel("(a) Single MM curve",
              expression(paste("Half-saturation identifies ", K[m])), TRUE)
  polygon(c(0, 1, 1, 0), c(0, 0, 0.5, 0.5), col = pale, border = NA)
  segments(0, 1, 8.4, 1, col = guide, lty = "33", lwd = 0.9)
  segments(0, 0.5, 1, 0.5, col = teal, lty = "33", lwd = 1)
  segments(1, 0, 1, 0.5, col = teal, lty = "33", lwd = 1)
  lines(s, mm(s), col = teal, lwd = 2.5)
  points(1, 0.5, pch = 21, bg = teal, col = "white", lwd = 1, cex = 1.1)
  text(4.45, 1.065, "Limiting velocity", col = "#596873", cex = 0.85)
  text(2.0, 0.38, expression(paste("At ", K[m], ": ", mu == V[plain(max)] / 2)),
       adj = 0, col = teal, cex = 0.9)

  start_panel("(b) Cluster-specific curves",
              expression(paste("Shared ", K[m], ", different plateaus")))
  segments(0, 1, 8.4, 1, col = guide, lty = "33", lwd = 0.9)
  segments(1, 0, 1, 0.625, col = guide, lty = "33", lwd = 0.9)
  lines(s, mm(s, 1.25), col = orange, lwd = 2.3, lty = "23")
  lines(s, mm(s), col = ink, lwd = 1.8, lty = 1)
  lines(s, mm(s, 0.75), col = teal, lwd = 2.3, lty = "62")
  points(rep(1, 3), c(0.625, 0.5, 0.375), pch = 21,
         bg = c(orange, ink, teal), col = "white", lwd = 0.8, cex = 0.85)
  text(4.1, 1.17, expression(paste("Cluster 1: ", b[1] > 0)),
       adj = 0, col = orange, cex = 0.85)
  text(0.35, 1.065, "Population plateau", adj = 0,
       col = "#596873", cex = 0.85)
  text(4.1, 0.91, "Population mean", adj = 0, col = ink, cex = 0.85)
  text(4.1, 0.51, expression(paste("Cluster 2: ", b[2] < 0)),
       adj = 0, col = teal, cex = 0.85)
}

raw_pdf <- tempfile(fileext = ".pdf")
pdf(raw_pdf, width = 7.4, height = 3.45, pointsize = 11,
    family = "Helvetica", useDingbats = FALSE,
    title = "Conceptual Michaelis-Menten mean models")
draw_mm_concept()
invisible(dev.off())
embedFonts(raw_pdf, outfile = output,
           options = c("-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true",
                       "-dSubsetFonts=true"))
unlink(raw_pdf)
message("Wrote ", normalizePath(output))
