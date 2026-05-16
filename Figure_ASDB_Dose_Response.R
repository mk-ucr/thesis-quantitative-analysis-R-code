## ============================================================
## Dose-Response Figure from ASDB dataset
## Panel A: 11D-ASC Insightfulness (LSD vs Psilocybin)
## Panel B: MEQ-30 Mystical (LSD vs Psilocybin)
##
## Doses are standardised:
##   LSD  -> ug base (tartrate * 0.813)
##   Psilocybin -> mg absolute (mg/kg * 70, ug/kg * 70/1000)
##
## USAGE: Place this script in the same folder as the ASDB
##        Excel file, then source() or run it in RStudio.
## ============================================================

library(readxl)
library(dplyr)

## ---- 0. Locate the ASDB file robustly ----
## Instead of hard-coding the filename (which can fail due to
## encoding mismatches between the script and the OS), we scan
## for any file matching "ASDB*.xlsx" in the working directory.

candidates <- list.files(
  path        = getwd(),
  pattern     = "^ASDB.*\\.xlsx$",
  full.names  = TRUE,
  ignore.case = TRUE
)

if (length(candidates) == 0) {
  stop(
    "Could not find an ASDB .xlsx file in the working directory:\n  ",
    getwd(),
    "\nPlease make sure the file is there, or use setwd() first."
  )
}

file <- normalizePath(candidates[1], mustWork = TRUE)
cat("Using data file:", file, "\n\n")

## ---- 1. Dose standardisation helpers ----

standardize_lsd_ug <- function(qty, unit) {
  dose <- suppressWarnings(as.numeric(qty))
  if (is.na(dose)) return(NA_real_)
  unit <- trimws(tolower(as.character(unit)))
  if (grepl("tartrate", unit))   return(dose * 0.813)
  if (grepl("base", unit))       return(dose)
  if (grepl("^\\s*\u00b5g\\s*$", unit) ||
      grepl("^\\s*\u03bcg\\s*$", unit)) return(dose)
  NA_real_
}

standardize_psy_mg <- function(qty, unit) {
  dose <- suppressWarnings(as.numeric(qty))
  if (is.na(dose)) return(NA_real_)
  unit <- trimws(tolower(as.character(unit)))
  if (unit %in% c("mg", "mg/70kg"))       return(dose)
  if (grepl("mg/kg", unit))                return(dose * 70)
  if (grepl("g/kg", unit))                 return(dose * 70 / 1000)
  if (unit == "g")                         return(dose * 1000)
  NA_real_
}

## ---- 2. Generic loader ----

load_sheet <- function(file, sheet, score_col) {
  df <- read_excel(file, sheet = sheet)

  ctrl_col <- grep("is_control", names(df), value = TRUE)[1]
  dr_col   <- grep("suited_for_dose", names(df), value = TRUE)[1]

  df <- df %>%
    filter(.data[[dr_col]] == 1, .data[[ctrl_col]] == 0) %>%
    mutate(
      Score  = as.numeric(.data[[score_col]]),
      method = trimws(as.character(induction_method_name)),
      n      = as.numeric(nr_of_subjects)
    ) %>%
    filter(!is.na(Score))

  lsd <- df %>%
    filter(method == "LSD") %>%
    rowwise() %>%
    mutate(Dose_std = standardize_lsd_ug(dosage_quantity, dosage_unit)) %>%
    ungroup() %>%
    filter(!is.na(Dose_std))

  psy <- df %>%
    filter(method == "Psilocybin") %>%
    rowwise() %>%
    mutate(Dose_std = standardize_psy_mg(dosage_quantity, dosage_unit)) %>%
    ungroup() %>%
    filter(!is.na(Dose_std))

  list(lsd = lsd, psy = psy)
}

## ---- 3. Load data ----

asc <- load_sheet(file, "11-ASC", "Insightfulness_mean")
meq <- load_sheet(file, "MEQ-30", "Mystical_mean")

cat("Panel A  11D-ASC Insightfulness\n")
cat(sprintf("  LSD: n = %d, r = %.2f\n", nrow(asc$lsd),
            cor(asc$lsd$Score, asc$lsd$Dose_std, use = "complete.obs")))
cat(sprintf("  Psy: n = %d, r = %.2f\n", nrow(asc$psy),
            cor(asc$psy$Score, asc$psy$Dose_std, use = "complete.obs")))
cat("Panel B  MEQ-30 Mystical\n")
cat(sprintf("  LSD: n = %d, r = %.2f\n", nrow(meq$lsd),
            cor(meq$lsd$Score, meq$lsd$Dose_std, use = "complete.obs")))
cat(sprintf("  Psy: n = %d, r = %.2f\n", nrow(meq$psy),
            cor(meq$psy$Score, meq$psy$Dose_std, use = "complete.obs")))

## ---- 4. Plot ----

col_lsd   <- "#9B2335"
col_psy   <- "#2D5F8A"
col_axis  <- "#2B2B2B"
col_lsd_a <- adjustcolor(col_lsd, alpha.f = 0.70)
col_psy_a <- adjustcolor(col_psy, alpha.f = 0.70)

out_file <- file.path(getwd(), "Figure_ASDB_Dose_Response.png")
png(out_file,
    width = 174, height = 110, units = "mm", res = 300, bg = "white")

layout(matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE),
       heights = c(1, 0.1))
par(oma = c(2, 0.5, 1, 0.5), family = "sans")

draw_dual <- function(dat_lsd, dat_psy, y_label, panel_letter,
                      ylim = c(0, 100), size_div = 4) {
  par(mar = c(3.5, 4, 3.5, 1))

  plot(dat_lsd$Dose_std, dat_lsd$Score,
       type = "n", axes = FALSE, xlab = "", ylab = "",
       ylim = ylim,
       xlim = range(dat_lsd$Dose_std, na.rm = TRUE) * c(0.9, 1.05))

  yticks <- pretty(ylim, n = 5)
  abline(h = yticks, col = "#ECECEC", lwd = 0.5)

  if (nrow(dat_lsd) >= 2) {
    abline(lm(Score ~ Dose_std, data = dat_lsd),
           col = col_lsd, lty = 2, lwd = 1.1)
  }
  cex_lsd <- sqrt(dat_lsd$n) / size_div
  cex_lsd <- pmax(cex_lsd, 0.4)
  points(dat_lsd$Dose_std, dat_lsd$Score,
         pch = 21, cex = cex_lsd, bg = col_lsd_a,
         col = "white", lwd = 0.4)

  axis(1, col = col_axis, col.axis = col_axis, tcl = -0.3,
       cex.axis = 0.7, mgp = c(3, 0.4, 0))
  axis(2, at = yticks, col = col_axis, col.axis = col_axis,
       las = 1, tcl = -0.3, cex.axis = 0.7, mgp = c(3, 0.5, 0))
  mtext(expression(bold("LSD dose") ~ "(\u00b5g base)"),
        side = 1, line = 1.8, cex = 0.55, col = col_lsd)
  mtext(y_label, side = 2, line = 2.5, cex = 0.55,
        col = col_axis, font = 2)

  par(new = TRUE)
  plot(dat_psy$Dose_std, dat_psy$Score,
       type = "n", axes = FALSE, xlab = "", ylab = "",
       ylim = ylim,
       xlim = range(dat_psy$Dose_std, na.rm = TRUE) * c(0.9, 1.05))

  if (nrow(dat_psy) >= 2) {
    abline(lm(Score ~ Dose_std, data = dat_psy),
           col = col_psy, lty = 2, lwd = 1.1)
  }
  cex_psy <- sqrt(dat_psy$n) / size_div
  cex_psy <- pmax(cex_psy, 0.4)
  points(dat_psy$Dose_std, dat_psy$Score,
         pch = 22, cex = cex_psy, bg = col_psy_a,
         col = "white", lwd = 0.4)

  axis(3, col = col_axis, col.axis = col_axis, tcl = -0.3,
       cex.axis = 0.7, mgp = c(3, 0.4, 0))
  mtext(expression(bold("Psilocybin dose") ~ "(mg)"),
        side = 3, line = 1.8, cex = 0.55, col = col_psy)

  r_psy <- cor(dat_psy$Score, dat_psy$Dose_std, use = "complete.obs")
  r_lsd <- cor(dat_lsd$Score, dat_lsd$Dose_std, use = "complete.obs")
  usr <- par("usr")
  text(usr[1] + 0.03 * diff(usr[1:2]), usr[4] - 0.04 * diff(usr[3:4]),
       bquote(italic(r) == .(sprintf("%.2f", r_psy))),
       col = col_psy, cex = 0.6, adj = c(0, 1))
  text(usr[1] + 0.03 * diff(usr[1:2]), usr[4] - 0.12 * diff(usr[3:4]),
       bquote(italic(r) == .(sprintf("%.2f", r_lsd))),
       col = col_lsd, cex = 0.6, adj = c(0, 1))

  box(bty = "l", col = col_axis, lwd = 0.6)
  mtext(panel_letter, side = 3, line = 2.5, adj = 0,
        cex = 0.85, font = 2, col = col_axis)
}

draw_dual(asc$lsd, asc$psy, "11D-ASC Insightfulness (%)", "A")
draw_dual(meq$lsd, meq$psy, "MEQ-30 Mystical (%)", "B")

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("LSD", "Psilocybin"),
       col    = c(col_lsd, col_psy),
       pt.bg  = c(col_lsd_a, col_psy_a),
       pch    = c(21, 22),
       pt.cex = 1.2,
       lty    = 2, lwd = 1.1,
       bty    = "n", horiz = TRUE,
       cex    = 0.75, text.col = col_axis,
       seg.len = 2, x.intersp = 0.8)

mtext("Point size proportional to sample size (n). LSD doses: \u00b5g base; Psilocybin doses: mg (70 kg assumed).",
      side = 1, outer = TRUE, line = 0.5, cex = 0.5, col = "grey50")

dev.off()
cat("\nFigure saved as:", out_file, "\n")
