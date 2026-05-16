################################################################################
# Fig 5B – Overlay: Subjective Effect Profiles at Reference Dose
#
# Method: WLS linear regression (lm, weights = n) per drug × subscale,
#   predicting score at reference dose. Ketamine = n-weighted mean, all data.
# Reference doses: Psilocybin=25mg, LSD=100µg, MDMA=125mg, Ketamine=all data.
# Data source: ASDB v.2025-12-31
#
# This version uses base R graphics (not ggplot2) to match the matplotlib
# polar projection output exactly.
#
# Required packages: readxl (for data loading only)
# install.packages("readxl")
################################################################################

library(readxl)

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD & ANALYSE DATA
# ══════════════════════════════════════════════════════════════════════════════

ASDB_PATH <- "ASDB_v.2025-12-31.xlsx"   # <-- adjust path if needed

dat_5d  <- read_excel(ASDB_PATH, sheet = "5D-ASC",  na = "NULL")
dat_11  <- read_excel(ASDB_PATH, sheet = "11-ASC",  na = "NULL")
dat_meq <- read_excel(ASDB_PATH, sheet = "MEQ-30",  na = "NULL")

# ── Dose normalisation ───────────────────────────────────────────────────────
normalise_dose <- function(qty, unit, drug) {
  qty  <- suppressWarnings(as.numeric(qty))
  unit <- trimws(tolower(as.character(unit)))
  if (is.na(qty)) return(NA_real_)
  if (drug == "Psilocybin") {
    if (unit %in% c("mg", "mg "))                                  return(qty)
    if (unit %in% c("mg/kg", "mg/kg ", "mg/kg body weight"))       return(qty * 70)
    if (unit %in% c("\u00b5g/kg", "\u00b5g/kg body weight"))       return(qty * 0.07)
    if (unit %in% c("mg/70 kg", "mg/70kg"))                        return(qty)
    if (grepl("minutes", unit))                                    return(NA_real_)
    return(NA_real_)
  }
  if (drug == "LSD") {
    if (unit %in% c("\u00b5g (base)", "\u03bcg (base)", "\u00b5g")) return(qty)
    if (unit == "\u00b5g (tartrate)")                               return(qty / 1.46)
    return(NA_real_)
  }
  if (drug == "MDMA") {
    if (unit %in% c("mg", "mg "))                                  return(qty)
    if (unit %in% c("mg/kg", "mg/kg ", "mg/kg body weight"))       return(qty * 70)
    return(NA_real_)
  }
  return(qty)
}

# ── WLS prediction at reference dose ─────────────────────────────────────────
predict_wls <- function(dat, mean_col, drug, ref_dose) {
  suited_col <- grep("suited", names(dat), value = TRUE)[1]
  d <- dat[dat$induction_method_name == drug & dat[[suited_col]] == 1, ]
  d$dose_std <- mapply(normalise_dose, d$dosage_quantity, d$dosage_unit, drug)
  d$n        <- suppressWarnings(as.numeric(d$nr_of_subjects))
  d$score    <- suppressWarnings(as.numeric(d[[mean_col]]))
  d <- d[!is.na(d$dose_std) & !is.na(d$n) & d$n > 0 & !is.na(d$score), ]
  if (nrow(d) < 3) return(NA_real_)
  fit  <- lm(score ~ dose_std, data = d, weights = n)
  pred <- predict(fit, newdata = data.frame(dose_std = ref_dose))
  max(0, min(100, as.numeric(pred)))
}

# ── Ketamine: n-weighted mean, all data (no suited filter) ───────────────────
ketamine_wmean <- function(dat, mean_col) {
  d <- dat[dat$induction_method_name == "Ketamine", ]
  d$n     <- suppressWarnings(as.numeric(d$nr_of_subjects))
  d$score <- suppressWarnings(as.numeric(d[[mean_col]]))
  d <- d[!is.na(d$n) & d$n > 0 & !is.na(d$score), ]
  if (nrow(d) == 0) return(NA_real_)
  weighted.mean(d$score, d$n)
}

# ── Define axes ──────────────────────────────────────────────────────────────
dim_info <- data.frame(
  dim_id    = c("OB","VR","DED","VigR","Ins","BS","EU","SE","Anx","Myst","PM","Trans"),
  sheet     = c("5d","5d","5d","5d","11","11","11","11","11","meq","meq","meq"),
  mean_col  = c("Oceanic_Boundlessness_mean","Visionary_Restructuralization_mean",
                "Dread_of_Ego_Dissolution_mean","Vigilance_Reduction_mean",
                "Insightfulness_mean","Blissful_state_mean",
                "Experience_of_unity_mean","Spiritual_experience_mean",
                "Anxiety_mean","Mystical_mean","Positive_mood_mean",
                "Transcendence_of_time_and_space_mean"),
  label     = c("Oceanic\nBoundlessness","Visionary\nRestructuraliz.",
                "Dread of Ego\nDissolution","Vigilance\nReduction",
                "Insightfulness","Blissful\nState",
                "Experience\nof Unity","Spiritual\nExperience",
                "Anxiety","Mystical\n(MEQ-30)","Positive\nMood",
                "Transcendence\n(MEQ-30)"),
  quest_col = c(rep("#1a6b3c", 4), rep("#8B4513", 5), rep("#7b1fa2", 3)),
  stringsAsFactors = FALSE
)

sheet_map <- list("5d" = dat_5d, "11" = dat_11, "meq" = dat_meq)
REF_DOSE  <- c(Psilocybin = 25, LSD = 100, MDMA = 125)

# ── Compute predictions ─────────────────────────────────────────────────────
drugs <- c("Psilocybin", "LSD", "MDMA", "Ketamine")
scores <- list()

for (drug in drugs) {
  vals <- numeric(nrow(dim_info))
  for (i in seq_len(nrow(dim_info))) {
    dat <- sheet_map[[ dim_info$sheet[i] ]]
    mc  <- dim_info$mean_col[i]
    if (drug == "Ketamine") {
      vals[i] <- ketamine_wmean(dat, mc)
    } else {
      vals[i] <- predict_wls(dat, mc, drug, REF_DOSE[drug])
    }
  }
  # Replace NA with 0 for plotting (Ketamine has NA for Mystical and Positive Mood)
  vals[is.na(vals)] <- 0
  scores[[drug]] <- vals
}

# Print for verification
cat("\n-- Predicted scores at reference dose --\n")
result_table <- do.call(cbind, scores)
rownames(result_table) <- dim_info$dim_id
print(round(result_table, 1))

# ══════════════════════════════════════════════════════════════════════════════
# 2. DRAW THE RADAR CHART (base R graphics)
# ══════════════════════════════════════════════════════════════════════════════

# ── Drug colours ─────────────────────────────────────────────────────────────
drug_cols <- c(Psilocybin = "#009E73", LSD = "#E69F00",
               MDMA = "#CC79A7", Ketamine = "#56B4E9")

# ── Polar geometry ───────────────────────────────────────────────────────────
N <- nrow(dim_info)
# Angles: start at top (pi/2), go clockwise
angles <- pi/2 - seq(0, 2*pi, length.out = N + 1)[1:N]

# Polar to cartesian
pol2cart <- function(r, theta) {
  list(x = r * cos(theta), y = r * sin(theta))
}

# ── Helper: draw a closed polygon from radius values ─────────────────────────
draw_poly <- function(vals, theta, col_line, col_fill = NA,
                      lwd = 2.5, alpha_fill = 0.10, lty = 1) {
  pts <- pol2cart(c(vals, vals[1]), c(theta, theta[1]))
  if (!is.na(col_fill)) {
    fill_rgb <- col2rgb(col_fill) / 255
    polygon(pts$x, pts$y,
            col = rgb(fill_rgb[1], fill_rgb[2], fill_rgb[3], alpha_fill),
            border = NA)
  }
  lines(pts$x, pts$y, col = col_line, lwd = lwd, lty = lty)
}

# ── Helper: draw a reference ring ────────────────────────────────────────────
draw_ring <- function(r, col = "grey80", lwd = 0.5, lty = 1) {
  th <- seq(0, 2*pi, length.out = 300)
  p  <- pol2cart(r, th)
  lines(p$x, p$y, col = col, lwd = lwd, lty = lty)
}

# ── Helper: draw spoke lines from center to edge ────────────────────────────
draw_spokes <- function(theta, r_max = 100, col = "#dddddd", lwd = 0.5) {
  for (th in theta) {
    p <- pol2cart(r_max, th)
    segments(0, 0, p$x, p$y, col = col, lwd = lwd)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. RENDER
# ══════════════════════════════════════════════════════════════════════════════

# Save to PDF, PNG, SVG, TIFF
for (fmt in c("pdf", "png", "svg", "tiff")) {
  fname <- paste0("Fig5B_Radar_Overlay.", fmt)

  if (fmt == "pdf") {
    pdf(fname, width = 10, height = 10)
  } else if (fmt == "png") {
    png(fname, width = 10, height = 10, units = "in", res = 300)
  } else if (fmt == "svg") {
    svg(fname, width = 10, height = 10)
  } else if (fmt == "tiff") {
    tiff(fname, width = 10, height = 10, units = "in", res = 600,
         compression = "lzw")
  }

  # Set up plot area: cartesian square, no axes
  par(mar = c(3, 1, 4, 1), bg = "white", family = "sans")
  plot(NULL, xlim = c(-140, 140), ylim = c(-140, 140),
       asp = 1, axes = FALSE, xlab = "", ylab = "")

  # ── Spoke lines ──────────────────────────────────────────────────────────
  draw_spokes(angles, r_max = 100, col = "#dddddd", lwd = 0.5)

  # ── Reference rings ──────────────────────────────────────────────────────
  draw_ring(100, col = "#cccccc", lwd = 0.8, lty = 1)   # outer boundary
  draw_ring(75,  col = "grey80",  lwd = 0.5, lty = 3)   # dotted
  draw_ring(50,  col = "grey60",  lwd = 0.9, lty = 2)   # dashed (50% ref)
  draw_ring(25,  col = "grey80",  lwd = 0.5, lty = 3)   # dotted

  # ── Ring labels ──────────────────────────────────────────────────────────
  # Place along the first spoke (top), slightly offset right
  for (rv in c(25, 50, 75, 100)) {
    p <- pol2cart(rv, angles[1])
    text(p$x + 3, p$y + 3, paste0(rv, "%"),
         cex = 0.55, col = "grey55", adj = c(0, 0))
  }

  # ── Drug polygons ────────────────────────────────────────────────────────
  for (drug in drugs) {
    v   <- scores[[drug]]
    col <- drug_cols[drug]

    # Filled polygon
    draw_poly(v, angles, col_line = col, col_fill = col,
              lwd = 2.5, alpha_fill = 0.10)

    # Points at vertices
    pts <- pol2cart(v, angles)
    points(pts$x, pts$y, pch = 16, col = col, cex = 1.1)
    # White edge ring around each point
    points(pts$x, pts$y, pch = 1, col = "white", cex = 1.1, lwd = 0.8)
  }

  # ── Axis labels (coloured by questionnaire) ─────────────────────────────
  label_r <- 113   # radial distance for labels
  for (i in seq_len(N)) {
    p <- pol2cart(label_r, angles[i])

    # Determine text alignment based on angle position
    angle_deg <- (angles[i] * 180 / pi) %% 360
    if (angle_deg > 80 & angle_deg < 100) {
      adj <- c(0.5, 0)      # top
    } else if (angle_deg > 260 & angle_deg < 280) {
      adj <- c(0.5, 1)      # bottom
    } else if (angle_deg >= 100 & angle_deg <= 260) {
      adj <- c(1, 0.5)      # left half
    } else {
      adj <- c(0, 0.5)      # right half
    }

    text(p$x, p$y, dim_info$label[i],
         col = dim_info$quest_col[i], font = 2, cex = 0.85,
         adj = adj)
  }

  # ── Legend ───────────────────────────────────────────────────────────────
  legend("topright", inset = c(-0.02, -0.02),
         legend = drugs,
         fill   = drug_cols[drugs],
         border = drug_cols[drugs],
         title  = "Drug",
         title.font = 2,
         cex    = 0.95,
         bty    = "o",
         box.col = "grey80",
         bg     = "white")

  # ── Title ────────────────────────────────────────────────────────────────
  title(main = expression(bold("Overlay: Subjective Effect Profiles at Reference Dose")),
        line = 2.5, cex.main = 1.15)
  mtext("Psilocybin=25mg \u00b7 LSD=100\u00b5g \u00b7 MDMA=125mg \u00b7 Ketamine=all data",
        side = 3, line = 1.0, cex = 0.85, col = "grey30")

  # ── Caption ──────────────────────────────────────────────────────────────
  mtext(paste0(
    "Lines = WLS-predicted score at reference dose. ",
    "Axis label colours: \u25a0 green = 5D-ASC \u00b7 ",
    "\u25a0 brown = 11-ASC \u00b7 \u25a0 purple = MEQ-30. ",
    "Dashed ring = 50% of maximum. Source: ASDB v.2025-12-31."),
    side = 1, line = 1.5, cex = 0.6, col = "grey50")

  dev.off()
  message("Saved: ", fname)
}

message("\nDone - all formats saved.")
