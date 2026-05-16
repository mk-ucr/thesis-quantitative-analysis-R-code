## ===========================================================
##  Figure & Summary Table: CAPS-Based PTSD Outcomes
##  MDMA-assisted therapy and ketamine RCTs
##
##  Data source : CAPS_PTSD_Outcomes_final.xlsx (every value
##                directly verified against the source papers;
##                see sheet Derivation_Notes)
##
##  Outputs     : Figure_CAPS_Outcomes.png / .pdf
##                Table_Demographics.png  / .pdf
##                stats_table.csv (reproduced statistical tests)
##
##  Author      : <your name>
##  Last update : 2026
## ===========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(readxl)
  library(scales)
  library(grid)
  library(gridExtra)
})

## ----------------------------------------------------------------
## 0.  Palette and theme
## ----------------------------------------------------------------
COL_ACTIVE_MDMA  <- "#8E3E6E"
COL_PLACEBO_MDMA <- "#D98CB3"
COL_ACTIVE_KET   <- "#2471A3"
COL_PLACEBO_KET  <- "#76BEDE"

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line          = element_line(colour = "black", linewidth = 0.4),
      axis.ticks         = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length  = unit(2.5, "pt"),
      axis.text          = element_text(colour = "black", size = base_size - 1),
      axis.title         = element_text(colour = "black", size = base_size,
                                        face = "bold"),
      axis.text.x        = element_text(angle = 0, hjust = 0.5),
      panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95", colour = NA),
      strip.text         = element_text(size = base_size - 0.5, face = "bold"),
      legend.position    = "right",
      legend.title       = element_text(size = base_size - 1, face = "bold"),
      legend.text        = element_text(size = base_size - 1.5),
      legend.key.size    = unit(8, "pt"),
      plot.margin        = margin(5, 8, 5, 5, "pt"),
      plot.title         = element_text(size = base_size + 1, face = "bold",
                                        hjust = 0),
      plot.subtitle      = element_text(size = base_size - 1.5,
                                        colour = "grey35", hjust = 0)
    )
}

## ----------------------------------------------------------------
## 1.  Load data
## ----------------------------------------------------------------
DATA_FILE <- "CAPS_PTSD_Outcomes_final.xlsx"
caps_raw  <- read_excel(DATA_FILE, sheet = "CAPS_Raw_Data")
demo_raw  <- read_excel(DATA_FILE, sheet = "Demographics")

label_map <- c(
  "Mitchell2021"      = "Mitchell\n2021",
  "Mitchell2023"      = "Mitchell\n2023",
  "Mithoefer2018_75"  = "Mithoefer\n2018\n(75 mg)",
  "Mithoefer2018_125" = "Mithoefer\n2018\n(125 mg)",
  "Feder2014"         = "Feder\n2014",
  "Feder2021"         = "Feder\n2021"
)

caps <- caps_raw |>
  rename(
    study_id      = `Study ID`,
    drug          = Drug,
    caps_version  = `CAPS Version`,
    placebo_type  = `Placebo Type`,
    n_active      = `n Active`,
    n_placebo     = `n Placebo`,
    bl_active     = `Baseline Active (mean)`,
    bl_placebo    = `Baseline Placebo (mean)`,
    bl_active_sd  = `Baseline Active (SD)`,
    bl_placebo_sd = `Baseline Placebo (SD)`,
    ch_active     = `Change Active (mean)`,
    ch_active_sd  = `Change Active (SD)`,
    ch_active_se  = `Change Active SE`,
    ch_placebo    = `Change Placebo (mean)`,
    ch_placebo_sd = `Change Placebo (SD)`,
    ch_placebo_se = `Change Placebo SE`,
    between_diff  = `Between-Group Diff`,
    between_se    = `Between-Group SE`,
    ci_lo_raw     = `Between CI Low`,
    ci_hi_raw     = `Between CI High`,
    cohen_d       = `Cohen's d`,
    d_lo          = `d CI Low`,
    d_hi          = `d CI High`,
    resp_active   = `Response Rate Active (%)`,
    resp_placebo  = `Response Rate Placebo (%)`,
    nodiag_active = `Loss of Diagnosis Active (%)`,
    nodiag_placebo= `Loss of Diagnosis Placebo (%)`,
    p_val         = `p value`
  ) |>
  mutate(
    label             = label_map[study_id],
    norm_factor       = ifelse(caps_version == "CAPS-IV", 80 / 136, 1),
    ch_active_norm    = ch_active   * norm_factor,
    ch_placebo_norm   = ch_placebo  * norm_factor,
    ch_active_se_n    = ch_active_se  * norm_factor,
    ch_placebo_se_n   = ch_placebo_se * norm_factor,
    between_norm      = between_diff * norm_factor,
    between_se_norm   = between_se   * norm_factor,
    ci_lo_norm        = ci_lo_raw    * norm_factor,
    ci_hi_norm        = ci_hi_raw    * norm_factor,
    study_label       = factor(label, levels = label),
    drug              = factor(drug, levels = c("MDMA", "Ketamine")),
    sig_label = case_when(
      p_val == "<0.0001" ~ "***",
      p_val == "<0.001"  ~ "***",
      p_val == "=0.0005" ~ "***",
      p_val == "=0.004"  ~ "**",
      p_val == "=0.20"   ~ "ns",
      TRUE               ~ ""
    )
  )

## ----------------------------------------------------------------
## 2.  Statistical re-analysis (published-summary t-tests + 2x2 chi-squared)
## ----------------------------------------------------------------
##
## RATIONALE
## ---------
## For every study where the source paper reports group means and SDs of
## CAPS change scores (Mitchell 2021, Mithoefer 2018 both doses), we run
## an independent-samples Welch t-test using the published summary
## statistics. This is an approximation -- we do not have individual
## patient data -- but it allows transparent recovery of test statistics,
## degrees of freedom and Hedges-corrected effect sizes from the
## reported summary data.
##
## For studies that report categorical outcomes (response, loss of
## diagnosis), we run 2 x 2 chi-squared tests on the published counts.
##
## Where the source paper reports an LS-mean difference from MMRM /
## mixed-effects ITT with an SE or 95% CI (Mitchell 2021/2023,
## Feder 2021), we accept the published Wald-style test statistic
## (Z = diff / SE) and use the published 95% CI directly. Mithoefer 2018
## did not report an SE; we derive it from Cohen's d and its 95% CI
## (see Derivation_Notes in the data file).
##
## Significance markers used in the figure:
##   ***  P < 0.001     **  P < 0.01     *  P < 0.05     ns  not significant

t_from_summary <- function(m1, sd1, n1, m2, sd2, n2) {
  if (any(is.na(c(m1, sd1, n1, m2, sd2, n2)))) return(NULL)
  se   <- sqrt(sd1^2 / n1 + sd2^2 / n2)
  tval <- (m1 - m2) / se
  df_w <- (sd1^2 / n1 + sd2^2 / n2)^2 /
          ((sd1^2 / n1)^2 / (n1 - 1) + (sd2^2 / n2)^2 / (n2 - 1))
  p    <- 2 * pt(-abs(tval), df = df_w)
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  d_est     <- (m1 - m2) / pooled_sd
  list(t = tval, df = df_w, p = p, d = d_est, se_diff = se)
}

stat_table <- caps |>
  rowwise() |>
  mutate(
    t_result = list(t_from_summary(ch_active, ch_active_sd, n_active,
                                   ch_placebo, ch_placebo_sd, n_placebo))
  ) |>
  ungroup() |>
  mutate(
    t_value      = sapply(t_result, function(x) if (is.null(x)) NA else round(x$t, 3)),
    df_welch     = sapply(t_result, function(x) if (is.null(x)) NA else round(x$df, 2)),
    p_recomputed = sapply(t_result, function(x) if (is.null(x)) NA else signif(x$p, 3)),
    d_recomputed = sapply(t_result, function(x) if (is.null(x)) NA else round(x$d, 2)),
    se_diff_t    = sapply(t_result, function(x) if (is.null(x)) NA else round(x$se_diff, 3))
  ) |>
  select(study_id, drug, n_active, n_placebo,
         ch_active, ch_active_sd, ch_placebo, ch_placebo_sd,
         t_value, df_welch, p_recomputed, d_recomputed, se_diff_t,
         published_p          = p_val,
         published_d          = cohen_d,
         published_between_dif = between_diff,
         published_between_se  = between_se)

## ---- 2x2 chi-square tests for response and loss-of-diagnosis -----
chi_categorical <- function(p_act, n_act, p_plc, n_plc) {
  if (any(is.na(c(p_act, n_act, p_plc, n_plc)))) return(NULL)
  k_act <- round(p_act * n_act / 100)
  k_plc <- round(p_plc * n_plc / 100)
  tab   <- matrix(c(k_act, n_act - k_act, k_plc, n_plc - k_plc),
                  nrow = 2, byrow = TRUE)
  test  <- suppressWarnings(chisq.test(tab, correct = FALSE))
  list(chi2 = unname(test$statistic),
       df   = unname(test$parameter),
       p    = test$p.value,
       k_act = k_act, k_plc = k_plc)
}

categorical_table <- caps |>
  rowwise() |>
  mutate(
    resp_chi = list(chi_categorical(resp_active,   n_active, resp_placebo,   n_placebo)),
    diag_chi = list(chi_categorical(nodiag_active, n_active, nodiag_placebo, n_placebo))
  ) |>
  ungroup() |>
  mutate(
    chi2_resp = sapply(resp_chi, function(x) if (is.null(x)) NA else round(x$chi2, 3)),
    p_resp    = sapply(resp_chi, function(x) if (is.null(x)) NA else signif(x$p, 3)),
    chi2_diag = sapply(diag_chi, function(x) if (is.null(x)) NA else round(x$chi2, 3)),
    p_diag    = sapply(diag_chi, function(x) if (is.null(x)) NA else signif(x$p, 3))
  ) |>
  select(study_id, drug,
         resp_active, resp_placebo, chi2_resp, p_resp,
         nodiag_active, nodiag_placebo, chi2_diag, p_diag)

## Write stats CSVs
write.csv(stat_table,        "stats_table_change_scores.csv", row.names = FALSE)
write.csv(categorical_table, "stats_table_categorical.csv",  row.names = FALSE)

cat("\n==== Reproduced t-tests on published change-score summaries ====\n")
print(as.data.frame(stat_table))
cat("\n==== Chi-squared tests on published categorical outcomes ====\n")
print(as.data.frame(categorical_table))

## Helper for asterisk labels from p-value -------------------
ast <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

## Apply chi-squared p-values to caps for use in panels C/D
caps <- caps |>
  left_join(categorical_table |>
              select(study_id, p_resp, p_diag),
            by = "study_id") |>
  mutate(
    sig_resp = sapply(p_resp, ast),
    sig_diag = sapply(p_diag, ast)
  )

## ----------------------------------------------------------------
## 3.  Panel A - mean CAPS change with within-group SE error bars
## ----------------------------------------------------------------
##
## Error bars represent within-group standard error of the mean
## change score (SE = SD / sqrt(n)) where SDs are reported in the
## source paper. Where SDs are not reported (Mitchell 2023 LS-means
## change scores), error bars are omitted (n.r. annotation).

## ----------------------------------------------------------------
## 3.  Panel A - mean CAPS change, with selectable error type
## ----------------------------------------------------------------
##
## error_type can be:
##   "sd" - within-arm SD of the change score on the normalised scale
##   "se" - within-arm SE of the mean change (default; LS-mean SE for
##          Mitchell 2023 because no SD is reported)
##   "ci" - 95% CI of the within-arm mean (= SE * 1.96, large-n approx)
## Mitchell 2023 has no SD reported, so the "sd" variant leaves that
## study's whiskers blank (the same way Feder trials are skipped here
## because they don't report within-group change at all).

build_panelA <- function(error_type = c("se", "sd", "ci")) {

  error_type <- match.arg(error_type)

  dat <- caps |>
    filter(!is.na(ch_active_norm)) |>
    mutate(
      act_val    = ch_active_norm,
      pla_val    = ch_placebo_norm,
      act_sd     = ch_active_sd  * norm_factor,
      pla_sd     = ch_placebo_sd * norm_factor,
      act_se     = ch_active_se_n,
      pla_se     = ch_placebo_se_n,
      caps_annot = caps_version
    ) |>
    pivot_longer(
      cols      = c(act_val, pla_val),
      names_to  = "arm",
      values_to = "change"
    ) |>
    filter(!is.na(change)) |>
    mutate(
      arm_label = ifelse(arm == "act_val", "Active", "Control"),
      arm_label = factor(arm_label, levels = c("Active", "Control")),
      err_sd    = ifelse(arm == "act_val", act_sd, pla_sd),
      err_se    = ifelse(arm == "act_val", act_se, pla_se),
      err_val   = switch(error_type,
                         "sd" = err_sd,
                         "se" = err_se,
                         "ci" = err_se * 1.96)
    )

  subtitle_text <- switch(
    error_type,
    "sd" = "Bars: mean change.  Error bars: +/-1 within-arm SD (where SD reported).",
    "se" = paste("Bars: mean change.  Error bars: +/-1 SEM within arm",
                 "(SD/sqrt(n) where SD reported; LS-mean SE from MMRM CI for Mitchell 2023)."),
    "ci" = "Bars: mean change.  Error bars: 95% CI of the within-arm mean."
  )

  ggplot(dat, aes(x = study_label, y = change,
                  fill   = interaction(arm_label, drug),
                  colour = interaction(arm_label, drug))) +
    geom_col(position = position_dodge(width = 0.72),
             width = 0.66, alpha = 0.88) +
    geom_errorbar(aes(ymin = change - err_val,
                      ymax = change + err_val),
                  position  = position_dodge(width = 0.72),
                  width = 0.22, linewidth = 0.4,
                  colour = "grey25", na.rm = TRUE) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    geom_text(
      data = dat |> filter(arm == "act_val"),
      aes(y = change / 2, label = sig_label),
      hjust = 0.5, size = 2.7, colour = "white",
      position = position_dodge(width = 0.72), fontface = "bold"
    ) +
    geom_text(
      data = dat |> filter(arm == "act_val") |>
             distinct(study_label, drug, caps_annot),
      aes(y = 1.5, label = caps_annot, x = study_label,
          group = study_label),
      hjust  = 0.5, vjust = 0, size = 1.9, colour = "grey55",
      inherit.aes = FALSE
    ) +
    facet_grid(. ~ drug, scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = c(
        "Active.MDMA"      = COL_ACTIVE_MDMA,
        "Control.MDMA"     = COL_PLACEBO_MDMA,
        "Active.Ketamine"  = COL_ACTIVE_KET,
        "Control.Ketamine" = COL_PLACEBO_KET),
      labels = c("MDMA (active)", "MDMA (control)",
                 "Ketamine (active)", "Ketamine (control)"),
      name  = "Drug & arm"
    ) +
    scale_colour_manual(
      values = c(
        "Active.MDMA"      = COL_ACTIVE_MDMA,
        "Control.MDMA"     = COL_PLACEBO_MDMA,
        "Active.Ketamine"  = COL_ACTIVE_KET,
        "Control.Ketamine" = COL_PLACEBO_KET),
      guide = "none"
    ) +
    scale_y_continuous(
      breaks = seq(-40, 5, 5),
      labels = function(x) as.character(x),
      expand = expansion(mult = c(0.05, 0.10))
    ) +
    labs(
      x        = NULL,
      y        = "Change in CAPS score\n(CAPS-5 normalised; negative = improvement)",
      title    = "A   Mean CAPS Change from Baseline (Active vs Control)",
      subtitle = subtitle_text
    ) +
    theme_pub() +
    theme(legend.position = "right")
}

## ----------------------------------------------------------------
## 4.  Panel B - forest plot of between-group difference (95% CI)
## ----------------------------------------------------------------

build_panelB <- function() {

  dat <- caps |>
    filter(!is.na(between_norm)) |>
    mutate(
      diff_val = between_norm,
      ci_lo    = ci_lo_norm,
      ci_hi    = ci_hi_norm,
      d_label  = ifelse(!is.na(cohen_d),
                        sprintf("d=%.2f", cohen_d), ""),
      ctrl_lbl = case_when(
        sig_label == "***" ~ paste("***", placebo_type),
        sig_label == "**"  ~ paste("**",  placebo_type),
        sig_label == "ns"  ~ paste("ns ", placebo_type),
        TRUE               ~ placebo_type
      )
    )

  # Plot order: top-to-bottom matches the data row order
  dat$study_label <- factor(dat$study_label,
                            levels = rev(levels(dat$study_label)))

  x_min <- min(dat$ci_lo, dat$diff_val, na.rm = TRUE) - 4
  x_max <- max(dat$ci_hi, dat$diff_val, na.rm = TRUE) + 4
  right_text_x <- x_max + 6

  ggplot(dat, aes(x = diff_val, y = study_label, colour = drug)) +
    geom_vline(xintercept = 0, colour = "grey40",
               linewidth = 0.5, linetype = "dashed") +
    geom_errorbar(
      aes(xmin = ci_lo, xmax = ci_hi),
      orientation = "y", width = 0.3,
      linewidth = 0.65, na.rm = TRUE
    ) +
    geom_point(aes(size = n_active), shape = 18) +
    geom_text(
      aes(x = diff_val, y = study_label, label = d_label),
      vjust = -1.3, hjust = 0.5, size = 2.6, colour = "grey20",
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(x = right_text_x, label = ctrl_lbl),
      hjust = 0, size = 2.3, colour = "grey35"
    ) +
    facet_grid(drug ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(
      values = c(MDMA = COL_ACTIVE_MDMA, Ketamine = COL_ACTIVE_KET),
      guide  = "none"
    ) +
    scale_size_continuous(
      name   = "n (active arm)",
      range  = c(2.5, 6),
      breaks = c(7, 12, 15, 22, 46, 53)
    ) +
    scale_x_continuous(
      limits = c(x_min, right_text_x + 30),
      breaks = seq(-80, 30, 10),
      labels = function(x) as.character(x)
    ) +
    scale_y_discrete() +
    labs(
      x        = paste("Between-group CAPS difference, active - control",
                       "\n(CAPS-5 normalised; negative = active arm improved more)"),
      y        = NULL,
      title    = "B   Between-Group Difference (Forest Plot)",
      subtitle = "Point: LS-mean diff (MMRM/mixed-effects ITT). Error bars: 95% CI. Diamond size proportional to n."
    ) +
    theme_pub() +
    theme(strip.text.y    = element_text(angle = 0),
          legend.position = "right")
}

## ----------------------------------------------------------------
## 5.  Panel C - response rate with proportion SE error bars
## ----------------------------------------------------------------
##
## Error bars represent +/- 1 standard error on the proportion, derived
## from the Wilson score interval:
##   wilson_se(p, n) = (CI_high - CI_low) / (2 * 1.96)
## This is conceptually equivalent to +/- 1 SE under the Wilson framework.
## It is narrower than the corresponding 95% Wilson CI by a factor of
## 1.96, gives non-zero error bars even at boundary proportions (p = 0
## or p = 1), and is visually consistent with the +/- 1 SEM bars in
## panel A. Statistical significance (asterisks) comes from the 2x2
## chi-squared test on the underlying counts, not from the SE bars.

wilson_ci <- function(p_pct, n) {
  if (is.na(p_pct) || is.na(n)) return(list(lo = NA, hi = NA))
  p <- p_pct / 100
  z <- qnorm(0.975)
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lo = max(0, centre - half) * 100,
       hi = min(1, centre + half) * 100)
}

wilson_se <- function(p_pct, n) {
  ci <- wilson_ci(p_pct, n)
  if (is.na(ci$lo) || is.na(ci$hi)) return(NA_real_)
  (ci$hi - ci$lo) / (2 * 1.96)
}

build_panelC <- function(error_type = c("se", "sd", "ci")) {

  error_type <- match.arg(error_type)

  dat <- caps |>
    filter(!is.na(resp_active)) |>
    select(study_label, drug, sig_resp,
           resp_active, resp_placebo, n_active, n_placebo)

  # SD of a single binary observation is sqrt(p(1-p)) in percent units.
  prop_sd <- function(p_pct) {
    if (is.na(p_pct)) return(NA_real_)
    p <- p_pct / 100
    100 * sqrt(p * (1 - p))
  }

  dat_long <- dat |>
    pivot_longer(
      cols      = c(resp_active, resp_placebo),
      names_to  = "arm",
      values_to = "pct"
    ) |>
    mutate(
      arm_n     = ifelse(arm == "resp_active", n_active, n_placebo),
      arm_label = ifelse(arm == "resp_active", "Active", "Control"),
      arm_label = factor(arm_label, levels = c("Active", "Control")),
      se_p      = mapply(wilson_se, pct, arm_n),
      sd_p      = sapply(pct, prop_sd),
      ci_obj    = mapply(wilson_ci, pct, arm_n, SIMPLIFY = FALSE),
      ci95_lo   = sapply(ci_obj, function(x) x$lo),
      ci95_hi   = sapply(ci_obj, function(x) x$hi)
    ) |>
    mutate(
      ci_lo = switch(error_type,
                     "sd" = pmax(0,   pct - sd_p),
                     "se" = pmax(0,   pct - se_p),
                     "ci" = ci95_lo),
      ci_hi = switch(error_type,
                     "sd" = pmin(100, pct + sd_p),
                     "se" = pmin(100, pct + se_p),
                     "ci" = ci95_hi)
    )

  subtitle_text <- switch(
    error_type,
    "sd" = "Response: >=10pt CAPS-5 / >=30% CAPS-IV reduction. Error bars: SD (binary).",
    "se" = "Response: >=10pt CAPS-5 / >=30% CAPS-IV reduction.",
    "ci" = "Response: >=10pt CAPS-5 / >=30% CAPS-IV reduction. Error bars: 95% CI (Wilson)."
  )

  ggplot(dat_long, aes(x = study_label, y = pct,
                       fill   = interaction(arm_label, drug),
                       colour = interaction(arm_label, drug))) +
    geom_col(position = position_dodge(width = 0.72),
             width = 0.66, alpha = 0.88) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  position = position_dodge(width = 0.72),
                  width = 0.22, linewidth = 0.4, colour = "grey25") +
    geom_text(
      aes(y = pmin(ci_hi + 3, 113),
          label = sprintf("%.0f%%", pct)),
      position = position_dodge(width = 0.72),
      vjust = 0, size = 2.3, colour = "grey20"
    ) +
    geom_text(
      data = dat_long |> filter(arm == "resp_active"),
      aes(y = 5, label = sig_resp),
      size = 2.5, colour = "white", fontface = "bold",
      position = position_dodge(width = 0.72)
    ) +
    facet_grid(. ~ drug, scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = c(
        "Active.MDMA"      = COL_ACTIVE_MDMA,
        "Control.MDMA"     = COL_PLACEBO_MDMA,
        "Active.Ketamine"  = COL_ACTIVE_KET,
        "Control.Ketamine" = COL_PLACEBO_KET),
      labels = c("MDMA (active)", "MDMA (control)",
                 "Ketamine (active)", "Ketamine (control)"),
      name   = "Drug & arm"
    ) +
    scale_colour_manual(
      values = c(
        "Active.MDMA"      = COL_ACTIVE_MDMA,
        "Control.MDMA"     = COL_PLACEBO_MDMA,
        "Active.Ketamine"  = COL_ACTIVE_KET,
        "Control.Ketamine" = COL_PLACEBO_KET),
      guide = "none"
    ) +
    scale_y_continuous(
      limits = c(0, 120),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x        = NULL,
      y        = "Participants (%)",
      title    = "C   Response Rate",
      subtitle = subtitle_text
    ) +
    theme_pub() +
    theme(legend.position = "none")
}

## ----------------------------------------------------------------
## 6.  Panel D - loss of PTSD diagnosis (% with binomial 95% CIs)
## ----------------------------------------------------------------

build_panelD <- function(error_type = c("se", "sd", "ci")) {

  error_type <- match.arg(error_type)

  dat <- caps |>
    filter(!is.na(nodiag_active)) |>
    select(study_label, drug, sig_diag,
           nodiag_active, nodiag_placebo, n_active, n_placebo)

  prop_sd <- function(p_pct) {
    if (is.na(p_pct)) return(NA_real_)
    p <- p_pct / 100
    100 * sqrt(p * (1 - p))
  }

  dat_long <- dat |>
    pivot_longer(
      cols      = c(nodiag_active, nodiag_placebo),
      names_to  = "arm",
      values_to = "pct"
    ) |>
    mutate(
      arm_n     = ifelse(arm == "nodiag_active", n_active, n_placebo),
      arm_label = ifelse(arm == "nodiag_active", "Active", "Control"),
      arm_label = factor(arm_label, levels = c("Active", "Control")),
      se_p      = mapply(wilson_se, pct, arm_n),
      sd_p      = sapply(pct, prop_sd),
      ci_obj    = mapply(wilson_ci, pct, arm_n, SIMPLIFY = FALSE),
      ci95_lo   = sapply(ci_obj, function(x) x$lo),
      ci95_hi   = sapply(ci_obj, function(x) x$hi)
    ) |>
    mutate(
      ci_lo = switch(error_type,
                     "sd" = pmax(0,   pct - sd_p),
                     "se" = pmax(0,   pct - se_p),
                     "ci" = ci95_lo),
      ci_hi = switch(error_type,
                     "sd" = pmin(100, pct + sd_p),
                     "se" = pmin(100, pct + se_p),
                     "ci" = ci95_hi)
    )

  subtitle_text <- switch(
    error_type,
    "sd" = "Participants no longer meeting DSM PTSD criteria. Error bars: SD (binary).",
    "se" = "Participants no longer meeting DSM PTSD criteria.",
    "ci" = "Participants no longer meeting DSM PTSD criteria. Error bars: 95% CI (Wilson)."
  )

  ggplot(dat_long, aes(x = study_label, y = pct,
                       fill   = interaction(arm_label, drug),
                       colour = interaction(arm_label, drug))) +
    geom_col(position = position_dodge(width = 0.72),
             width = 0.66, alpha = 0.88) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  position = position_dodge(width = 0.72),
                  width = 0.22, linewidth = 0.4, colour = "grey25") +
    geom_text(
      aes(y = pmin(ci_hi + 3, 113),
          label = sprintf("%.0f%%", pct)),
      position = position_dodge(width = 0.72),
      vjust = 0, size = 2.3, colour = "grey20"
    ) +
    geom_text(
      data = dat_long |> filter(arm == "nodiag_active"),
      aes(y = 5, label = sig_diag),
      size = 2.5, colour = "white", fontface = "bold",
      position = position_dodge(width = 0.72)
    ) +
    facet_grid(. ~ drug, scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = c(
        "Active.MDMA"      = COL_ACTIVE_MDMA,
        "Control.MDMA"     = COL_PLACEBO_MDMA,
        "Active.Ketamine"  = COL_ACTIVE_KET,
        "Control.Ketamine" = COL_PLACEBO_KET),
      labels = c("MDMA (active)", "MDMA (control)",
                 "Ketamine (active)", "Ketamine (control)"),
      name   = "Drug & arm"
    ) +
    scale_colour_manual(
      values = c(
        "Active.MDMA"      = COL_ACTIVE_MDMA,
        "Control.MDMA"     = COL_PLACEBO_MDMA,
        "Active.Ketamine"  = COL_ACTIVE_KET,
        "Control.Ketamine" = COL_PLACEBO_KET),
      guide = "none"
    ) +
    scale_y_continuous(
      limits = c(0, 120),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x        = NULL,
      y        = "Participants (%)",
      title    = "D   Loss of PTSD Diagnosis",
      subtitle = subtitle_text
    ) +
    theme_pub() +
    theme(legend.position = "none")
}

## ----------------------------------------------------------------
## 7.  Assemble figure - three error-bar variants
## ----------------------------------------------------------------
cat("\nBuilding figures (3 error-bar variants)...\n")

# Forest plot (Panel B) is always 95% CI - the alternatives are not
# sensible for a between-group difference (see methods).
pB <- build_panelB()

variants <- list(
  list(suffix = "SE",  label = "+/- 1 standard error of the mean"),
  list(suffix = "SD",  label = "+/- 1 standard deviation (within-arm)"),
  list(suffix = "CI",  label = "95% confidence interval of the mean")
)

for (v in variants) {

  et <- tolower(v$suffix)
  pA <- build_panelA(error_type = et)
  pC <- build_panelC(error_type = et)
  pD <- build_panelD(error_type = et)

  bottom_row <- (pC | pD) + plot_layout(widths = c(1.4, 1))

  fig_main <- (pA / pB / bottom_row) +
    plot_annotation(
      title    = "CAPS-Based PTSD Outcomes: MDMA-AT and Ketamine RCTs",
      subtitle = paste0(
        "CAPS-IV normalised to CAPS-5 scale (x80/136).  ",
        "Significance: *** P<0.001  ** P<0.01  * P<0.05  ns: not significant.\n",
        "Error bars in panels A, C, D: ", v$label, "."),
      caption  = paste0(
        "Sources: Mitchell 2021 (Nat Med); Mitchell 2023 (Nat Med); ",
        "Mithoefer 2018 (Lancet Psychiatry); Feder 2014 (JAMA Psychiatry); ",
        "Feder 2021 (Am J Psychiatry)."),
      theme = theme(
        plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5, colour = "grey25",
                                     lineheight = 1.2),
        plot.caption  = element_text(size = 6.5, colour = "grey35",
                                     hjust = 0, lineheight = 1.3))
    )

  fname_png <- sprintf("Figure_CAPS_Outcomes_%s.png", v$suffix)
  fname_pdf <- sprintf("Figure_CAPS_Outcomes_%s.pdf", v$suffix)

  png(fname_png, width = 220, height = 260, units = "mm", res = 300)
  print(fig_main)
  invisible(dev.off())
  cat("Saved:", fname_png, "\n")

  pdf(fname_pdf, width = 220 / 25.4, height = 260 / 25.4)
  print(fig_main)
  invisible(dev.off())
  cat("Saved:", fname_pdf, "\n")
}

## ----------------------------------------------------------------
## 8.  Demographic table (gridExtra tableGrob) - PNG + PDF
## ----------------------------------------------------------------

demo_tbl_data <- demo_raw |>
  mutate(across(everything(), ~ ifelse(is.na(.x), "-", as.character(.x))))

# Shorten the longer Mithoefer arm label so the column doesn't blow out
demo_tbl_data$Arm <- gsub("30 mg \\(active control\\)", "30 mg (ctrl)",
                         demo_tbl_data$Arm)

nice_cols <- c("Study","Arm","n","Age (y)","Female","White","Hisp/Latino",
               "BMI","PTSD dur. (y)","Dissoc.","MDD",
               "Veteran","Baseline CAPS")
names(demo_tbl_data) <- nice_cols

demo_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = 8, hjust = 0, x = 0.04),
    bg_params = list(fill = rep(c("white", "grey97"),
                                length.out = nrow(demo_tbl_data)))
  ),
  colhead = list(
    fg_params = list(fontsize = 8.5, fontface = "bold",
                     hjust = 0, x = 0.04),
    bg_params = list(fill = "grey92", col = NA)
  ),
  padding = unit(c(2.5, 2.5), "mm")
)

demo_grob <- tableGrob(demo_tbl_data, rows = NULL, theme = demo_theme)

# Wrap in a ggplot via patchwork so we can use plot_annotation
# for clean title/subtitle/caption.
tbl_plot <- patchwork::wrap_elements(full = demo_grob) +
  plot_annotation(
    title    = "Table 1.  Demographic and clinical characteristics at baseline",
    subtitle = "MDMA-AT and ketamine RCTs included in the CAPS-based outcomes synthesis",
    caption  = paste0(
      "Values are mean (SD) or n (%) as reported in the source paper. ",
      "NR = not reported.  Dissoc. = dissociative-subtype PTSD; ",
      "MDD = comorbid major depressive disorder.\n",
      "Mithoefer 2018 demographics are arm-stratified (whole cohort: 22 veterans, ",
      "3 firefighters, 1 police officer).\n",
      "Feder 2014 reports first-period demographics for the crossover design."
    ),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, colour = "grey25", hjust = 0,
                                   margin = margin(b = 4)),
      plot.caption  = element_text(size = 7, colour = "grey35",
                                   hjust = 0, face = "italic",
                                   margin = margin(t = 4),
                                   lineheight = 1.2),
      plot.margin   = margin(6, 6, 6, 6)
    )
  )

# Compute height from the actual table content (so we don't get whitespace).
# Table height (mm) ~ rows * 7 + header 8 + margins
tbl_h_mm <- 7 * (nrow(demo_tbl_data) + 1) + 38   # rows + header rows + chrome
tbl_w_mm <- 330

png("Table_Demographics.png", width = tbl_w_mm, height = tbl_h_mm,
    units = "mm", res = 300)
print(tbl_plot)
invisible(dev.off())
cat("Saved: Table_Demographics.png\n")

pdf("Table_Demographics.pdf", width = tbl_w_mm / 25.4, height = tbl_h_mm / 25.4)
print(tbl_plot)
invisible(dev.off())
cat("Saved: Table_Demographics.pdf\n")

## ----------------------------------------------------------------
## 9.  Summary-statistics table (efficacy)
## ----------------------------------------------------------------
summary_tbl <- caps |>
  select(label, drug, caps_version, placebo_type,
         n_active, n_placebo, bl_active, bl_placebo,
         between_diff, ci_lo_raw, ci_hi_raw, cohen_d,
         resp_active, resp_placebo, p_resp,
         nodiag_active, nodiag_placebo, p_diag, p_val) |>
  mutate(label = gsub("\n", " ", label))

write.csv(summary_tbl, "summary_statistics.csv", row.names = FALSE)
cat("Saved: summary_statistics.csv\n")

cat("\nAll outputs written to working directory.\n")
