# ==============================================================================
# Plotting helpers shared across analysis scripts
# Consistent colour scheme: ACPL purple, MST red, neutral grey for ties.
# ==============================================================================

#' ACPL colour palette
#' @export
acpl_colours <- list(
  acpl     = "#7b5ea7",
  acpl_dk  = "#512DA8",
  mst      = "#cc6655",
  sling    = "#4caf50",
  paga     = "#2a7ab5",
  neutral  = "#aaaaaa",
  G1       = "#D1C4E9",
  S        = "#9575CD",
  G2M      = "#512DA8"
)

#' Ridge plot of phase distributions along pseudotime
#'
#' @param df Data frame with columns: pseudotime (numeric), Phase (G1/S/G2M)
#' @param method_name Title for the plot
#' @param sw_sa SW-SA value to display in subtitle
#' @return ggplot object
#' @export
phase_ridge <- function(df, method_name, sw_sa) {
  if (!requireNamespace("ggridges", quietly = TRUE))
    stop("Install ggridges: install.packages('ggridges')")
  ggplot2::ggplot(
      df[!is.na(df$Phase), ],
      ggplot2::aes_string(x = "pseudotime", y = "Phase", fill = "Phase")) +
    ggridges::geom_density_ridges(
      scale = 1.4, alpha = 0.85,
      colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c(
      G1  = acpl_colours$G1,
      S   = acpl_colours$S,
      G2M = acpl_colours$G2M)) +
    ggridges::theme_ridges() +
    ggplot2::labs(
      title    = method_name,
      subtitle = sprintf("SW-SA: %.1f%%", sw_sa),
      x        = "Pseudotime") +
    ggplot2::theme(legend.position = "none",
                   plot.title    = ggplot2::element_text(face = "bold"))
}

#' Forest plot of bootstrap CIs
#' @export
bootstrap_forest <- function(boot_df) {
  pd <- boot_df
  pd$Label <- paste0(pd$Dataset, "\n(", pd$Method_A, " vs ", pd$Method_B, ")")
  ggplot2::ggplot(pd,
      ggplot2::aes(x = Obs_pp,
                   y = stats::reorder(Label, Obs_pp),
                   colour = Sig)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                         colour = "grey50", linewidth = 0.8) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = CI_lo_bca, xmax = CI_hi_bca),
      height = 0.22, linewidth = 0.9) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = acpl_colours$acpl_dk,
                 "FALSE" = acpl_colours$mst),
      labels = c("TRUE" = "CI excludes 0",
                 "FALSE" = "CI includes 0"),
      name = NULL) +
    ggplot2::labs(
      title    = "Pairwise Bootstrap CIs (95% BCa)",
      subtitle = "10,000 resamples per comparison",
      x        = "SW-SA difference (pp, Method A minus B)",
      y        = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "top")
}
