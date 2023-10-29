# Generates a figure with deltas of "adjusted proportion of singletons" scores

library(dplyr)
library(ggplot2)
library(ggprism)
library(latex2exp)

scores <- read.table(snakemake@input[["scores"]],
  header = TRUE,
  sep = "\t"
)

if (!is.null(snakemake@params[["xlab_labels_set"]])) scores <- filter(scores, variable_value %in% snakemake@params[["xlab_labels_set"]])

if (!is.null(snakemake@params[["xlab_labels"]]) || !is.null(snakemake@params[["new_xlab_labels"]])) {
  if (!is.null(snakemake@params[["xlab_labels"]]) && !is.null(snakemake@params[["new_xlab_labels"]])) {
    scores[["variable_value"]] <- factor(scores[["variable_value"]],
      levels = snakemake@params[["xlab_labels"]],
      labels = snakemake@params[["new_xlab_labels"]]
    )
  } else {
    stop("Either original or new labels were not provided")
  }
}

# TODO: make this filtering specific: scores[!is.na(...),]
if (!is.null(snakemake@params[["NA_omit"]]) && snakemake@params[["NA_omit"]]) {
  scores <- na.omit(scores)
}

scores <- scores %>%
  mutate(
    diff = !!sym(snakemake@params[["score_name"]]) - scores[scores$variable_value == snakemake@params[["ref_group"]], ][[snakemake@params[["score_name"]]]],
    diff_sem = sqrt(
      (!!sym(paste0(snakemake@params[["score_name"]], "_sem")))^2 + (scores[scores$variable_value == snakemake@params[["ref_group"]], ][[paste0(snakemake@params[["score_name"]], "_sem")]])^2
    )
  ) %>%
  mutate(
    diff_uconf = diff + 1.96 * diff_sem,
    diff_lconf = diff - 1.96 * diff_sem
  )

scores[[snakemake@params[["score_name"]]]] <- scores$diff
scores[[paste0(snakemake@params[["score_name"]], "_sem")]] <- scores$diff_sem
scores[[paste0(snakemake@params[["score_name"]], "_uconf")]] <- scores$diff_uconf
scores[[paste0(snakemake@params[["score_name"]], "_lconf")]] <- scores$diff_lconf

pdf(snakemake@output[["plot"]])
ggplot(scores) +
  {
    if (snakemake@params[["reorder_xlab_by_score"]]) {
      aes(x = reorder(factor(variable_value), !!sym(snakemake@params[["score_name"]])), y = !!sym(snakemake@params[["score_name"]]))
    } else {
      aes(x = factor(variable_value), y = !!sym(snakemake@params[["score_name"]]))
    }
  } +
  ylab(ifelse(snakemake@params[["score_name"]] %in% c("maps", "traps", "caps", "caps_pdd"), case_when(
    (snakemake@params[["score_name"]] == "maps") ~ "MAPS",
    (snakemake@params[["score_name"]] == "traps") ~ "TRAPS",
    (snakemake@params[["score_name"]] == "caps") ~ "CAPS",
    (snakemake@params[["score_name"]] == "caps_pdd") ~ "CAPS-PDD"
  ), stop("Score name error"))) +
  {
    if (!is.null(snakemake@params[["ylab"]])) {
      ylab(
        TeX(snakemake@params[["ylab"]])
      )
    }
  } +
  xlab(snakemake@params[["xlab"]]) +
  geom_pointrange(
    aes(
      ymin = !!sym(paste0(snakemake@params[["score_name"]], "_lconf")),
      ymax = !!sym(paste0(snakemake@params[["score_name"]], "_uconf"))
    ),
    size = 1.3, linewidth = 1.8,
  ) +
  {
    if (!is.null(snakemake@params[["ylim_min"]]) &&
      !is.null(snakemake@params[["ylim_max"]])) {
      ylim(
        snakemake@params[["ylim_min"]],
        snakemake@params[["ylim_max"]]
      )
    }
  } +
  theme_classic() +
  theme(
    aspect.ratio = ifelse(!is.null(snakemake@params[["aspect_ratio"]]), snakemake@params[["aspect_ratio"]], 1),
    legend.direction = "horizontal",
    text = element_text(size = 24),
    axis.text.x = element_text(
      size = ifelse(is.null(snakemake@params[["xlab_size"]]), 24, snakemake@params[["xlab_size"]]),
      vjust = snakemake@params[["xlab_vjust"]],
      hjust = snakemake@params[["xlab_hjust"]],
      angle = ifelse(is.null(snakemake@params[["xlab_angle"]]), 0, snakemake@params[["xlab_angle"]])
    ),
    plot.margin = margin(0, 5.5, 0, 5.5)
  ) +
  {
    if (!is.null(snakemake@params[["add_pvalues"]]) && snakemake@params[["add_pvalues"]] == TRUE) {
      add_pvalue(pvals, label.size = 4.5)
    }
  }
dev.off()
