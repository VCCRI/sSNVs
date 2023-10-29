# Analysis of different distributions

library(dplyr)
library(ggplot2)
library(ggprism)
library(magrittr)
library(stringr)

MO2LO <- "Optimality reduced"
LO2MO <- "Optimality increased"

two_codon_aa <- c("N", "K", "F", "C", "Q", "E", "Y", "H", "D")
two_codon_aa_full <- c("Asn", "Lys", "Phe", "Cys", "Gln", "Glu", "Tyr", "His", "Asp")

dist_analysis <- function(gnomad,
                          dist_name = "GERP",
                          dist = "gerp",
                          violin_scale = "count",
                          caption = NULL,
                          legend_title = NULL,
                          AA_full = FALSE,
                          AF_upper = 1,
                          AF_lower = 0,
                          label = NULL,
                          lower_bound = 0,
                          log10_trans = FALSE,
                          plot = "violin",
                          add_pvalues = NULL,
                          pvalue_stars = NULL,
                          colors = NULL,
                          xlab = NULL,
                          plot_name = NULL,
                          only_2codon = TRUE,
                          tests = TRUE,
                          bonferroni = NULL,
                          tests_output = NULL,
                          ks = FALSE) {
  if (is.character(gnomad)) {
    gnomad <- read.table(gnomad, header = TRUE)
  }
  if (only_2codon) {
    gnomad <- gnomad %>%
      filter(AA %in% two_codon_aa) %>%
      mutate(
        Ref = stringr::str_extract(ref_codon, "[A-Z]"),
        Alt = stringr::str_extract(alt_codon, "[A-Z]")
      ) %>%
      mutate(optimality = case_when(
        # For Asn, aaC->aaT creates a less optimal codon
        (AA == "N" & Ref == "C" & Alt == "T") ~ MO2LO,
        (AA == "N" & Ref == "T" & Alt == "C") ~ LO2MO,
        # For Lys, aaG->aaA creates a less optimal codon
        (AA == "K" & Ref == "G" & Alt == "A") ~ MO2LO,
        (AA == "K" & Ref == "A" & Alt == "G") ~ LO2MO,
        # For Phe, ttC->ttT creates a less optimal codon
        (AA == "F" & Ref == "C" & Alt == "T") ~ MO2LO,
        (AA == "F" & Ref == "T" & Alt == "C") ~ LO2MO,
        # For Cys, tgC->tgT creates a less optimal codon
        (AA == "C" & Ref == "C" & Alt == "T") ~ MO2LO,
        (AA == "C" & Ref == "T" & Alt == "C") ~ LO2MO,
        # For Gln, caG->caA creates a less optimal codon
        (AA == "Q" & Ref == "G" & Alt == "A") ~ MO2LO,
        (AA == "Q" & Ref == "A" & Alt == "G") ~ LO2MO,
        # For Glu, gaG->gaA creates a less optimal codon
        (AA == "E" & Ref == "G" & Alt == "A") ~ MO2LO,
        (AA == "E" & Ref == "A" & Alt == "G") ~ LO2MO,
        # For Tyr, taC->taT creates a less optimal codon
        (AA == "Y" & Ref == "C" & Alt == "T") ~ MO2LO,
        (AA == "Y" & Ref == "T" & Alt == "C") ~ LO2MO,
        # For His, caC->caT creates a less optimal codon
        (AA == "H" & Ref == "C" & Alt == "T") ~ MO2LO,
        (AA == "H" & Ref == "T" & Alt == "C") ~ LO2MO,
        # For Asp, gaC->gaT creates a less optimal codon
        (AA == "D" & Ref == "C" & Alt == "T") ~ MO2LO,
        (AA == "D" & Ref == "T" & Alt == "C") ~ LO2MO,
      ))
  }

  if (tests) {
    lapply(two_codon_aa, function(aa) {
      x <- gnomad %>%
        filter(AA == aa & optimality == MO2LO) %>%
        select(!!ensym(dist)) %>%
        unlist()
      y <- gnomad %>%
        filter(AA == aa & optimality == LO2MO) %>%
        select(!!ensym(dist)) %>%
        unlist()
      if (ks) {
        data.frame(
          AA = ifelse(AA_full, two_codon_aa_full[match(aa, two_codon_aa)], aa),
          med.diff = median(x) - median(y),
          ks.stat = suppressWarnings(ks.test(x, y)$statistic),
          p = wilcox.test(x, y, alternative = "greater")$p.value
        )
      } else {
        data.frame(
          AA = ifelse(AA_full, two_codon_aa_full[match(aa, two_codon_aa)], aa),
          med.diff = median(x) - median(y),
          p = wilcox.test(x, y, alternative = "greater")$p.value
        )
      }
    }) %>%
      data.table::rbindlist() %>%
      arrange(desc(med.diff)) -> tests
    if (bonferroni) {
      tests <- mutate(tests, p = p.adjust(p, "bonferroni"))
    }
    if (!is.null(tests_output)) {
      if (ks) {
        tests_out <- tests
        colnames(tests_out) <- c("Amino acid", "Difference in medians", "KS test statistic", "Wilcoxon p-value")
        print(xtable::xtable(tests_out, caption = caption, label = label, display = c("s", "s", "g", "g", "g")), include.rownames = FALSE) %>% write(tests_output)
      } else {
        tests_out <- tests
        colnames(tests_out) <- c("Amino acid", "Difference in medians", "Wilcoxon p-value")
        print(xtable::xtable(tests_out, caption = caption, label = label, display = c("s", "s", "g", "g")), include.rownames = FALSE) %>% write(tests_output)
      }
    }
  }

  if (plot == "AF_dist") {
    gnomad %>%
      filter(AC == 1) %>%
      select(AF) %>%
      max() -> last_singleton_AF
    gnomad %>%
      filter(AF >= AF_lower & AF <= AF_upper) %>%
      ggplot() +
      {
        if (legend_title == "Chromosome X?") aes(x = AF, color = (locus.contig == "X"))
      } +
      {
        if (legend_title == "AC") aes(x = AF, color = as.factor(AC))
      } +
      {
        if (legend_title == "Singleton or doubleton?") aes(x = AF, color = as.factor(AC %in% c(1, 2)))
      } +
      {
        if (legend_title == "Singleton?") aes(x = AF, color = as.factor(AC == 1))
      } +
      {
        if (legend_title == "Singleton?") geom_vline(xintercept = last_singleton_AF)
      } +
      {
        if (legend_title == "Singleton?") annotate(x = last_singleton_AF, y = +Inf, label = "Last singleton", vjust = 2, geom = "label")
      } +
      ylab("count") +
      geom_histogram(bins = 200) +
      theme_minimal() +
      theme(
        text = element_text(size = 24),
        legend.position = "bottom",
        legend.direction = "horizontal"
      ) +
      labs(color = legend_title) -> distribution
    if (is.null(plot_name)) {
      print(distribution)
    } else {
      ggsave(plot = distribution, filename = plot_name, device = "pdf", width = 24, height = 12)
    }
  }

  if (plot == "violin") {
    gnomad %>%
      `if`(AA_full, mutate(., AA = case_when(
        (AA == two_codon_aa[1]) ~ two_codon_aa_full[1],
        (AA == two_codon_aa[2]) ~ two_codon_aa_full[2],
        (AA == two_codon_aa[3]) ~ two_codon_aa_full[3],
        (AA == two_codon_aa[4]) ~ two_codon_aa_full[4],
        (AA == two_codon_aa[5]) ~ two_codon_aa_full[5],
        (AA == two_codon_aa[6]) ~ two_codon_aa_full[6],
        (AA == two_codon_aa[7]) ~ two_codon_aa_full[7],
        (AA == two_codon_aa[8]) ~ two_codon_aa_full[8],
        (AA == two_codon_aa[9]) ~ two_codon_aa_full[9],
      )), .) %>%
      ggplot(aes(x = AA, y = !!ensym(dist), color = optimality)) +
      ylab(dist_name) +
      geom_violin(scale = violin_scale, size = 1.5, fill = "gray") +
      geom_boxplot(size = 1.2, outlier.shape = NA, fill = "white", width = 0.15, position = position_dodge(width = 0.9)) +
      geom_hline(yintercept = 0, color = "gray", size = 2, alpha = 0.6) +
      facet_grid(. ~ AA, scales = "free", space = "free") +
      theme_minimal() +
      {
        if (log10_trans == TRUE) scale_y_continuous(trans = "log10")
      } +
        scale_color_manual(values = colors) +
      theme(
        text = element_text(size = 48),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = unit(1, "in")
      ) -> violin
    if (add_pvalues && !is.null(pvalue_stars) && pvalue_stars) {
        if (any(tests$p >= 0.0001)) stop("Some p-values above 0.0001 (****)")
        violin <- violin + annotate("text", x = 1, y = 7, label = "****", size = 18)
    }

    if (is.null(plot_name)) {
      print(violin)
    } else {
      ggsave(plot = violin, filename = plot_name, device = "pdf", width = 24, height = 12)
    }
  }

  if (plot == "histogram") {
    gnomad %>%
      ggplot(aes(color = optimality, x = !!ensym(dist))) +
      geom_histogram(bins = 200) +
      {
        if (log10_trans == TRUE) scale_x_continuous(trans = "log10")
      } +
      theme_minimal() +
      theme(
        text = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = "bottom"
      ) -> histogram
    if (is.null(plot_name)) {
      print(histogram)
    } else {
      ggsave(plot = histogram, filename = plot_name, device = "pdf", width = 24, height = 12)
    }
  }

  if (plot == "proportion") {
    lapply(two_codon_aa, function(x) {
      gnomad %>%
        filter(AA == x) %>%
        ggplot(aes(x = !!ensym(dist), color = optimality)) +
        stat_bin(data = subset(gnomad, AA == x & optimality == MO2LO), aes(y = sqrt((..count..) / sum(..count..))), geom = "step", bins = 300) +
        stat_bin(data = subset(gnomad, AA == x & optimality == LO2MO), aes(y = sqrt((..count..) / sum(..count..))), geom = "step", bins = 300) +
        {
          if (log10_trans == TRUE) scale_x_continuous(trans = "log10")
        } +
        theme_minimal() +
        ylab(dist_name) +
        xlab(xlab) +
        ggtitle(x) +
        scale_color_manual(values = colors) +
        scale_y_continuous(labels = scales::percent, limits = c(lower_bound, 1)) +
        theme(
          text = element_text(size = 30),
          legend.title = element_blank(),
          legend.position = "bottom"
        )
    }) -> plots
    ggsave(plot = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE), filename = plot_name, device = "pdf", width = 24, height = 12)
  }

  if (plot == "cdf_all") {
    gnomad %>%
      ggplot(aes(x = !!ensym(dist))) +
      stat_ecdf(pad = FALSE) +
      theme_minimal() +
      ylab(dist_name) +
      xlab(xlab) +
      {
        if (log10_trans == TRUE) scale_x_continuous(trans = "log10")
      } +
      scale_color_manual(values = colors) +
      scale_y_continuous(labels = scales::percent, limits = c(lower_bound, 1)) +
      theme(
        text = element_text(size = 30),
        legend.title = element_blank(),
        legend.position = "bottom"
      ) -> cdf
    if (is.null(plot_name)) {
      print(cdf)
    } else {
      ggsave(plot = cdf, filename = plot_name, device = "pdf", width = 24, height = 12)
    }
  }

  if (plot == "cdf_by_AA") {
    lapply(two_codon_aa, function(x) {
      gnomad %>%
        filter(AA == x) %>%
        ggplot(aes(x = !!ensym(dist), color = optimality)) +
        stat_ecdf(pad = FALSE, size = 2) +
        theme_minimal() +
        ylab(dist_name) +
        xlab(xlab) +
        ggtitle(x) +
        {
          if (log10_trans == TRUE) scale_x_continuous(trans = "log10")
        } +
        scale_color_manual(values = colors) +
        scale_y_continuous(labels = scales::percent, limits = c(lower_bound, 1)) +
        theme(
          text = element_text(size = 40),
          legend.title = element_blank(),
          legend.key = element_blank()
        )
    }) -> plots

    ggsave(plot = ggpubr::ggarrange(plots, common.legend = TRUE), filename = plot_name, device = "pdf", width = 40, height = 20)
    ggsave(plot = ggpubr::ggarrange(plots[[1]], plots[[5]], common.legend = TRUE), filename = paste0("mini_", plot_name), device = "pdf", width = 28, height = 10)
  }

  if (plot == "cdf_by_GC") {
    gnomad <- mutate(gnomad, GC = (str_detect(alt_codon, "..C") | str_detect(alt_codon, "..G")))
    lapply(c(TRUE, FALSE), function(x) {
      gnomad %>%
        filter(GC == x) %>%
        ggplot(aes(x = !!ensym(dist), color = optimality)) +
        stat_ecdf(pad = FALSE) +
        theme_minimal() +
        ylab(dist_name) +
        xlab(xlab) +
        ggtitle(x) +
        {
          if (log10_trans == TRUE) scale_x_continuous(trans = "log10")
        } +
        scale_color_manual(values = colors) +
        scale_y_continuous(labels = scales::percent, limits = c(lower_bound, 1)) +
        theme(
          text = element_text(size = 30),
          legend.title = element_blank(),
          legend.position = "bottom"
        )
    }) -> plots

    ggsave(plot = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE), filename = plot_name, device = "pdf", width = 24, height = 12)
  }

  if (plot == "cumsum") {
    lapply(two_codon_aa, function(x) {
      gnomad %>%
        filter(AA == x) %>%
        ggplot(aes(x = !!ensym(dist), color = optimality)) +
        stat_bin(data = subset(gnomad, AA == x & optimality == MO2LO), aes(y = cumsum(..count..) / sum(..count..)), geom = "step", bins = 100) +
        stat_bin(data = subset(gnomad, AA == x & optimality == LO2MO), aes(y = cumsum(..count..) / sum(..count..)), geom = "step", bins = 100) +
        theme_minimal() +
        ylab(dist_name) +
        xlab(xlab) +
        ggtitle(x) +
        {
          if (log10_trans == TRUE) scale_x_continuous(trans = "log10")
        } +
        scale_color_manual(values = colors) +
        scale_y_continuous(labels = scales::percent, limits = c(lower_bound, 1)) +
        theme(
          text = element_text(size = 30),
          legend.title = element_blank(),
          legend.position = "bottom"
        )
    }) -> plots
    ggsave(plot = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE), filename = plot_name, device = "pdf", width = 24, height = 12)
  }
}
dist_analysis(snakemake@input[["variants"]],
  plot_name = snakemake@output[["plot"]],
  tests_output = snakemake@output[["tests"]],
  bonferroni = snakemake@params[["bonferroni"]],
  colors = snakemake@params[["colors"]],
  add_pvalues = snakemake@params[["add_pvalues"]],
  pvalue_stars = snakemake@params[["pvalue_stars"]]
)
