library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)

mo2lo <- read.table(snakemake@input[["ratios_table"]], header = TRUE, sep = "\t")

confint_upper <- mo2lo$gnomad_confint_upper
confint_lower <- mo2lo$gnomad_confint_lower

by_gnomad <- mo2lo %>%
  select(AA, mo2lo_gnomad) %>%
  arrange(mo2lo_gnomad) %>%
  select(AA) %>%
  unlist()
mo2lo$AA <- factor(mo2lo$AA, levels = by_gnomad)
mo2lo %>%
  ggplot(aes(
    x = factor(AA), y = mo2lo_gnomad
  )) +
  ylab("Percentage of\noptimality-reducing variants") +
  xlab("Amino acids") +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = confint_lower, ymax = confint_upper),
    width = 0.2
  ) +
  geom_vline(linetype = "dotted", xintercept = seq(1.5, 8.5, by = 1.0)) +
  theme_classic() +
  theme(
    text = element_text(size = 18),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    panel.spacing = unit(3, "lines")
  ) -> mo2lo_gnomad

ggsave(mo2lo_gnomad, width = 9, height = 6, filename = snakemake@output[["plot"]])
