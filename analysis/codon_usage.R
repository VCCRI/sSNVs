library(magrittr)
library(dplyr)
library(ggplot2)
library(vroom)
library(epitools)

cu_human <- vroom(snakemake@input[["codon_usage_db"]]) %>% filter(
  Species == "Homo sapiens",
  Organelle == "genomic"
)

opt <- c("MO", "LO")
asn <- c("AAC", "AAT")
names(asn) <- opt
lys <- c("AAG", "AAA")
names(lys) <- opt
phe <- c("TTC", "TTT")
names(phe) <- opt
cys <- c("TGC", "TGT")
names(cys) <- opt
gln <- c("CAG", "CAA")
names(gln) <- opt
glu <- c("GAG", "GAA")
names(glu) <- opt
tyr <- c("TAC", "TAT")
names(tyr) <- opt
his <- c("CAC", "CAT")
names(his) <- opt
asp <- c("GAC", "GAT")
names(asp) <- opt

perc_asn <- cu_human[[asn[["MO"]]]] /
  (cu_human[[asn[["MO"]]]] + cu_human[[asn[["LO"]]]])
perc_lys <- cu_human[[lys[["MO"]]]] /
  (cu_human[[lys[["MO"]]]] + cu_human[[lys[["LO"]]]])
perc_phe <- cu_human[[phe[["MO"]]]] /
  (cu_human[[phe[["MO"]]]] + cu_human[[phe[["LO"]]]])
perc_cys <- cu_human[[cys[["MO"]]]] /
  (cu_human[[cys[["MO"]]]] + cu_human[[cys[["LO"]]]])
perc_gln <- cu_human[[gln[["MO"]]]] /
  (cu_human[[gln[["MO"]]]] + cu_human[[gln[["LO"]]]])
perc_glu <- cu_human[[glu[["MO"]]]] /
  (cu_human[[glu[["MO"]]]] + cu_human[[glu[["LO"]]]])
perc_tyr <- cu_human[[tyr[["MO"]]]] /
  (cu_human[[tyr[["MO"]]]] + cu_human[[tyr[["LO"]]]])
perc_his <- cu_human[[his[["MO"]]]] /
  (cu_human[[his[["MO"]]]] + cu_human[[his[["LO"]]]])
perc_asp <- cu_human[[asp[["MO"]]]] /
  (cu_human[[asp[["MO"]]]] + cu_human[[asp[["LO"]]]])

lower_asn <- binom.wilson(cu_human[[asn[["MO"]]]], (cu_human[[asn[["MO"]]]] + cu_human[[asn[["LO"]]]]))$lower
lower_lys <- binom.wilson(cu_human[[lys[["MO"]]]], (cu_human[[lys[["MO"]]]] + cu_human[[lys[["LO"]]]]))$lower
lower_phe <- binom.wilson(cu_human[[phe[["MO"]]]], (cu_human[[phe[["MO"]]]] + cu_human[[phe[["LO"]]]]))$lower
lower_cys <- binom.wilson(cu_human[[cys[["MO"]]]], (cu_human[[cys[["MO"]]]] + cu_human[[cys[["LO"]]]]))$lower
lower_gln <- binom.wilson(cu_human[[gln[["MO"]]]], (cu_human[[gln[["MO"]]]] + cu_human[[gln[["LO"]]]]))$lower
lower_glu <- binom.wilson(cu_human[[glu[["MO"]]]], (cu_human[[glu[["MO"]]]] + cu_human[[glu[["LO"]]]]))$lower
lower_tyr <- binom.wilson(cu_human[[tyr[["MO"]]]], (cu_human[[tyr[["MO"]]]] + cu_human[[tyr[["LO"]]]]))$lower
lower_his <- binom.wilson(cu_human[[his[["MO"]]]], (cu_human[[his[["MO"]]]] + cu_human[[his[["LO"]]]]))$lower
lower_asp <- binom.wilson(cu_human[[asp[["MO"]]]], (cu_human[[asp[["MO"]]]] + cu_human[[asp[["LO"]]]]))$lower

upper_asn <- binom.wilson(cu_human[[asn[["MO"]]]], (cu_human[[asn[["MO"]]]] + cu_human[[asn[["LO"]]]]))$upper
upper_lys <- binom.wilson(cu_human[[lys[["MO"]]]], (cu_human[[lys[["MO"]]]] + cu_human[[lys[["LO"]]]]))$upper
upper_phe <- binom.wilson(cu_human[[phe[["MO"]]]], (cu_human[[phe[["MO"]]]] + cu_human[[phe[["LO"]]]]))$upper
upper_cys <- binom.wilson(cu_human[[cys[["MO"]]]], (cu_human[[cys[["MO"]]]] + cu_human[[cys[["LO"]]]]))$upper
upper_gln <- binom.wilson(cu_human[[gln[["MO"]]]], (cu_human[[gln[["MO"]]]] + cu_human[[gln[["LO"]]]]))$upper
upper_glu <- binom.wilson(cu_human[[glu[["MO"]]]], (cu_human[[glu[["MO"]]]] + cu_human[[glu[["LO"]]]]))$upper
upper_tyr <- binom.wilson(cu_human[[tyr[["MO"]]]], (cu_human[[tyr[["MO"]]]] + cu_human[[tyr[["LO"]]]]))$upper
upper_his <- binom.wilson(cu_human[[his[["MO"]]]], (cu_human[[his[["MO"]]]] + cu_human[[his[["LO"]]]]))$upper
upper_asp <- binom.wilson(cu_human[[asp[["MO"]]]], (cu_human[[asp[["MO"]]]] + cu_human[[asp[["LO"]]]]))$upper

codon_usage_df <- data.frame(AA = c("N", "K", "F", "C", "Q", "E", "Y", "H", "D")) %>% mutate(
  percentage_of_optimal_codons = case_when(
    AA == "N" ~ perc_asn,
    AA == "K" ~ perc_lys,
    AA == "F" ~ perc_phe,
    AA == "C" ~ perc_cys,
    AA == "Q" ~ perc_gln,
    AA == "E" ~ perc_glu,
    AA == "Y" ~ perc_tyr,
    AA == "H" ~ perc_his,
    AA == "D" ~ perc_asp,
  ),
  confint_upper = case_when(
    AA == "N" ~ upper_asn,
    AA == "K" ~ upper_lys,
    AA == "F" ~ upper_phe,
    AA == "C" ~ upper_cys,
    AA == "Q" ~ upper_gln,
    AA == "E" ~ upper_glu,
    AA == "Y" ~ upper_tyr,
    AA == "H" ~ upper_his,
    AA == "D" ~ upper_asp,
  ),
  confint_lower = case_when(
    AA == "N" ~ lower_asn,
    AA == "K" ~ lower_lys,
    AA == "F" ~ lower_phe,
    AA == "C" ~ lower_cys,
    AA == "Q" ~ lower_gln,
    AA == "E" ~ lower_glu,
    AA == "Y" ~ lower_tyr,
    AA == "H" ~ lower_his,
    AA == "D" ~ lower_asp,
  )
)

# Sort all amino acids by codon usage bias
by_cu <- codon_usage_df %>%
  arrange(desc(percentage_of_optimal_codons)) %>%
  select(AA) %>%
  unlist()
codon_usage_df$AA <- factor(codon_usage_df$AA, levels = by_cu)

codon_usage_df %>%
  ggplot(aes(
    x = factor(AA), y = percentage_of_optimal_codons
  )) +
  ylab("Percentage of optimal codons") +
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
  ) -> mo2lo_cu

ggsave(mo2lo_cu, width = 9, height = 6, filename = snakemake@output[["plot"]])
