# packages and data
library(tidyverse)
library(ggrepel)
library(kableExtra)
aminoacids <- read_tsv("./aminoacidsTab.tsv")

# Tabela resumindo as características dos aminoácidos
# summary of aminoacids table
aminoacids %>% 
    select(Abbr, Name, Hydropathy, pKa, pKb, pKx, Monoisotopic_Mass) %>%
    arrange(Hydropathy) %>%
    kbl(caption = "Summary of amino acids properties)",
      align = "c") %>% 
  kable_classic(full_width = F, html_font = "Cambria")

# pKa dos diferentes aminoácidos para cada região (COOH, NH2 e cadeia lateral)
# pKa of different amino acids for each group (COOH, NH2 and R-group)
pKa <- ggplot(aminoacids, aes(x = Abbr)) +
     geom_line(aes(y = pKa, group = 1), color = "#0073C2FF") +
        annotate("text", x = 16, y = 2.7, label = "-COOH",
            color = "#0073C2FF", size = 6) +
     geom_line(aes(y = pKb, group = 2), color = "#CD534CFF") +
        annotate("text", x = 15, y = 11, label = "-NH2",
            color = "#CD534CFF", size = 6) +
     geom_col(aes(y = pKx, group = 3), fill = "#4A6990FF", alpha = 0.4) +
        annotate("text", x = 9, y = 6.8, label = "R-group",
            color = "#4A6990FF", size = 6) +
     scale_y_continuous(limits = c(0, 12.5), breaks = seq(0, 12.5, 2)) +
        labs(x = "Aminoacids", y = "pKa",
             title = "Aminoacids properties") +
    theme_classic(base_size = 20)

ggsave("pKa.png", pKa, width = 10, height = 6, dpi = 300,
    path = "/Volumes/Expansion/proteomics_basics/")

# distribuição das massas monoisotópicas dos aminoácidos
monoiso <- ggplot(aminoacids, aes(x = Abbr, y = Monoisotopic_Mass,
                    size = Monoisotopic_Mass, color = Hydropathy)) +
                scale_color_manual(values = c("#0073C2FF", "#CD534CFF", "#8F7700FF")) +
    geom_point(alpha = 0.5) +
    guides(colour = "legend", size = FALSE) +
    geom_segment(aes(x = Abbr, xend = Abbr,
                     y = 50, yend = Monoisotopic_Mass),
                     linewidth = 0.2, linetype = "dashed",
                     show.legend = FALSE) +
    geom_text_repel(aes(label = Monoisotopic_Mass),
        color = "#868686FF", size = 5, nudge_x = 2,
            nudge_y = 2, show.legend = FALSE) +
    scale_size(range = c(6, 20)) +
        scale_y_continuous(limits = c(50, 190),
        breaks = seq(50, 190, 10)) +
    labs(x = "Aminoacid", y = "Monoisotopic mass",
         title = "Monoisotopic mass and charge of aminoacids") +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90,
                        vjust = 0.5, hjust = 1), legend.position = "bottom") +
                        facet_wrap(~ Hydropathy, ncol = 3, scales = "free_x")

ggsave("Monoisotopic_Mass.png", monoiso, width = 9, height = 8, dpi = 300,
    path = "/Volumes/Expansion/proteomics_basics/")
