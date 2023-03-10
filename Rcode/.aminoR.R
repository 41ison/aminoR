# packages and data
library(tidyverse)
library(ggrepel)
library(kableExtra)
library(seqinr)
aminoacids <- read_tsv("./AAdata.tsv")

# Tabela resumindo as características dos aminoácidos
# summary of aminoacids table
aminoacids %>% 
    select(Name,Abbr,AA, Hydropathy, pKa, pKb, pKx, pl,
            Monoisotopic_Mass, frequency_in_proteins) %>%
    arrange(Hydropathy) %>%
     rename("Aminoacid" = Name,
            "pKa (COOH)" = pKa,
            "pKb (NH2)" = pKb,
            "pKx (side chain)" = pKx,
            "Monoisotopic Mass" = Monoisotopic_Mass,
            "Frequency in proteins (%)" = frequency_in_proteins) %>%
    kbl(caption = "Summary of amino acids properties",
      align = "c") %>% 
  kable_classic(full_width = F, html_font = "Cambria")

# Functions to calculate the pKa of a peptide, given the one letter code for each amino acid
pKa <- function(x) {
    if (x == "A" || x =="G") {
        return(2.34)
    } else if (x == "V") {
        return(2.32)
    } else if (x == "L" || x == "I") {
        return(2.36)
    } else if (x == "P") {
        return(1.99)
    } else if (x == "F") {
        return(1.83)
    } else if (x == "Y") {
        return(2.20)
    } else if (x == "W") {
        return(2.83)
    } else if (x == "C") {
        return(1.96)
    } else if (x == "M") {
        return(2.28)
    } else if (x == "S") {
        return(2.21)
    } else if (x == "T") {
        return(2.09)
    } else if (x == "N") {
        return(2.02)
    } else if (x == "Q") {
        return(2.17)
    } else if (x == "D") {
        return(1.88)
    } else if (x == "E") {
        return(2.19)
    } else if (x == "K") {
        return(2.18)
    } else if (x == "R") {
        return(2.17)
    } else if (x == "H") {
        return(1.82)
    }
}

pKb <- function(x) {
    if (x == "A") {
        return(9.69)
    } else if (x == "G" || x == "L" ||
               x == "I" || x == "D") {
        return(9.60)
    } else if (x == "C") {
        return(10.28)
    } else if (x == "M") {
        return(9.21)
    } else if (x == "S") {
        return(9.15)
    } else if (x == "T") {
        return(9.10)
    } else if (x == "N") {
        return(8.80)
    } else if (x == "Q" || x == "F") {
        return(9.13)
    } else if (x == "E") {
        return(9.67)
    } else if (x == "K") {
        return(8.95)
    } else if (x == "R") {
        return(9.04)
    } else if (x == "H") {
        return(9.17)
    } else if (x == "V") {
        return(9.62)
    } else if (x == "P") {
        return(10.60)
    } else if (x == "Y") {
        return(9.11)
    } else if (x == "W") {
        return(9.39)
    }
}

pKx <- function(x) {
    if (x == "A" || x == "G" || x == "V" ||
        x == "L" || x == "I" || x == "P" ||
        x == "F" || x == "W" || x == "N" ||
        x == "Q" || x == "S" || x == "T" ||
        x == "M") {
        return(NA)
    } else if (x == "C") {
        return(8.18)
    } else if (x == "Y") {
        return(10.07)
    } else if (x == "D") {
        return(3.65)
    } else if (x == "E") {
        return(4.25)
    } else if (x == "K") {
        return(10.53)
    } else if (x == "R") {
        return(12.48)
    } else if (x == "H") {
        return(6.00)
    }
}

# Function to calculate the monoisotopic mass of a peptide, given the one letter code amino acid
monoisotopic_mass <- function(x) {
    if (x == "A") {
        return(71.03711)
    } else if (x == "C") {
        return(103.00919)
    } else if (x == "D") {
        return(115.02694)
    } else if (x == "E") {
        return(129.04259)
    } else if (x == "F") {
        return(147.06841)
    } else if (x == "G") {
        return(57.02146)
    } else if (x == "H") {
        return(137.05891)
    } else if (x == "I") {
        return(113.08406)
    } else if (x == "K") {
        return(128.09496)
    } else if (x == "L") {
        return(113.08406)
    } else if (x == "M") {
        return(131.04049)
    } else if (x == "N") {
        return(114.04293)
    } else if (x == "P") {
        return(97.05276)
    } else if (x == "Q") {
        return(128.05858)
    } else if (x == "R") {
        return(156.10111)
    } else if (x == "S") {
        return(87.03203)
    } else if (x == "T") {
        return(101.04768)
    } else if (x == "V") {
        return(99.06841)
    } else if (x == "W") {
        return(186.07931)
    } else if (x == "Y") {
        return(163.06333)
    }
}

# function to calculate the isoelectric point for the amino acid, given the one letter code
pI <- function(x) {
    if (x == "A") {
        return(6.00)
    } else if (x == "G") {
        return(5.97)
    } else if (x == "V") {
        return(5.96)
    } else if (x == "L") {
        return(5.98)
    } else if (x == "I") {
        return(6.02)
    } else if (x == "P") {
        return(6.30)
    } else if (x == "F") {
        return(5.48)
    } else if (x == "Y") {
        return(5.66)
    } else if (x == "W") {
        return(5.89)
    } else if (x == "C") {
        return(5.07)
    } else if (x == "M") {
        return(5.74)
    } else if (x == "T") {
        return(5.60)
    } else if (x == "N") {
        return(5.41)
    } else if (x == "D") {
        return(2.77)
    } else if (x == "E") {
        return(3.22)
    } else if (x == "Q") {
        return(5.65)
    } else if (x == "S") {
        return(5.68)
    } else if (x == "R") {
        return(10.76)
    } else if (x == "K") {
        return(9.74)
    } else if (x == "H") {
        return(7.59)
    }
}

# calculate the monoisotopic mass of a peptide
monoisotopic_pepetide <- function(x) {
    sum <- 0
    for (i in 1:length(x)) {
        sum <- sum + monoisotopic_mass(x[i])
    }
    return(sum)
}

# how to use the monoisotopic_pepetide function
# Natriuretic peptides B (P16860) sequence for testing
AA <- c("M","D","P","Q","T","A","P","S","R","A","L","L",
"L","L","L","F","L","H","L","A","F","L","G","G","R","S","H",
"P","L","G","S","P","G","S","A","S","D","L","E","T","S","G","L","Q","E","Q",
"R","N","H","L","Q","G","K","L","S","E","L","Q","V","E","Q","T","S","L","E",
"P","L","Q","E","S","P","R","P","T","G","V","W","K","S","R","E","V","A","T",
"E","G","I","R","G","H","R","K","M","V","L","Y","T","L","R","A","P","R","S",
"P","K","M","V","Q","G","S","G","C","F","G","R","K","M","D","R","I","S",
"S","S","S","G","L","G","C","K","V","L","R","R","H")
monoisotopic_pepetide(AA)

# plot the pKa values for a given sequence of amino acids
pKa_plot <- function(x) {
    ggplot(aminoacids, aes(x = Abbr)) +
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
        labs(x = "Aminoacids", y = "pKa") +
    theme_classic(base_size = 20)
}

# using the function to plot the pKa
# it is necessary to have a data.frame containing the variables pKa, pKb, and pKx
pKa_plot(aminoacids)

# function to plot the monoisotopic mass of the aminoacids coloring by hydropathy
monoiso_plot <- function(x) {
    ggplot(aminoacids, aes(x = Abbr, y = Monoisotopic_Mass,
                    size = Monoisotopic_Mass, color = Hydropathy)) +
            scale_color_manual(values = c("#0073C2FF", "#CD534CFF", "#8F7700FF")) +
    geom_point(alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size=6)), size = FALSE) +
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
}

# using the function for plot monoisotopic masses and hydropathy
monoiso_plot(aminoacids)
