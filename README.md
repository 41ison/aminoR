# aminoR
This repository aims to build a complete table of information about amino acids and share some useful codes for summary and graphical presentations. There are a couple of small functions for quick information about each amino acid.

pKa() returns the pKa for COOH group of a given amino acid qeued as one letter code.
pKb() returns the pKa for NH2 group of a given amino acid qeued as one letter code.
pKx() returns the value for R-group of a given amino acid qeued as one letter code.
pI() returns the value of isoelectric point for a given amino acid qeued as one letter code.
monoisotopic_mass() returns the monoisotopic mass for a given amino acid qeued as one letter code.
monoisotopic_pepetide() returns the monoisotopic mass for a given peptide qeued a vector of strings containing individual amino acids.
pKa_plot() produces a plot containing the pKa distribution for each amino acid in the table.
monoiso_plot() produces a plot with monoisotopic mass of the aminoacids coloring by hydropathy.
