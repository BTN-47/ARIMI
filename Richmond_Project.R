#CRITICAL, YOU MUST CHANGE THE NAME OF THE INPUT FILE,
#     INPUT FILE NAME IS CURRENTLY "RASPBERRY.TXT" on LINE #27

#Challenge E, 1 protein sequence -> hydropathy values & plot
#     Hydropathy values exported as CSV, "Hydropathy_values.csv"
#     Hydropathy plot exported as PNG, "Hydropathy_plot.png"
#         1. Exported Hydropathy plot gets ugly when AA seq >1000, ggsave did
#         not recognize working directory as a path to save to without error.
#         2. Code is sufficient to handle one FASTA formatted sequence. Any 
#         more and it will forcefully and erroneously splice them together.

#Line purpose is annotated above the lines.
#qiagenbioinformatics.com Kyte-Doolittle values.

#Presupposes the packages are already installed, hopefully not a grave sin.
library(tidyverse)
library(ggplot2)

#Sets up a vector of Kyte-Doolittle Hydropathy values (qiagenbioinformatics.com)
# https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/current/index.php?manual=BE_Protein_hydrophobicity.html
ResidueHydropathy = c("A" =1.80, "C" =2.50, "D" =-3.50, "E" =-3.50, "F" =2.80,
                      "G" =-0.40, "H" =-3.20, "I" =4.50, "K"=-3.90, "L" =3.80,
                      "M" =1.90, "N" =-3.50, "P" =-1.60, "Q" =-3.50, "R" =-4.50,
                      "S" =-0.80, "T" =-0.70, "V" =4.20, "W" =-0.90, "Y" =-1.30)

#Stores simple txt file, MAKE SURE TO CHANGE TO FILE OF CHOICE.
input <- "RASPBERRY.txt"
#Reads all lines, stored as var.
my_data <- readLines(input)

#Any lines that do not start with ">", i.e. a header, will be stored 1:1.
#This simple solution is why multiple sequences will be erroneously spliced.
my_seq <- grep("^[^>].{1,}",my_data, value = TRUE)
#This will combine all seq lines together.
combinedseq <- paste(my_seq,collapse = "")
#This will split all lines apart residue by residue
split_seq <- strsplit(combinedseq, split = NULL)
#Sets up a data frame, column of residues "seq", individual residue per row/cell
seq_frame <- data.frame(seq = unlist(split_seq, use.names = FALSE))

#Residues have Hydropathy values looked up from Hydropathy vector and applied.
#These Hydropathy values are stored in the $hydropathy column of the data frame.
#Any erroneous entries will be assigned NA as Hydropathy value.
seq_frame$hydropathy <- sapply(seq_frame$seq, function(x) {
  sapply(x, function(y) {
    if (y %in% names(ResidueHydropathy)) {
      ResidueHydropathy[[y]]
    } else {
      NA
    }
  })
})
#Sets up a variable range, from 1 to total number of residues
x_range <- 1:nrow(seq_frame)

#exports Hydropathy values as CSV.
write.csv(seq_frame, file = "Hydropathy_values.csv")
#constructs a linegraph of Hydropathy, labelled with input file name, x-axis
#runs the length of the sequence, y-axis are the corresponding hydropathy values (KD)
#Any non-matches of AA are listed as NA (not-plotted in graph, just blank at the spot)
Hydropathy_plot <- ggplot(seq_frame, aes(x = x_range, y=hydropathy)) + 
  geom_path(mapping = aes( x = x_range, y = hydropathy,)) +
  xlab("Position") + ylab("Kyte-Doolittle Hydropathy Scores") +
  labs(title = input)
#Plot is saved as png, ggsave disputed a save path, so png() saves to a rigid
#dimension, not my favorite, gets ugly past 1000 residues.
#dev.off terminates the wizardry that is png()/print().
png("Hydropathy_plot.png",
    width = 1000, height = 500)
print(Hydropathy_plot)
dev.off()
#Check working directory for exported hydropathy values (.csv) and plot (.png)
