suppressWarnings({
  suppressPackageStartupMessages({
    library(readxl)
    library(writexl)
  })
})


args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])


proteins <- readxl::read_xlsx("final_protein.xlsx")
proteins_species <- read.csv('final_protein_species.csv', header = FALSE)
colnames(proteins_species) <- c('Assembly', 'Contig', 'coordinates', 'strand', 'start', 'end', 'Scientific name')
proteins_species$Assembly <- sub( '>', '', proteins_species$Assembly)

final_protein <- merge(proteins, proteins_species, by = c('Assembly', 'Contig','coordinates', 'strand', 'start','end'))
write_xlsx(final_protein, "final_protein_species.xlsx", col_names = TRUE)
