
suppressWarnings({
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(writexl)
    library(stringr)
    library(data.table)
  })
})

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

################ FUNCTION TO TAKE OFF THE REDUNDANCY ######################
# function_redundancy <- function(df) {
#   for (i in 2:nrow(df)) {
#     if (df[i, "Contig"] == df[i-1, "Contig"] && 
#         df[i, "coordinates"] == df[i-1, "coordinates"] && 
#         df[i, "strand"] == df[i-1, "strand"]) {
#       # Take the coordinates of start & end the prediction for i (start_2, end_2) and i-1
#       start_1 <- df[i-1,"start"] # ex: 1000
#       end_1 <- df[i-1, "end"] # ex: 2000 
#       start_2 <- df[i,"start"] # ex: 500
#       end_2 <- df[i,"end"]
#       # Just an expression to check if the prediction overlap between them
#       if (start_2 < end_1-200 && end_2 > start_1+200) {
#         # We have to check the e-value before because we have ordered by start of the prediction.
#         e_value_1 <- df[i-1,"e-value"]
#         e_value_2 <- df[i,"e-value"]
#         # Put true if we want to filter out the results
#         if (e_value_1 < e_value_2) {
#           df$compare[i] <- TRUE
#           # In the case both e-value are 0 we keep them
#         } else if (e_value_1 == 0 && e_value_2 == 0) {
#           df$compare[i-1] <- FALSE
#           df$compare[i] <- FALSE
#         } else {
#           df$compare[i-1] <- TRUE
#           df$compare[i] <- FALSE
#         }  
#       } else {
#         df$compare[i] <- FALSE
#       }  
#     } else {
#       df$compare[i] <- FALSE
#     }
#   }
#   return(df)
# }


function_redundancy_dt <- function(dt) {
  colnames(dt) <- gsub("-", "_", colnames(dt))
  setDT(dt)  # Ensure it's a data.table
  
  # Explicitly order by desired columns
  setorder(dt, Contig, coordinates, strand, start)
  
  # Adding the shift/lag columns
  dt[, `:=` (
    start_lag = shift(start, type = "lag"),
    end_lag = shift(end, type = "lag"),
    e_value_lag = shift(e_value, type = "lag")
  ), by = .(Contig, coordinates, strand)]
  
  # Perform the redundancy check
  dt[, compare := {
    compare <- rep(FALSE, .N)  # Initialize with FALSE
    
    if (.N > 1) {
      for (i in 2:.N) {
        if (!is.na(start[i]) && !is.na(end[i]) && !is.na(start_lag[i]) && !is.na(end_lag[i])) {
          if (start[i] < end_lag[i] - 200 && end[i] > start_lag[i] + 200) {
            if (e_value[i] > e_value_lag[i]) {
              compare[i] <- TRUE
            } else if (e_value[i] == 0 && e_value_lag[i] == 0) {
              compare[i] <- FALSE
            } else {
              compare[i-1] <- TRUE
            }
          }
        }
      }
    }
    compare
  }]
  
  # Remove the lag columns before returning
  dt[, c("start_lag", "end_lag", "e_value_lag") := NULL]
  
  return(as.data.frame(dt))
}


#########################################################################################


protein_from_fasta <- read.delim2("filtered_protein.tsv", header = FALSE, sep = "\t")

# Take only the header (grepl is useful to use regex expression into filter command)
predicted_protein <- protein_from_fasta %>% 
  filter(grepl('^>GCA', V1))

colnames(predicted_protein) <- c('Assembly', 'Contig', 'coordinates', 'strand', 'e-value',
                                 'start','end')

# predicted_protein <- predicted_protein[,-c(1)]

# substitute > with nothing
predicted_protein$Assembly <- gsub( '>', '', predicted_protein$Assembly)

# Make a new dataframe to work with it and don't lose the original one
predicted_protein_simplificado <- predicted_protein
predicted_protein_simplificado <- predicted_protein_simplificado %>% 
  mutate_at(c('start', 'end', 'e-value'), as.numeric)

# Sort the data. In this case is important sort for e-value after take into account other 
# variables
# Take unique entry for start and stop of the prediction
# ORder the data of the new dataframe. In this case we can sort just for the start to run the
# function (the function take into account i and i-1 meaning that check one row with the
# previous one)

predicted_protein_simplificado <- predicted_protein_simplificado %>%
  arrange(Assembly, strand, coordinates, `e-value`) %>%
  distinct(coordinates, Assembly, Contig, strand, start, .keep_all = TRUE) %>%
  distinct(coordinates, Assembly, Contig, strand, end, .keep_all = TRUE) %>%
  arrange(Assembly, strand, coordinates, start)

# testing
# predicted_protein_simplificado_2 <- predicted_protein_simplificado
# for (i in 1:10000){
#   start = sample(10000:20000, 1)
#   end = start + sample(2000:5000, 1)
#   strand = sample(c("+", "-"), 1)
#   eval = runif(1, min = 0, max = 1e-2)
#   df <- data.frame(Assembly = paste0("GCA_", sample(1000:1005, 1)), 
#                    Contig = paste0("FR", sample(1000:1005, 1)), 
#                    coordinates = paste0(sample(10000:10003, 1),"-", sample(20000:20003, 1)), strand = strand,
#                    "e-value" = eval, start = start, end = end)
#   colnames(df) <- gsub("\\.", "-", colnames(df))
#   predicted_protein_simplificado_2 <- rbind(predicted_protein_simplificado_2, df)
# }

# Function run
predicted_protein_final <- predicted_protein_simplificado
predicted_protein_final_1 <- function_redundancy_dt(predicted_protein_final) %>%
  filter(compare != TRUE) %>%
  arrange(Assembly, strand, coordinates, start)

# Loop until no more redundancy
while (nrow(predicted_protein_final) != nrow(predicted_protein_final_1)) {
  predicted_protein_final <- predicted_protein_final_1
  predicted_protein_final_1 <- function_redundancy_dt(predicted_protein_final) %>%
    filter(compare != TRUE) %>%
    arrange(Assembly, strand, coordinates, start)
}

predicted_protein_final_1 <- predicted_protein_final_1 %>%
  dplyr::select(-compare)
write_xlsx(predicted_protein_final_1, "final_protein.xlsx", col_names = TRUE)


# Change all the columns as.character to join them again (for seqkit grep)
predicted_protein_final_1 <- predicted_protein_final_1 %>% 
  mutate_all(as.character)

# Save the results as text file
text_file <- apply(predicted_protein_final_1[, c(1:7)], 1, FUN = function(i) paste(i, collapse = "_"))
write.table(text_file, 'filtered_protein.txt', row.names = FALSE, col.names = FALSE, quote = F)
