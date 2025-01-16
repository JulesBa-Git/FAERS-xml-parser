##### USAGE ####
# example :
# Rscript csv_to_R_data.R ./FAERSParser_result_csv_file.csv desired_RData_filename(opt.) 

if (length(args) < 1) {
    stop("Csv data file path to be converted should be passed as an argument,
         .Rdata file name desired should also be passed as an argument, default is csv filename", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)

# Read the CSV file into a dataframe
df <- read.csv(args[1], sep = ";", header = TRUE)

# Split the 'patientATC' column by ':'
df <- df %>% 
  mutate(patientATC = strsplit(as.character(patientATC), ":"))

# Convert each list element to a numeric vector
df$patientATC <- lapply(df$patientATC, function(x) as.numeric(x))

# Convert the 'patientADR' column to boolean
df$patientADR <- as.logical(df$patientADR)

RData_name <- if (length(args) >= 2){
    filename <- args[2]
}  else {
    filename_with_extension <- sub(".*/", "", args[1])

    filename <- sub("\\..*", "", filename_with_extension)
}

save.image(paste("./", filename, ".RData", sep = ""))
print("Succesfully converted.")