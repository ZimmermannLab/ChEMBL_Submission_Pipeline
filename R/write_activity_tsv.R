write_activity_tsv <- function(output_dir, input_csv, compound_tsv, assay_tsv,
                               ridx, category = c("bacteria", "enzyme", "community"), type = "Biotransformation"){
  # INPUT CSV SHOULD HAVE ATLEAST THE FOLLOWING COLUMNS: compound names, bacteria/strain (when applied), binary_activity, activity comment;
  # For additional ChEMBL defined columns, please use exactly the same name of the column as mentioned by ChEMBL!



}


input_csv <- "/g/zimmermann/Members/zulfiqar/ZM000_PistoiaPrj/ChEMBL_Submission_Pipeline/ChEMBL_Submission_Pipeline/inputs/activity_input.csv"
input_activity <- read.csv(input_csv)
compound_tsv <- "/g/zimmermann/Members/zulfiqar/ZM000_PistoiaPrj/ChEMBL_Submission_Pipeline/ChEMBL_Submission_Pipeline/outputs/COMPOUND_RECORD.tsv"
input_compound <- read.csv(compound_tsv, sep = "\t")
assay_tsv <- "/g/zimmermann/Members/zulfiqar/ZM000_PistoiaPrj/ChEMBL_Submission_Pipeline/ChEMBL_Submission_Pipeline/outputs/ASSAY.tsv"
input_assay <- read.csv(assay_tsv)

ACTIVITY <- data.frame(matrix(ncol = 20, nrow = nrow(input_activity)))
colnames(ACTIVITY) <- c("CIDX", "CRIDX", "SRC_ID_CIDX", "AIDX", "SRC_ID_AIDX", "RIDX", "TEXT_VALUE",
                     "RELATION", "VALUE", "UPPER_VALUE", "UNITS", "SD_MINUS",
                     "SD_PLUS", "ACTIVITY_COMMENT", "CRIDX_CHEMBLID", "CRIDX_DOCID",
                     "ACT_ID", "TEOID", "TYPE", "ACTION_TYPE")
ACTIVITY$CRIDX <- ridx
ACTIVITY$TYPE <- "Biotransformation"
ACTIVITY$ACTIVITY_COMMENT <- input_activity$Comment
ACTIVITY[which(input_activity$Biotransformation == "Compound Metabolized"), "ACTION_TYPE"] <- "Compound Metabolized"
ACTIVITY$CIDX <- input_compound$CIDX[match(input_activity$DrugName, input_compound$COMPOUND_NAME)]


# Pre-allocate with NA
input_activity$matched_AIDX <- NA

# Correctly splitting and extracting the first two parts of the Sample column
input_activity$first_split <- sapply(input_activity$Species, function(x) {
  if (grepl(" ", x)) {
    # Split by space and extract the first two parts
    parts <- strsplit(x, " ")[[1]]  # Split the string
    if (length(parts) >= 2) {
      return(paste(parts[1:2], collapse = " "))  # Return the first two parts as a single string
    } else {
      return(parts[1])  # Return the only part if there's less than two
    }
  } else {
    return(x[1])  # Return the original string if no spaces are found
  }
})



input_activity$first_split



# Use sapply to vectorize the matching process
input_activity$matched_AIDX <- sapply(input_activity$first_split, function(split_string) {

  # Find matching AIDX in input_assay for the current split_string
  match_idx <- grepl(split_string, input_assay$AIDX)

  # If a match is found, return the first matching AIDX, otherwise return NA
  if (any(match_idx)) {
    return(input_assay$AIDX[match_idx][1])  # Return the first match
  } else {
    return(NA)
  }
})

# Print the updated input_activity
print(input_activity)






# Correctly splitting and extracting the first two parts of the Sample column
input_activity$first_split <- sapply(input_activity$Species, function(x) {
  print(x)  # Print each value of 'Sample' for debugging

  if (grepl(" ", x)) {
    # Split by space and extract the first two parts
    parts <- strsplit(x, " ")[[1]]  # Split the string

    if (length(parts) >= 2) {
      # Return the first two parts joined by a space
      result <- paste(parts[1:2], collapse = " ")
      return(result)
    } else {
      # If there's only one part, return that
      return(parts[1])
    }
  } else {
    # No spaces found, return the original string
    print(x)  # Print the original string
    return(x)
  }
})

# Print the resulting first_split column to verify
print(input_activity$first_split)



