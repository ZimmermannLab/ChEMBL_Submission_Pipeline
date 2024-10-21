

# function to generate reference.tsv
write_reference_tsv <- function(output_file, RIDX, DOI, TITLE, AUTHORS, ABSTRACT, REF_TYPE, YEAR,
                                JOURNAL_NAME = c(), VOLUME = c(), ISSUE = c(),
                                FIRST_PAGE = c(), LAST_PAGE = c(), PATENT_ID = c(), CONTACT = c()){

  # Validate mandatory fields
  if (any(is.na(c(RIDX, TITLE, YEAR, ABSTRACT, AUTHORS, REF_TYPE, DOI)))) {
    stop("RIDX, TITLE, YEAR, ABSTRACT, AUTHORS, REF_TYPE and  DOI are mandatory fields.")
  }

  if (REF_TYPE == "Publication" && is.na(FIRST_PAGE)) {
    stop("FIRST_PAGE is mandatory for if the dataset is a Publication.")
  }

  ReferenceTSV_ChEMBL <- cbind(RIDX, DOI, TITLE, AUTHORS, ABSTRACT, REF_TYPE, YEAR,
                               JOURNAL_NAME, VOLUME, ISSUE,
                               FIRST_PAGE, LAST_PAGE, PATENT_ID, CONTACT)

  write.table(ReferenceTSV_ChEMBL, file=paste(output_dir, "/REFERENCE.tsv", sep = ""), sep='\t', row.names=FALSE)
  return(ReferenceTSV_ChEMBL)
}
