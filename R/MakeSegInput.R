#' @title MakeSegInput()
#'
#' @description
#' This function renames user-specified columns to a standard format and computes additional variables useful for copy number analysis.
#'
#' @param Segments A data frame of allele-specific segments for all samples.
#' @param SampleID_col Name of the column which contains sample IDs.
#' @param Chr_col Name of the column which contains chromosome names.
#' @param Start_col Name of the column which contains start position of segments.
#' @param End_col Name of the column which contains end position of segments.
#' @param MajorCN_col Name of the column which contains major copy number of segments.
#' @param MinorCN_col Name of the column which contains minor copy number of segments.
#'
#' @importFrom magrittr %>%
#'
#' @return A tibble with standardised columns: SampleID, Chr, Start, End, MinorCN, CopyNumber, len_bp, len_mb.
#' @export
MakeSegInput <- function(Segments,
                         SampleID_col,
                         Chr_col,
                         Start_col,
                         End_col,
                         MajorCN_col,
                         MinorCN_col) {
  # check that df provided
  if (missing(Segments) || !is.data.frame(Segments)) {
    stop("'Segments' must be a data frame")
  }

  # collect user-provided column names
  cols <- c(SampleID_col, Chr_col, Start_col, End_col, MajorCN_col, MinorCN_col)
  names(cols) <- c("SampleID", "Chr", "Start", "End", "MajorCN", "MinorCN")

  # check if all exist
  missing_cols <- cols[!cols %in% colnames(Segments)]
  if (length(missing_cols) > 0) {
    stop("These columns were not found in Segments: ",
         paste(unique(missing_cols), collapse = ", "))
  }

  # standardise to consistent names
  seg <- Segments %>%
    dplyr::rename(
      SampleID = !!rlang::sym(SampleID_col),
      Chr      = !!rlang::sym(Chr_col),
      Start    = !!rlang::sym(Start_col),
      End      = !!rlang::sym(End_col),
      MajorCN  = !!rlang::sym(MajorCN_col),
      MinorCN  = !!rlang::sym(MinorCN_col)
    ) %>%
    dplyr::mutate(
      CopyNumber = MajorCN + MinorCN,
      len_bp     = End - Start,
      len_mb     = (End - Start)/1e6
    )

  return(seg)
}
