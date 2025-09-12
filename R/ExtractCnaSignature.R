#' @title ExtractCnaSignature()
#'
#' @description
#' Function to extract CNA signatures based on df created by MakeSegInput() function.
#'
#' @param Input.Segments Segment data frame created by MakeSegInput()
#' @param Focal_th Threshold of length for focal alterations (default = 3 Mb)
#' @param Broad_th Threshold of length for broad alterations (default = 10 Mb)
#' @param Art.len Threshold for artifact segments (default = 0.05 Mb)
#' @param k_range K range for NMF
#' @param n_run number of runs for NMF
#'
#' @importFrom magrittr %>%
#'
#' @return A list of Signatures with their annotation and contribution to each sample
#' @export
ExtractCnaSignature <- function(Input.Segments,
                                Focal_th = 3,
                                Broad_th = 10,
                                Art.len = 0.05,
                                k_range = 2:10,
                                n_run = 100) {

  # Tokenisation of segments' features
  segment_features <- Input.Segments %>%
    dplyr::filter(len_mb > Art.len) %>% # Exclude artifact segments
    dplyr::mutate(
      # Categorisation of lenghts
      len_cat = dplyr::case_when(
        len_mb <= Focal_th ~ "Focal",
        len_mb > Focal_th & len_mb < Broad_th ~ "Medium",
        TRUE ~ "Broad"
      ),
      # Categorisation of Copy numbers
      cn_cat = dplyr::case_when(
        CopyNumber == 0 ~ "DEL",
        CopyNumber == 1 ~ "LOH",
        CopyNumber == 2 & MinorCN == 1 ~ "NEU",
        CopyNumber == 2 & MinorCN == 0 ~ "LOH",
        CopyNumber %in% c(3, 4, 5, 6) & MinorCN != 0 ~ "Gain",
        CopyNumber %in% c(3, 4, 5, 6) & MinorCN == 0 ~ "Gain-LOH",
        CopyNumber %in% c(7, 8, 9) & MinorCN != 0 ~ "ModAMP",
        CopyNumber %in% c(7, 8, 9) & MinorCN == 0 ~ "ModAMP-LOH",
        CopyNumber > 9 & MinorCN != 0 ~ "HiAMP",
        TRUE ~ "HiAMP-LOH"
      ),
      # Make token per segment
      token = paste0(len_cat, "_", cn_cat),
      token = dplyr::if_else(grepl("NEU", token), "NEU", token)
    )

  # Prepare matrix for NMF
  mat_features <- segment_features %>%
    dplyr::group_by(SampleID, token) %>%
    dplyr::summarise(len_mb = sum(len_mb), .groups = "drop") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(frac = len_mb / sum(len_mb)) %>%
    dplyr::ungroup() %>%
    dplyr::select(.,-len_mb) %>%
    tidyr::pivot_wider(names_from = token, values_from = frac, values_fill = 0) %>%
    tibble::column_to_rownames("SampleID") %>%
    t() %>% as.matrix()

  # run NMF algorithm
  set.seed(1)
  # Test different Ks to find best one
  nmf_test <- NMF::nmf(mat_features, k_range, nrun = n_run, method = "brunet", .opt = "v")

  # Chose best k based on residual sum of squares (reconstruction error)
  print("Estimating the best rank ...")
  rss <- nmf_test[["measures"]][["rss"]]
  diffs_k <- data.frame(i = integer(), d1 = numeric(), d2 = numeric())
  for(i in 1:(length(rss) - 2)) {
    d1 <- rss[i] - rss[i + 1]
    d2 <- rss[i + 1] - rss[i + 2]
    diffs_k <- rbind(diffs_k, data.frame(i = i, d1 = d1, d2 = d2))
  }
  diffs_k <- diffs_k %>%
    dplyr::mutate(Diff = d1 - d2) %>%
    dplyr::filter(Diff == stats::median(Diff))

  best_k <- diffs_k$i + 1

  # Run NMF based on best k
  nmf_fit <- NMF::nmf(mat_features, rank = best_k, method = "brunet", nrun = n_run)

  signatures_tokens <- NMF::basis(nmf_fit)
  samples_signatures <- NMF::coef(nmf_fit)

  # Normalise signatures features and signatures per sample
  s <- colSums(signatures_tokens)
  signatures_tokens_normalised <- sweep(signatures_tokens, 2, s, "/") %>% as.data.frame()
  samples_signatures_normalised <- sweep(samples_signatures, 1, s, "*") %>% t() %>% as.data.frame()

  # Set signature names
  colnames(signatures_tokens_normalised) <- paste("Signature", 1:best_k, sep = "")
  colnames(samples_signatures_normalised) <- paste("Signature", 1:best_k, sep = "")

  # Annotation of signatures
  annotate_signatures <- function(sig_df) {
    # Catch pattern of tokens
    addToken <- function(pat) rowSums(sig_df[ , grepl(pat, colnames(sig_df)), drop = FALSE])
    addToken2 <- function(pat1, pat2) {
      sel <- grepl(pat1, colnames(sig_df)) & grepl(pat2, colnames(sig_df))
      rowSums(sig_df[ , sel, drop = FALSE])}

    # Find patterns using defined functions
    BroadHiAMP <- addToken("Broad_HiAMP"); MediumHiAMP <- addToken("Medium_HiAMP"); FocalHiAMP <- addToken("Focal_HiAMP")
    BroadModAMP <- addToken("Broad_ModAMP"); MediumModAMP <- addToken("Medium_ModAMP"); FocalModAMP <- addToken("Focal_ModAMP")
    BroadGain <- addToken("Broad_Gain"); MediumGain <- addToken("Medium_Gain"); FocalGain <- addToken("Focal_Gain")
    BroadDEL <- addToken("Broad_DEL"); MediumDEL <- addToken("Medium_DEL"); FocalDEL <- addToken("Focal_DEL")
    BroadLOH <- addToken2("Broad", "LOH"); MediumLOH <- addToken2("Medium", "LOH"); FocalLOH <- addToken2("Focal", "LOH")
    NEU <- addToken("NEU")

    # Make annotation data frame
    anno_df <- cbind(`Broad Hi-AMP` = BroadHiAMP, `Medium Hi-AMP` = MediumHiAMP, `Focal Hi-AMP` = FocalHiAMP,
                     `Broad Mod-AMP` = BroadModAMP, `Medium Mod-AMP` = MediumModAMP, `Focal Mod-AMP` = FocalModAMP,
                     `Broad Gain` = BroadGain, `Medium Gain` = MediumGain, `Focal Gain` = FocalGain,
                     `Broad DEL` = BroadDEL, `Medium DEL` = MediumDEL, `Focal DEL` = FocalDEL,
                     `Broad LOH` = BroadLOH, `Medium LOH` = MediumLOH, `Focal LOH` = FocalLOH,
                     NEU = NEU)
    out <- anno_df %>%
      as.data.frame() %>%
      dplyr::mutate(Signature = rownames(sig_df)) %>%
      dplyr::relocate(Signature, .before = 1)
  }

  signatures_annotation <- signatures_tokens_normalised %>%
    t() %>% as.data.frame() %>%
    annotate_signatures()

  # Gather results as a list
  SignaturesResults <- list(
    Sample_Weights = samples_signatures_normalised,
    Signature_Annotation = signatures_annotation
  )

  return(SignaturesResults)

}
