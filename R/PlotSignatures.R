#' @title PlotSignatures()
#'
#' @description
#' Function to generate overview plot of identified signatures
#'
#' @param Sig.Results List of results created by ExtractCnaSignature()
#'
#' @importFrom magrittr %>%
#'
#' @return List of bar plots for each identified signatures
#' @export
PlotSignatures <- function(Sig.Results) {
  # Extract features per signature
  features <- Sig.Results$Signature_Annotation %>%
    dplyr::select(.,-Signature) %>%
    t() %>% as.data.frame()

  features$Features <- rownames(features)
  features <- dplyr::relocate(features, Features)
  features <- dplyr::mutate(features, dplyr::across(-Features, ~ as.numeric(.)))

  # Plot each signature and store in a list of graphs
  plots_list <- list()

  for(sig in colnames(features %>% dplyr::select(.,-Features))) {
    s <- sig
    sig_plot <- ggplot2::ggplot(features,
                                ggplot2::aes(x = Features,
                                             y = .data[[s]],
                                             fill = Features)) +
      ggplot2::geom_col() +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  breaks = c(0, 0.25, 0.5, 0.75, 1),
                                  expand = c(0, 0)) +
      ggplot2::labs(title = s, y = "Weight") +
      ggplot2::theme_light(base_size = 12) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white"),
        panel.border = ggplot2::element_rect(color = "black"),
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(color = "black", angle = 90, size = 8, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(color = "black", size = 10),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_line(color = "black"),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(color = "black", size = 11),
        plot.title = ggplot2::element_text(color = "black", size = 13, hjust = 0.5),
        legend.position = "none"
      )

    plots_list[[s]] <- sig_plot
  }

  plots_list
}
