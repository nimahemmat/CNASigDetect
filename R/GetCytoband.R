#' Get cytoband table (packaged)
#' @export
GetCytoband <- function(genome = c("hg38","hg19")) {
  genome <- match.arg(genome)
  if (genome == "hg38") return(hg38_cytoband)
  else return(hg19_cytoband)
}
