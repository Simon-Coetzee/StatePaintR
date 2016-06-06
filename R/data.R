#' Default Segmentation Definitions
#'
#' A table describing the default state definitions used with the default
#' translation layer
#'
#' load with \code{data(Cedars.BFG.states)}
#'
#' @return \code{matrix} of states
#' @format A matrix with 20 states and 8 marks:
#' \describe{
#'   \item{TRS}{Transcribed Regions}
#'   \item{HET}{Heterochromatin Regions}
#'   \item{SCR}{Silenced Regions}
#'   \item{PAR}{Promotor Active Regions}
#'   \item{PARC}{Promotor Active Regions, Core}
#'   \item{PPR}{Promotor Poised Regions}
#'   \item{PPRC}{Promotor Poised Regions, Core}
#'   \item{EAR}{Enhancer Active Regions}
#'   \item{EARC}{Enhancer Active Regions, Core}
#'   \item{EPR}{Enhancer Poised Regions}
#'   \item{EPRC}{Enhancer Poised Regions, Core}
#'   \item{ER}{Enhancer Regions}
#'   \item{ERC}{Enhancer Regions, Core}
#'   \item{PR}{Promoter Regions}
#'   \item{PRC}{Promoter Regions, Core}
#'   \item{CTCF}{CTCF Regions}
#'   \item{CTCFC}{CTCF Regions, Core}
#'   \item{RPS}{Regulatory Putative Region (DNase1 Hypersensitivity, ATAC-Seq, or TF binding)}
#'   \item{AR}{Active Region}
#'   \item{ARC}{Active Region, Core}
#'   }
#' @details Active vs. Poised Regions are called based upon the presence or absence of
#'   H3K27ac, if that mark is present in the data set.
#'   If H3K27ac is not present in the data set, Active vs. Poised is not called.
#'
#' @source \url{"http://github.com/Simon-Coetzee/StateHub/"}
#' @examples
#' data(Cedars.BFG.states)
#' Cedars.BFG.states

"Cedars.BFG.states"
