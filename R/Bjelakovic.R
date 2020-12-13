#' The meta-analysis conducted by Bjelakovic et al. (2014) to assess beneficial and harmful
#' effects of vitamin D supplementation for the prevention of mortality in healthy adults and 
#' adults in a stable phase of disease
#'
#' A dataset containing 53 studies, with the number of patients in treated and control groups.
#' The number of events for each group is also reported for each study. The log odds ratio and
#' its corresponding variance are reported or caculated for each study.
#'
#' @format A data frame with 53 rows and 6 variables:
#' \describe{
#'   \item{n1}{the number of people received Vitamin D}
#'   \item{n2}{the number of people received placebo}
#'   \item{r1}{the number of all-cause mortality in the treated group}
#'   \item{r2}{the number of all-cause mortality in the control group}
#'   \item{y}{the log odds ratio}
#'   \item{s2}{the corresponding within-study variance}
#' }
#' @source \url{https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD007470.pub3/epdf/full}
"Bjelakovic"




