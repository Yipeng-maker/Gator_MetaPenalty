#' The meta-analysis by Bohren et al. (2017) to investigate the effects of spontaneous vaginal
#' births, and continuous, one-to-one intrapartum support compared with usual care.
#'
#' A dataset containing 21 studies, with the number of patients in treated and control groups.
#' The number of events for each group is also reported for each study. The log odds ratio and
#' its corresponding variance are reported or caculated for each study.
#'
#' @format A data frame with 21 rows and 6 variables:
#' \describe{
#'   \item{n1}{the number of women received continuous support during childbirth}
#'   \item{n2}{the number of women received usual care during childbirth}
#'   \item{r1}{the number of spontaneous vaginal birth in the treated group}
#'   \item{r2}{the number of spontaneous vaginal birth in the control group}
#'   \item{y}{the log odds ratio}
#'   \item{s2}{the corresponding within-study variance}
#' }
#' @source \url{https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD003766.pub6/epdf/full}
"Bohren"









