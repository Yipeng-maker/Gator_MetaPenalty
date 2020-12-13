#' The meta-analysis by Carless et al. (2011) to assess the efficacy of platelet-rich-plasmpheresis
#' in reducing peri-operative allogeneic red blood cell transfusion in cardiac surgery.
#'
#' A dataset containing 20 studies, with the number of patients in treated and control groups.
#' The number of events for each group is also reported for each study. The log odds ratio and
#' its corresponding variance are reported or caculated for each study.
#'
#' @format A data frame with 20 rows and 6 variables:
#' \describe{
#'   \item{n1}{the number of patients received PRP}
#'   \item{n2}{the number of patients received placebo}
#'   \item{r1}{the number of patients exposed to allogeneic RBC transfusion in the treated group}
#'   \item{r2}{the number of patients exposed to allogeneic RBC transfusion in the control group}
#'   \item{y}{the log odds ratio}
#'   \item{s2}{the corresponding within-study variance}
#' }
#' @source \url{https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD004172.pub2/epdf/full}
"Carless"
