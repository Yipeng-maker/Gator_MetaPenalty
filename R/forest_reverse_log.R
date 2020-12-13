##' Produce a forest plot of a given meta-analysis when effects are binary outcomes
##'
##' The result plot summarises the point estimate and the 95% confidecne interval for the overall effect size of four models in meta-analysis: the common-effect model, the random-effects model, the penalization method by tuning lambda, and the penalization method by tuning tau. The effect sizes of binary outcomes are consided in this function: log odds ratio and log relative risk. This function reverses the results of effects on log scale to orignal scale 
##' @title Plot a forestplot including the penalization methods
##' @param y an observed effect size vector that is measured on log scale
##' @param s2 the corresponding within-study varaince vector
##' @param out the summary outcome of the metagen( ) in the R package meta
##' @param adjust_font_size a scalar used to control the font size of labels and ticks
##' @param llimit the lower limit for clipping confidence intervals to arrows
##' @param ulimit the upper limit for clipping confidence intervals to arrows
##' @param user_define_x_ticks if it is true, then user-specified x-axis tick marks are used. Otherwise, default ticks are used
##' @param x_ticks the user-specified x-axis tick marks, it must be provided when user_define_x_ticks is TRUE. x_ticks is not specified when user_define_x_ticks is FALSE, and the plot will use the default values
##' @param effect_type the type of effect sizes that are not on the log scale, which can be 'OR' (odds ratio)
##' or 'RR' (relative risk)
##' @author Yipeng
##' @import meta forestplot
##' @export
forest_reverse_log <- function(y, s2, out, adjust_font_size = 0, llimit = 0.25, ulimit = 16, user_define_x_ticks = FALSE, 
    x_ticks, effect_type = "OR") {
    if (length(y) != length(s2) | any(s2 < 0)) 
        stop("Errors in the input data.")
    if (abs(adjust_font_size) > 1) 
        stop("The adjusted font size cannot be unreasonable")
    if (missing(out)) {
        out <- metagen(y, sqrt(s2), sm = effect_type, method.tau = "ML")
    }
    
    I <- length(y)
    
    eff.size <- exp(y)
    eff.upper <- exp(out$upper)
    eff.lower <- exp(out$lower)
    
    tradeoff.lambda <- metapen.lamb(y, s2, n.lambda = 100, lam.c = 1.2)
    lower.lambda <- exp(tradeoff.lambda$mu.opt - 1.96 * (tradeoff.lambda$se.opt))
    upper.lambda <- exp(tradeoff.lambda$mu.opt + 1.96 * (tradeoff.lambda$se.opt))
    
    tradeoff.tau <- metapen.tau(y, s2)
    lower.tau <- exp(tradeoff.tau$mu.opt - 1.96 * (tradeoff.tau$se.opt))
    upper.tau <- exp(tradeoff.tau$mu.opt + 1.96 * (tradeoff.tau$se.opt))
    
    gap <- as.integer(-(I + 10))
    
    units <- structure(list(mean = c(NA, NA, eff.size, NA, exp(out$TE.fixed), NA, exp(out$TE.random), 
        NA, exp(tradeoff.lambda$mu.opt), NA, exp(tradeoff.tau$mu.opt)), lower = c(NA, NA, eff.lower, 
        NA, exp(out$lower.fixed), NA, exp(out$lower.random), NA, lower.lambda, NA, lower.tau), upper = c(NA, 
        NA, eff.upper, NA, exp(out$upper.fixed), NA, exp(out$upper.random), NA, upper.lambda, NA, 
        upper.tau)), .Names = c("mean", "lower", "upper"), row.names = c(NA, gap), class = "data.frame")
    
    eff.round <- format(round(eff.size, digits = 2), nsmall = 2)
    fix.round <- format(round(exp(out$TE.fixed), digits = 2), nsmall = 2)
    ran.round <- format(round(exp(out$TE.random), digits = 2), nsmall = 2)
    trdlam.round <- format(round(exp(tradeoff.lambda$mu.opt), digits = 2), nsmall = 2)
    trdtau.round <- format(round(exp(tradeoff.tau$mu.opt), digits = 2), nsmall = 2)
    upper.round <- format(round(eff.upper, digits = 2), nsmall = 2)
    lower.round <- format(round(eff.lower, digits = 2), nsmall = 2)
    fix.upper.round <- format(round(exp(out$upper.fixed), digits = 2), nsmall = 2)
    fix.lower.round <- format(round(exp(out$lower.fixed), digits = 2), nsmall = 2)
    ran.upper.round <- format(round(exp(out$upper.random), digits = 2), nsmall = 2)
    ran.lower.round <- format(round(exp(out$lower.random), digits = 2), nsmall = 2)
    trdlam.upper.round <- format(round(upper.lambda, digits = 2), nsmall = 2)
    trdlam.lower.round <- format(round(lower.lambda, digits = 2), nsmall = 2)
    trdtau.upper.round <- format(round(upper.tau, digits = 2), nsmall = 2)
    trdtau.lower.round <- format(round(lower.tau, digits = 2), nsmall = 2)
    
    CI <- paste("[", lower.round, ",", " ", upper.round, "]", sep = "")
    CI.fixed <- paste("[", fix.lower.round, ",", " ", fix.upper.round, "]", sep = "")
    CI.random <- paste("[", ran.lower.round, ",", " ", ran.upper.round, "]", sep = "")
    CI.tradeoff.lambda <- paste("[", trdlam.lower.round, ",", " ", trdlam.upper.round, "]", sep = "")
    CI.tradeoff.tau <- paste("[", trdtau.lower.round, ",", " ", trdtau.upper.round, "]", sep = "")
    
    tabletext <- cbind(c("Study", NA, as.character(out$studlab), NA, "CE", NA, "RE", NA, paste("PRE", 
        "(", intToUtf8(955), ")", sep = ""), NA, paste("PRE", "(", intToUtf8(964), ")", sep = "")), 
        c(effect_type, NA, eff.round, NA, fix.round, NA, ran.round, NA, trdlam.round, NA, trdtau.round), 
        c("95% CI", NA, CI, NA, CI.fixed, NA, CI.random, NA, CI.tradeoff.lambda, NA, CI.tradeoff.tau))
    
    space <- I + 2
    
    if (effect_type == "OR") {
        x_label <- "Odds ratio (log scale)"
    }
    if (effect_type == "RR") {
        x_label <- "Risk ratio (log scale)"
    }
    
    if (user_define_x_ticks == FALSE) {
        forestplot(tabletext, hrzl_lines = list(`2` = gpar(lty = 1)), units, new_page = F, is.summary = c(TRUE, 
            rep(FALSE, space), TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE), txt_gp = fpTxtGp(label = gpar(fontfamily = "", 
            cex = 0.6 + adjust_font_size), ticks = gpar(fontfamily = "", cex = 0.65 + adjust_font_size), 
            xlab = gpar(fontfamily = "", cex = 0.7 + adjust_font_size)), xlog = T, clip = c(llimit, 
            ulimit), colgap = unit(3, "mm"), line.margin = 0.35, col = fpColors(box = "royalblue", 
            lines = "darkblue", summary = "royalblue"), xlab = x_label, mar = unit(c(0.01, 3, 2, 3.5), 
            "mm"))
    } else {
        forestplot(tabletext, hrzl_lines = list(`2` = gpar(lty = 1)), units, new_page = F, is.summary = c(TRUE, 
            rep(FALSE, space), TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE), txt_gp = fpTxtGp(label = gpar(fontfamily = "", 
            cex = 0.6 + adjust_font_size), ticks = gpar(fontfamily = "", cex = 0.65 + adjust_font_size), 
            xlab = gpar(fontfamily = "", cex = 0.7 + adjust_font_size)), xlog = T, xticks = x_ticks, 
            clip = c(llimit, ulimit), colgap = unit(3, "mm"), line.margin = 0.35, col = fpColors(box = "royalblue", 
                lines = "darkblue", summary = "royalblue"), xlab = x_label, mar = unit(c(0.01, 3, 
                2, 3.5), "mm"))
    }
    
}




