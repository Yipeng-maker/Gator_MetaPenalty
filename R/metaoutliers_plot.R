##' Plot study-specific standardized residuals for a given meta-analysis to check potential outliers
##'
##' @title A study-specific standardized residuals plot of a meta-analysis
##' @param est.eff an observed effect size vector
##' @param est.var the corresponding within-study varaince vector
##' @author Yipeng
##' @export
metaoutliers.plot <- function(est.eff, est.var) {
    if (length(est.eff) != length(est.var)) 
        stop("error in input data.")
    I <- length(est.eff)
    
    y <- est.eff
    s2 <- est.var
    
    ## Compute the standardized residuals under the common-effect setting
    w <- 1/s2
    mu.hat.i <- e <- e.tilde.fe <- numeric(I)
    for (i in 1:I) {
        w.temp <- w[-i]
        y.temp <- y[-i]
        mu.hat.i[i] <- sum(y.temp * w.temp)/sum(w.temp)
        e[i] <- y[i] - mu.hat.i[i]
        sig2.e.i <- 1/sum(w.temp) + 1/w[i]
        e.tilde.fe[i] <- e[i]/sqrt(sig2.e.i)
    }
    
    ## Compute the standardized residuals under the random-effects setting
    mu.hat.i <- e <- e.tilde.re <- numeric(I)
    for (i in 1:I) {
        s2.temp <- s2[-i]
        y.temp <- y[-i]
        tau2.temp <- metapen.tau(y.temp, s2.temp)$tau2.re  ##
        w.temp <- 1/(s2.temp + tau2.temp)
        mu.hat.i[i] <- sum(y.temp * w.temp)/sum(w.temp)
        e[i] <- y[i] - mu.hat.i[i]
        var.e.i <- 1/sum(w.temp) + s2[i] + tau2.temp
        e.tilde.re[i] <- e[i]/sqrt(var.e.i)
    }
    
    std.res <- c(e.tilde.fe, e.tilde.re)
    
    shapes <- NULL
    for (i in 1:I) {
        if (std.res[i] <= -5) {
            shapes[i] = 3
            std.res[i] = -5
        } else if (std.res[i] > -5 & std.res[i] < 5) {
            shapes[i] = 17
        } else {
            shapes[i] = 3
            std.res[i] = 5
        }
    }
    
    for (i in (I + 1):(2 * I)) {
        if (std.res[i] <= -5) {
            shapes[i] = 4
            std.res[i] = -5
        } else if (std.res[i] > -5 & std.res[i] < 5) {
            shapes[i] = 1
        } else {
            shapes[i] = 4
            std.res[i] = 5
        }
    }
    
    par(mar = c(2.4, 2.8, 1, 1) + 0.1)
    plot(std.res[1:I], I:1, ylim = c(1, I), xlim = c(-5, 5), yaxt = "n", xaxt = "n", pch = shapes[1:I], 
        cex = 0.6)
    title(xlab = "Standardized residuals", ylab = "Study", cex.lab = 0.8, line = 1.1)
    axis(side = 2, at = 1:I, labels = I:1, tck = -0.03, cex.axis = 0.6, las = 1, mgp = c(3, 0.5, 0))
    axis(side = 1, at = -5:5, tck = -0.03, cex.axis = 0.6, mgp = c(3, 0.2, 0))
    abline(v = c(-3, 3), lwd = 1, lty = 2)
    abline(v = 0, lwd = 1, lty = 3, col = "dimgrey")
    
    points(std.res[(I + 1):(2 * I)], I:1, pch = shapes[(I + 1):(2 * I)], cex = 0.9)
    
}

