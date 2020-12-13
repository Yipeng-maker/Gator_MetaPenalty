##' This function is used to obtain the results of the penalization method by tuning tau
##'
##' Many values are returned as outcomes, which are similar as the results of the function metapen_lamb( ). Except two values, tau.cand represent the candidates of standard deviations used in the penalization method, where tau2.cand represent the corresponding candidates of between-study varainces
##' @title Obtain the results of the penalization method by tuning tau
##' @param y an observed effect size vector
##' @param s2 the corresponding within-study varaince vector
##' @param tau2.re the between-study variance estimate of the random-effects model, if it is not specified,
##' the ML estimate will be used as default
##' @param upp a scalar used to multiply tau2.re. Controlling its magnitude can adjust the range of tau
##' @param n.cand the number of candidates for the standard deviation considered in the penalization method
##' @param tol the relative convergence tolerence in optimization
##' @return values of different parameters of the penalization approach by tuning tau
##' @author Yipeng
##' @export
metapen.tau <- function(y, s2, tau2.re, upp, n.cand = 100, tol = 10^(-10)) {
    if (length(y) != length(s2) | any(s2 < 0)) 
        stop("Errors in the input data.")
    I <- length(y)
    
    w <- 1/s2
    y.bar <- sum(w * y)/sum(w)
    Q <- sum(w * (y - y.bar)^2)
    if (Q < I - 1) 
        Q <- I - 1
    I2 <- (Q - (I - 1))/Q
    
    if (missing(tau2.re)) {
        tau2.ml <- metaml(y, s2, tol)
        tau2.re <- tau2.ml
    }
    
    mu.fe <- y.bar
    se.fe <- sqrt(1/sum(w))
    w.re <- 1/(s2 + tau2.re)
    mu.re <- sum(w.re * y)/sum(w.re)
    se.re <- sqrt(1/sum(w.re))
    
    if (missing(upp)) 
        upp <- 1
    tau2.upp <- upp * tau2.re
    
    tau.cand <- seq(from = 0, to = sqrt(tau2.upp), length.out = n.cand)
    tau2.cand <- tau.cand^2
    
    errs <- matrix(0, I, n.cand)
    for (i in 1:I) {
        y.train <- y[-i]
        s2.train <- s2[-i]
        tau2.re.train <- metaml(y.train, s2.train)
        for (j in 1:n.cand) {
            tau2.temp <- tau2.cand[j]
            w.temp <- 1/(s2.train + tau2.temp)
            mu_i <- sum(w.temp * y.train)/sum(w.temp)
            var_i <- sum(w.temp^2 * (s2.train + tau2.re.train))/(sum(w.temp))^2
            errs[i, j] <- (mu_i - y[i])^2/(var_i + s2[i] + tau2.re.train)
        }
    }
    loss <- sqrt(colMeans(errs))
    opt.idx <- which(loss == min(loss))
    opt.idx <- opt.idx[1]  ## If several cases have the same minimize loss, select the first case
    tau2.opt <- tau2.cand[opt.idx]
    tau.opt <- tau.cand[opt.idx]
    w.opt <- 1/(s2 + tau2.opt)
    mu.opt <- sum(w.opt * y)/sum(w.opt)
    var.opt <- sum(w.opt^2 * (s2 + tau2.re))/(sum(w.opt))^2
    
    out <- NULL
    out$n.study <- I
    out$tau2.re <- tau2.re
    out$I2 <- I2
    out$mu.fe <- mu.fe
    out$se.fe <- se.fe
    out$mu.re <- mu.re
    out$se.re <- se.re
    out$tau2.cand <- tau2.cand
    out$tau.cand <- tau.cand
    out$loss <- loss
    out$tau2.opt <- tau2.opt
    out$tau.opt <- tau.opt  ## The estimate of the between-study standard deviation
    out$mu.opt <- mu.opt  ## The estimate of the overall effect size
    out$se.opt <- sqrt(var.opt)  ## The standard error of the overall effect size estimate
    return(out)
}
