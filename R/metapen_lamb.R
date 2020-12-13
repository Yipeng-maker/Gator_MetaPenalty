##' Conduct the penalization method by tuning lambda
##'
##' Many values are returned in the outcome list. n.study is the number of studies for a given meta-analysis; tau2.re is the estimated between-study variance of the random-effects model; I2 is the value of the I^2 statistic; mu.fe represents the overall effect size estimate of the common-effect model; se.fe is the estimate of the standard deviation of the common-effect model; mu.re is the overall effect size estimate of the random-effects model; se.re is the estimated standard deviation of the random-effects model; tau2 (tau) is a vector of the between-study variance (standard deviation) estimates correspond to candidate values of lambda; loss is a vector of losses for the loss function of lambda; tau.opt (tau2.opt) is the optimal between-study standard deviation (variance) using the penalization method by tuning lambda; mu.opt and se.opt are the overall effect size and the between-study standard deviation correspond to optimal lambda
##' @title Obtain all necessary results of the penalization method by tuning lambda
##' @param est.eff an observed effect size vector
##' @param est.var the corresponding within-study varaince vector
##' @param tau2.re the estimated between-study variance of the random-effects model, the default value is the ML estimate
##' @param penalty the penalty function used in the penalized likelihood, the default value is 'tau2'
##' @param lambda a vector of candidate values for lambda
##' @param n.lambda the number of candidate values of lambda, the default number is 100
##' @param lambda.scale the scale of lambda will display, 'log' represent the log scale, 'linear' represent the linear scale
##' @param tau2.ml the ML estimate of the between-study variance. If it is not specified, the function will automatically compute the value
##' @param tol the relative convergence tolerence in optimization
##' @param lam.c a scalar used to adjust the upper bound of candidate values for lambda
##' @return a list consist of the results of the fixed-effect model, the random-effects model and the penalization approach by tuning lambda
##' @author Yipeng
##' @export
metapen.lamb <- function(est.eff, est.var, tau2.re, penalty = "tau2", lambda, n.lambda, lambda.scale = "log", 
    tau2.ml, tol = 10^(-10), lam.c = 1.2) {
    if (length(est.eff) != length(est.var) | any(est.var < 0)) 
        stop("error in input data.")
    I <- length(est.eff)
    
    if (!is.element(lambda.scale, c("log", "linear"))) 
        stop("lambda.scale must be either 'linear' or 'log'.")
    
    if (!missing(lambda) & !missing(n.lambda)) {
        if (n.lambda != length(lambda)) 
            stop("mismatched lambda and n.lambda.")
    }
    
    if (!missing(lambda)) {
        n.lambda <- length(lambda)
    }
    
    ## If lam.c = 1, the range of lambda is between 0 and lambda_max
    if (missing(lambda)) {
        if (missing(n.lambda)) 
            n.lambda <- 100
        lambda <- determine.lambda(est.eff = est.eff, est.var = est.var, penalty = penalty, n.lambda = n.lambda, 
            lambda.scale = lambda.scale, tol = tol, lam.c = lam.c)
    }
    
    y <- est.eff
    s2 <- est.var
    
    w <- 1/s2
    y.bar <- sum(w * y)/sum(w)
    Q <- sum(w * (y - y.bar)^2)
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
    
    diff <- matrix(NA, n.lambda, I)
    for (i in 1:I) {
        y_i <- y[-i]
        s2_i <- s2[-i]
        w_i <- 1/s2_i
        tau2.re_i <- metaml(y_i, s2_i, tol)
        out <- metapen.lambda(est.eff = y_i, est.var = s2_i, lambda = lambda, penalty = penalty, tol = tol)
        mu_i <- out[, "mu"]
        tau2_i <- out[, "tau2"]
        mu.hat.var <- lapply(X = 1:n.lambda, FUN = function(lam) {
            sum((s2_i + tau2.re_i)/(s2_i + tau2_i[lam])^2)/(sum(1/(s2_i + tau2_i[lam])))^2
        })
        mu.hat.var <- unlist(mu.hat.var)
        diffi <- (y[i] - mu_i)^2/(s2[i] + tau2.re_i + mu.hat.var)
        diff[, i] <- diffi
    }
    
    ## With the full dataset, using metapen.lambda( ) to estimate tau^2 and mu for a given lambda
    if (missing(tau2.ml)) {
        full <- metapen.lambda(est.eff = y, est.var = s2, lambda = lambda, penalty = penalty, tol = tol)
    } else {
        full <- metapen.lambda(est.eff = y, est.var = s2, lambda = lambda, penalty = penalty, tau2.ml = tau2.ml, 
            tol = tol)
    }
    mu <- full[, "mu"]
    tau2 <- full[, "tau2"]
    tau <- sqrt(tau2)
    
    loss <- apply(diff, 1, mean)
    loss <- sqrt(loss)
    opt.idx <- which(loss == min(loss))
    opt.idx <- opt.idx[1]  ## If several cases have the same minimize loss, select the first case
    tau2.opt <- tau2[opt.idx]
    tau.opt <- tau[opt.idx]
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
    out$tau2 <- tau2
    out$tau <- tau
    out$loss <- loss
    out$tau2.opt <- tau2.opt
    out$tau.opt <- tau.opt  ## The estimate of the between-study standard deviation
    out$mu.opt <- mu.opt  ## The estimate of the overall effect size
    out$se.opt <- sqrt(var.opt)  ## The standard error of the overall effect size estimate
    return(out)
}
