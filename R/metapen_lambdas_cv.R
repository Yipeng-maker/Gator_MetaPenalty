##' This function is used to save values to visualize the cross-validation process of the penalization method by tuning lambda
##'
##' All augements are the same in metapen.lamb function
##' @title Obtain values to show relationships among different parameters in the penalization method by tuning lambda
##' @return values of different parameters of the penalization approach by tuning lambda
##' @param est.eff an observed effect size vector
##' @param est.var the corresponding within-study variance vector
##' @param penalty the penalty function used in the penalized likelihood, the default value is 'tau2'
##' @param lambda a vector of candidate values for lambda
##' @param n.lambda the number of elements in the lambda vector
##' @param lambda.scale the scale of lambda will display, 'log' represent the log scale, 'linear' represent the linear scale
##' @param tau2.ml the ML estimate of the between-study variance. If it is not specified, the function will automatically compute the value
##' @param tol the relative convergence tolerence in optimization
##' @param lam.c a scalar used to adjust the upper bound of candidate values for lambda
##' @author Yipeng
##' @export
metapen.lambdas.cv <- function(est.eff, est.var, penalty = "tau2", lambda, n.lambda, lambda.scale = "log", 
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
    
    ## With the full dataset, using metapen.lambda function again to find tau^2 and mu for a given
    ## lambda
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
    est <- cbind(lambda, loss, mu, tau2, tau)
    rownames(est) <- rep("", n.lambda)
    
    out <- list(est = est, penalty = penalty, n.lambda = n.lambda, lambda.scale = lambda.scale, tol = tol)
    return(out)
}

