##' Obtain candidate values of the tuning parameter lambda
##'
##' A threshold value of lambda is used to define the upper bound of candidate values for lambda, and the lower bound is zero
##' @title Obtain candidate values for the tuning parameter lambda
##' @param est.eff an observed effect size vector
##' @param est.var the corresponding within-study varaince vector
##' @param penalty the penalty function used in the penalized likelihood, the default value is 'tau2'
##' @param n.lambda the number of candidate values of lambda, the default number is 100
##' @param lambda.scale the transformed scale of lambda, the default value is 'log'
##' @param tol the relative convergence tolerence in optimization
##' @param lam.c a scalar used to mutiply lambda_max, by controlling its magnitude we can specify the upper bound of candidate values for lambda
##' @return a threshold value for tuning parameters lambda
##' @author Yipeng
##' @export
determine.lambda <- function(est.eff, est.var, penalty = "tau2", n.lambda = 100, lambda.scale = "log", 
    tol, lam.c = 1.2) {
    if (missing(tol)) 
        stop("tolerence in optimization must be specified")
    y <- est.eff
    s2 <- est.var
    
    mu.hat <- function(tau2) {
        sum(y/(s2 + tau2))/sum(1/(s2 + tau2))
    }
    mu.hat <- Vectorize(mu.hat)
    
    
    target.ml <- function(tau2) {
        mean(log(s2 + tau2) + (y - mu.hat(tau2))^2/(s2 + tau2))
    }
    target.ml <- Vectorize(target.ml)
    
    target1.ml <- function(tau2) {
        p1 <- mean(1/(s2 + tau2))
        p2 <- mean((y - mu.hat(tau2))^2/(s2 + tau2)^2)
        out <- p1 - p2
        return(out)
    }
    target1.ml <- Vectorize(target1.ml)
    
    
    tau2.ml <- optim(par = var(y), fn = target.ml, gr = target1.ml, method = "Brent", lower = 0, upper = 100 * 
        var(y), control = list(reltol = tol))$par
    
    if (penalty == "tau2") {
        fcn <- function(tau2) {
            sum((y - mu.hat(tau2))^2/(s2 + tau2)^2) - sum(1/(s2 + tau2))
        }
    }
    if (penalty == "tau") {
        fcn <- function(tau2) {
            (sum((y - mu.hat(tau2))^2/(s2 + tau2)^2) - sum(1/(s2 + tau2))) * (2 * sqrt(tau2))
        }
    }
    fcn <- Vectorize(fcn)
    lam.max <- optimize(f = fcn, lower = 0, upper = tau2.ml * 1.1, maximum = TRUE)$objective
    if (lambda.scale == "linear") 
        lam.max <- max(c(0, lam.max)) * lam.c
    if (lambda.scale == "log") 
        lam.max <- exp(log(max(c(0, lam.max)) + 1) + log(lam.c)) - lam.c
    
    if (lambda.scale == "linear") 
        lambda <- seq(from = 0, to = lam.max, length.out = n.lambda)
    if (lambda.scale == "log") 
        lambda <- exp(seq(from = 0, to = log(lam.max + 1), length.out = n.lambda)) - 1
    
    return(lambda)
}
