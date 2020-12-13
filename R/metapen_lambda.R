##' Compute the corresponding overall effect size estimates and the between-study variances for a vector of lambda
##'
##' For a vector of candidate values of lambda, calculate estimates for the overall effect size and the between-study varaince
##' @title Implement the penalization method by tuning lambda
##' @param est.eff an observed effect size vector
##' @param est.var the corresponding within-study variance vector
##' @param penalty the penalty function used in the penalized likelihood, the default value is 'tau2'
##' @param lambda a vector of candidate values for lambda
##' @param tau2.ml the ML estimate of the between-study variance
##' @param tol the relative convergence tolerence in optimization
##' @return a list consist of candidate values of lambda, estimates of the corresponding overall effect size and the between-study variance
##' @author Yipeng
##' @export
metapen.lambda <- function(est.eff, est.var, penalty = "tau2", lambda, tau2.ml, tol = 10^(-10)) {
    if (length(est.eff) != length(est.var) | any(est.var < 0)) 
        stop("error in input data.")
    I <- length(est.eff)
    
    y <- est.eff
    s2 <- est.var
    
    if (!is.element(penalty, c("tau", "tau2"))) 
        stop("penalty must be specified as 'tau' or 'tau2'.")
    if (penalty == "tau2") {
        pen <- function(tau2) {
            tau2
        }
        pen1 <- function(tau2) {
            1
        }
    } else {
        pen <- function(tau2) {
            sqrt(tau2)
        }
        pen1 <- function(tau2) {
            1/(2 * sqrt(tau2))
        }
    }
    pen <- Vectorize(pen)
    pen1 <- Vectorize(pen1)
    
    mu.hat <- function(tau2) {
        sum(y/(s2 + tau2))/sum(1/(s2 + tau2))
    }
    mu.hat <- Vectorize(mu.hat)
    
    if (missing(tau2.ml)) {
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
        
        tau2.ml <- optim(par = var(y), fn = target.ml, gr = target1.ml, method = "Brent", lower = 0, 
            upper = 100 * var(y), control = list(reltol = tol))$par
    }
    
    tau2.hat <- function(lambda) {
        target <- function(tau2) {
            mean(log(s2 + tau2) + (y - mu.hat(tau2))^2/(s2 + tau2)) + lambda * pen(tau2)/I
        }
        target <- Vectorize(target)
        
        target1 <- function(tau2) {
            p1 <- mean(1/(s2 + tau2))
            p2 <- mean((y - mu.hat(tau2))^2/(s2 + tau2)^2)
            out <- p1 - p2 + lambda * pen1(tau2)/I
            return(out)
        }
        
        out <- optim(par = tau2.ml, fn = target, gr = target1, method = "Brent", lower = 0, upper = tau2.ml * 
            1.1, control = list(reltol = tol))$par
        return(out)
    }
    tau2.hat <- Vectorize(tau2.hat)
    
    out.tau2 <- tau2.hat(lambda)
    out.mu <- mu.hat(out.tau2)
    out <- cbind(lambda = lambda, mu = out.mu, tau2 = out.tau2)
    rownames(out) <- rep("", length(lambda))
    return(out)
}
