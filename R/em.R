em <- function(x, k = 2, tol = 1e-6, max_iter = 1000) {

    if (identical(k,'auto')){
        bic <- find_optimal_k(x, tol = tol, max_iter = max_iter)
        k <- bic$best_k
    }
    
    n <- length(x)
    
    mu <- seq(from=min(x),to=max(x),length.out=k)
    sigma <- rep(stats::var(x), k)
    prop <- rep(1/k, k)
  
    log_likelihood <- 0
    log_likelihoods <- numeric(max_iter)
    
    for (iter in 1:max_iter) {
        
        # E-step: calculate responsibilities
        gam <- matrix(0, n, k)
        for (j in 1:k) {
            gam[, j] <- prop[j] * stats::dnorm(x, mean = mu[j], sd = sqrt(sigma[j]))
        }
        gam <- gam / rowSums(gam)
    
        # M-step: update parameters
        n_k <- colSums(gam)
        for (j in 1:k) {
            mu[j] <- sum(gam[, j] * x) / n_k[j]
            sigma[j] <- sum(gam[, j] * (x - mu[j])^2) / n_k[j]
            prop[j] <- n_k[j] / n
        }
    
        # Check for convergence
        LL <- sapply(1:k,
                     function(j)
                         prop[j] * stats::dnorm(x, mean = mu[j], sd = sqrt(sigma[j]))
                     )
        new_log_likelihood <- sum(log(rowSums(LL)))
        log_likelihoods[iter] <- new_log_likelihood

        if (abs(new_log_likelihood - log_likelihood) < tol) {
            log_likelihood <- new_log_likelihood
            break
        }

        log_likelihood <- new_log_likelihood
    }
    
    list(mu = mu, sigma = sigma, prop = prop,
         log_likelihood = log_likelihood,
         log_likelihoods = log_likelihoods[1:iter])
}

calculate_bic <- function(log_likelihood, n, k) {
    num_params <- k * 3 - 1
    bic <- -2 * log_likelihood + num_params * log(n)
    return(bic)
}

find_optimal_k <- function(x, max_k = 5, tol = 1e-6, max_iter = 1000) {
    n <- length(x)
    best_k <- 1
    best_bic <- Inf
    bic_values <- numeric(max_k)
    
    for (k in 1:max_k) {
        model <- em(x, k = k, tol = tol, max_iter = max_iter)
        
        bic <- calculate_bic(model$log_likelihood, n, k)
        bic_values[k] <- bic
        
        if (bic < best_bic) {
            best_bic <- bic
            best_k <- k
        } else {
            break
        }
    }
    
    list(best_k = best_k, best_bic = best_bic, bic_values = bic_values)
}
