source("helper_functions.R")



### Swendsen-Wang sampler

# We will use `weak` function from this library to find connected components
library(pooh)

draw_random_cluster_model <- function(X, b, m, n) {
    Y <- list(from=c(), to=c())
    p <- 1 - exp(-b)
    for (i in seq_len(m)) {
        for (j in seq_len(n)) {
            vec_idx <- lat_to_vec(i, j, m, n)
            if (j < n && X[i, j] == X[i, j+1]) {
                U <- rbinom(1, 1, p)
                if (U == 1) {
                    Y$from <- c(Y$from, vec_idx)
                    Y$to <- c(Y$to, vec_idx + 1)
                }
            }
            if (i < m && X[i, j] == X[i+1, j]) {
                U <- rbinom(1, 1, p)
                if (U == 1) {
                    Y$from <- c(Y$from, vec_idx)
                    Y$to <- c(Y$to, vec_idx + n)
                }
            }
        }
    }
    return(Y)
}


find_connected_components <- function(Y, m, n) {
    C <- matrix(0, nrow=m, ncol=n)
    num_nodes <- m*n
    # See documentation for `weak`
    # https://www.rdocumentation.org/packages/pooh/versions/0.3-1/topics/weak
    C_vec <- weak(Y$from, Y$to, domain=c(1:num_nodes), markers=T)
    for (node in seq_along(C_vec)) {
        lat_idx <- vec_to_lat(node, m, n)
        C[lat_idx[1], lat_idx[2]] <- C_vec[node]
    }
    return(C)
}

color_connected_components <- function(C, m, n, h, K) {
    if (K == 2) {
        colors <- c(-1, 1)
    }
    else {
        colors <- seq(1, K)
    }
    X <- matrix(0, nrow=m, ncol=n)
    num_C <- max(C)
    for (comp_id in seq_len(num_C)) {
        p_idx <- which(C == comp_id, arr.ind=T)
        
        # If no external field, then sample uniformly
        if (h == 0) {
            color <- sample(colors, 1)
        }
        # Otherwise, it is slightly more complicated
        else {
            prob_color <- rep(0, K)
            C_size <- length(which(C == comp_id))
            for (i in seq_len(K)) {
                prob_color[i] <- exp(C_size * colors[i])
                if (is.infinite(prob_color[i])) {
                    prob_color[i] <- 1
                }
            }
            prob_color <- prob_color / sum(prob_color)
            color <- sample(colors, 1, prob=prob_color)
        }
        X[p_idx] <- color
    }
    return(X)
}

swendsen_wang <- function(X, b, h=0, K=2, iters=1000, burnin=100) {
    m <- nrow(X)
    n <- ncol(X)
    X_post <- c()
    for (i in seq_len(iters)) {
        Y <- draw_random_cluster_model(X, b, m, n)
        C <- find_connected_components(Y, m, n)
        X <- color_connected_components(C, m, n, h, K)
        
        if (i > burnin) {
            X_post <- c(X_post, list(X))
        }
    }
    
    return(X_post)
}


## Sampling phi for Question 5
estimate_density <- function(A, X, phi, log=T) {
    q <- phi * A %*% X
    if (log == FALSE) {
        q <- exp(q)
    }
    return(q)
}


sample_phi <- function(Z, K, phi_0, alpha_phi, beta_phi) {
    m <- nrow(Z)
    n <- ncol(Z)
    A <- lattice_to_adjacency(m, n)
    
    Z_vec <- c(t(Z))
    
    # Propose a new value
    phi_t <- rgamma(1, shape=phi_0, rate=1)
    
    # Conditional proposal densities
    h_t_0 <- dgamma(phi_t, shape=phi_0, rate=1, log=T)
    h_0_t <- dgamma(phi_0, shape=phi_t, rate=1, log=T)
    
    # Priors
    prior_phi_t <- dgamma(phi_t, alpha_phi, 1/beta_phi, log=T)
    prior_phi_0 <- dgamma(phi_0, alpha_phi, 1/beta_phi, log=T)
    
    # Density without normalization constant
    qz_phi_0 <- estimate_density(A, Z_vec, phi_0)
    qz_phi_t <- estimate_density(A, Z_vec, phi_t)
    
    # Estimate ratio of normalization constants
    # Set N > 1 for noisy exchange and N = 1 for exchange algorithm
    # Setting N = 1 for quick computation
    N <- 1
    ratio <- -log(N)
    for (i in seq_len(N)) {
        # Ideally niters and burnin should be larger. But that slows down
        # the Gibbs sampler.
        niters <- 20
        burnin <- 10
        Zi_samples <- swendsen_wang(Z, phi_t, 0, K, iters=niters, burnin=burnin)
        Zi <- Zi_samples[[niters-burnin]]
        Zi_vec <- c(t(Zi))
        
        q_zi_phi_0 <- estimate_density(A, Zi_vec, phi_0)
        q_zi_phi_t <- estimate_density(A, Zi_vec, phi_t)
        
        ratio <- ratio + q_zi_phi_0 - q_zi_phi_t
    }
    
    # Final step to sample phi_t
    alpha <- qz_phi_t + prior_phi_t + h_0_t
    alpha <- alpha - (qz_phi_0 + prior_phi_0 + h_t_0)
    alpha <- alpha + ratio
    acceptance_prob <- min(1, exp(alpha))
    U <- runif(1, 0, 1)
    if (U < acceptance_prob) {
        phi_0 <- phi_t
    }
    
    return(phi_0)
}
