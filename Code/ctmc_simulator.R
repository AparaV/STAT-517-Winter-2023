# Helper functions for Lab 2
# Author: Apara Venkat

simulate_ctmc <- function(X0, R, time) {
    state_space <- nrow(R)
    Xt <- c()
    t <- c()
    
    clock <- 0
    X_prev <- X0
    sojourn_time <- 0
    while (clock < time) {
        Xt <- c(Xt, X_prev)
        t <- c(t, sojourn_time)
        
        # Simulate the time it takes to move from X_prev to all adjacent states
        transition_rates <- R[X_prev, ]
        sojourn_times <- c()
        for (i in seq_along(transition_rates)) {
            q_i <- transition_rates[i]
            if (q_i > 0) {
                sojourn_i <- rexp(1, q_i)
            } else {
                sojourn_i <- Inf
            }
            sojourn_times <- c(sojourn_times, sojourn_i)
        }
        
        # The next state will be adjacent state we first jump to
        sojourn_time <- min(sojourn_times)
        X_prev <- which.min(sojourn_times)
        
        clock <- clock + sojourn_time
    }
    
    return(list(states=Xt, times=t))
}


estimate_transition_matrix <- function(R, t, k) {
    I <- diag(nrow=nrow(R))
    
    n <- 2^k
    Pt <- I + R*t / n 
    for (i in seq_len(k-1)) {
        Pt <- Pt %*% Pt
    }
    return(Pt)
}