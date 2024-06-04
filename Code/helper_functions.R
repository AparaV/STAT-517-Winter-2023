
# From lattice position, get a unique node ID
# Unrolls the matrix into a vector
lat_to_vec <- function(i, j, m, n) {
    k <- n*(i-1) + j
    return(k)
}

# Convert unique node ID to lattice position
# Roll a vector into a matrix
vec_to_lat <- function(k, m, n) {
    i <- floor((k-1) / n) + 1
    j <- (k %% n)
    if (j == 0) {
        j <- n
    }
    return(c(i, j))
}

# Get locations of neighbors of a node on regular finite lattice
get_neighbors_idx <- function(i, j, m, n) {
    neighbors <- c()
    if (j - 1 > 0) {
        neighbors <- c(neighbors, list(c(i, j-1)))
    }
    if (j + 1 <= n) {
        neighbors <- c(neighbors, list(c(i, j+1)))
    }
    if (i - 1 > 0) {
        neighbors <- c(neighbors, list(c(i-1, j)))
    }
    if (i + 1 <= m) {
        neighbors <- c(neighbors, list(c(i+1, j)))
    }
    return(neighbors)
}

# Function to obtain neighbors of a node in a lattice
get_neighbors <- function(X, i, j) {
    M <- nrow(X)
    N <- ncol(X)
    neighbors_idx <- get_neighbors_idx(i, j, M, N)
    neighbors <- c()
    for (k in seq_along(neighbors_idx)) {
        i1 <- neighbors_idx[[k]][1]
        j1 <- neighbors_idx[[k]][2]
        neighbors <- c(neighbors, X[i1, j1])
    }
    return(neighbors)
}


# Extract adjacency matrix of a regular 2D lattice
lattice_to_adjacency <- function(m, n) {
    num_nodes <- m*n
    A <- matrix(0, ncol=num_nodes, nrow=num_nodes)
    for (node in seq_len(num_nodes)) {
        lat_loc <- vec_to_lat(node, m, n)
        neighbors_idx <- get_neighbors_idx(lat_loc[1], lat_loc[2], m, n)
        for (k in seq_along(neighbors_idx)) {
            i <- neighbors_idx[[k]][1]
            j <- neighbors_idx[[k]][2]
            neighbor_loc <- lat_to_vec(i, j, m, n)
            A[node, neighbor_loc] <- A[neighbor_loc, node] <- 1
        }
    }
    return(A)
}



### Finding modes of lattices
find_mode <- function(v) {
    unique_val <- unique(v)
    return(unique_val[which.max(tabulate(match(v, unique_val)))])
}

find_marginal_mode <- function(X_posterior) {
    X <- X_posterior[[1]]
    m <- nrow(X)
    n <- ncol(X)
    locations <- vector("list", length = m*n)
    for (k in seq_along(X_posterior)) {
        X_k <- X_posterior[[k]]
        for (i in seq_len(m)) {
            for (j in seq_len(n)) {
                idx <- lat_to_vec(i, j, m, n)
                locations[[idx]] <- c(locations[[idx]], X_k[i, j])
            }
        }
    }
    for (k in seq_along(locations)) {
        lat_idx <- vec_to_lat(k, m, n)
        X[lat_idx[1], lat_idx[2]] <- find_mode(locations[[k]])
    }
    return(X)
}
