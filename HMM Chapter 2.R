# Zucchini Chapter 2
library(tidyverse)

#### Problem 1 ####

# Consider STATIONARY 2-State Poisson-HMM
G <- matrix(c(0.1, 0.9, 0.4, 0.6), nrow=2, ncol=2, byrow = TRUE)
lambda <- c(1, 3)

n = 3
m = nrow(G)
ones <- rep(1, m)
I <- diag(ones)
U <- matrix(rep(1, length(G)), nrow=2, ncol=2)

# Calculate Delta (see Chapter 1 p18)
delta <- c(1,1) %*% solve(I-G+U)

# In each of the following ways, compute the probability that the first
# three observations from this model are 0, 2, 1.

# (a) Consider all possible sequences of states of the Markov chain that
# could have occurred. Compute the probability of each sequence,
# and the probability of the observations given each sequence.

data <- expand.grid(rep(list(1:2), 3)) %>% as_tibble() %>% 
  rename(i = Var1, j = Var2, k = Var3) %>% 
  mutate(`p_i(0)` = dpois(0, lambda = lambda[i]),
         `p_j(2)` = dpois(2, lambda = lambda[j]),
         `p_k(1)` = dpois(1, lambda = lambda[k]),
         delta_i = delta[i],
         gamma_ij = 0,
         gamma_jk = 0)

for(n in 1:nrow(data)){
  data$gamma_ij[n] = G[[data$i[[n]], data$j[[n]]]]
  data$gamma_jk[n] = G[[data$j[[n]], data$k[[n]]]]
}

data <- data %>% 
  mutate(Product = `p_i(0)` * 
                   `p_j(2)` *
                   `p_k(1)` * 
                    delta_i * 
                   gamma_ij * 
                   gamma_jk)

# P(X_1 = 0, X_2 = 2, X_3 = 1)
answerA <- sum(data$Product)
print(paste(answerA*100, "%", sep=""))

# (b) Apply the formula
# Pr(X1 = 0;X2 = 2;X3 = 1) = delta * P(0) * GAMMA * P(2) * GAMMA * P(1) * 1'
# P is diagonal matrix of dpois
P <- function(x, l){
  # x is data vale, l is vector of lambdas
  a <- numeric()
  
  for(n in 1:length(l)){
    a <- c(a,   dpois(x, l[[n]]))
  }
  
  a <- diag(a, nrow = length(l))
  
  return(a)
}

answerB <- delta %*% P(0, lambda) %*% 
  G %*% P(2, lambda) %*% G %*% P(1, lambda) %*% rep(1, length(lambda))

answerA == answerB
# Hooray