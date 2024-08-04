# main script

# Import libraries
library(tidyverse)
library(here)
library(MASS)
library(glmnet)
library(dirmult)

# Reset environment
set.seed(1)
rm(list = ls())

# Import data sets
demog <- read_csv(here("data/demog.csv"))
sales <- read_csv(here("data/sales_visit.csv"))

demog <- demog |> dplyr::select(-...1)
sales <- sales |> dplyr::select(-...1)


# Descriptive statistics
ggplot(data = sales) +
  geom_point(mapping = aes(x = visit_num, y = sales_num))
ggsave(here("output/scatter1.png"))

df <- sales |> 
  group_by(id) |> 
  summarize(
    mean_sales = mean(sales_num),
    mean_visits = mean(visit_num)
  ) |> 
  dplyr::select(id, mean_sales, mean_visits) |> 
  right_join(demog, by = "id")
ggplot(data = df) +
  geom_point(mapping = aes(x = mean_visits, y = mean_sales))
ggsave(here("output/scatter2.png"))

# Define X and y
X <- cbind(1, df$mean_visits, df$exp10, df$avg_samp)
y <- df$mean_sales

# Regression
# Predetermined
K <- 3
p <- 3

# Settings
n = nrow(X)

# parameters of prior distributions
# beta_k ~ N(b0, B0)
b0 <- matrix(0, nrow=p+1)
B0 <- 100*diag(p+1)
# sigma2_k ~ IG(nu0/2, nu0*sigma0_2/2)
nu0 <- 0.1
sigma0_2 <- 0.1
alpha <- matrix(1, ncol=K)

# storage for MCMC draw
betas_MCMC <- NULL 
sigma2s_MCMC <- NULL
pi_MCMC <- NULL

# starting values
sigma2s <- matrix(10, nrow=1, ncol=K)
beta_OLS <- solve(t(X)%*%X)%*%t(X)%*%y
betas <- matrix(beta_OLS, nrow=p+1, ncol=K)
pi <- matrix(1,nrow=1,ncol=K)/K

## Start MCMC ##
n_MCMC <- 5000 # number of MCMC chain
n_bin <- 2000 #number of burn-in period

for (mc in 1:n_MCMC) {
  
  # 1. sampling z_i
  pp <- matrix(0,n,K)
  for (k in 1:K){
    pp_k <- dnorm(y, X%*%as.matrix(betas[,k]), sqrt(sigma2s[k]) )
    pp[,k] <- pi[k]*pp_k
  }
  p_z_i <- pp / apply(pp, 1, sum)
  z_i <- rmult(p_z_i)
  
  z_i_mat <- matrix(0,n,K)
  for (i in 1:n){
    z_i_mat[i,z_i[i]] <- 1
  } # matrix of latent class index
  N_k <- apply(z_i_mat, 2, sum) # number of class components
  
  # 2. sampling pi
  pi <- rdirichlet(1, N_k+alpha)
  
  # 3. sampling beta 
  for (k in 1:K){
    y_k <- y[ which(z_i_mat[,k] == 1)]
    X_k <- X[ which(z_i_mat[,k] == 1),]
    B_k <- solve((t(X_k)%*%(X_k))/sigma2s[k] + solve(B0))
    b_k <- B_k %*% ( (t(X_k)%*%y_k)/sigma2s[k]+ solve(B0)%*%b0) 
    beta_k <- as.matrix( mvrnorm(1, b_k, B_k))
    betas[,k] <- beta_k
  }
  
  # 4. sampling sigma2
  for (k in 1:K){
    y_k <- y[ which(z_i_mat[,k] == 1)]
    X_k <- X[ which(z_i_mat[,k] == 1),]
    nu1_k <- nu0 + N_k[k]
    s1_k <- nu0*sigma0_2 + t(y_k-X_k%*%betas[,k]) %*% (y_k-X_k%*%betas[,k])
    sigma2_k <- 1/(rgamma(1, nu1_k/2, s1_k/2))
    sigma2s[k] <- sigma2_k
  }
  
  
  # store the sample 
  betas_MCMC <- rbind(betas_MCMC, c(betas))
  sigma2s_MCMC <- rbind(sigma2s_MCMC, c(sigma2s) ) 
  pi_MCMC <- rbind(pi_MCMC,c(pi))
  
  print(mc)
}

# result table
poseterior.mean <- apply(betas_MCMC[(n_bin+1):n_MCMC,], 2, mean)  # posterior mean of betas
poseterior.sd <- apply(betas_MCMC[(n_bin+1):n_MCMC,], 2, sd) # posterior s.d. of betas
CI.lower <- apply(betas_MCMC[(n_bin+1):n_MCMC,], 2, quantile, probs=0.05)
CI.upper <- apply(betas_MCMC[(n_bin+1):n_MCMC,], 2, quantile, probs=0.95)
print(res_table <- as.matrix(t(rbind(poseterior.mean, poseterior.sd, CI.lower, CI.upper))) )
print(poseterior.mean_pi <- apply(pi_MCMC[(n_bin+1):n_MCMC,], 2, mean) )

# plot
min_pi_MCMC <- min(pi_MCMC)
max_pi_MCMC <- max(pi_MCMC)
ggplot() +
  labs(
    x = "MCMC iteration",
    y = "Simulated pi"
  ) +
  geom_line(aes(x = 1:n_MCMC, y = pi_MCMC[,1], color = "pi1")) +
  geom_line(aes(x = 1:n_MCMC, y = pi_MCMC[,2], color = "pi2")) +
  geom_line(aes(x = 1:n_MCMC, y = pi_MCMC[,3], color = "pi3"))
ggsave(here("output", "pi_MCMC.png"))

ggplot() +
  labs(
    x = "MCMC iteration",
    y = "Simulated constant"
  ) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,1], color = "Class 1")) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,5], color = "Class 2")) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,9], color = "Class 3"))
ggsave(here("output", "constant_MCMC.png"))

ggplot() +
  labs(
    x = "MCMC iteration",
    y = "Simulated beta 1"
  ) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,2], color = "Class 1")) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,6], color = "Class 2")) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,10], color = "Class 3"))
ggsave(here("output", "beta1_MCMC.png"))

ggplot() +
  labs(
    x = "MCMC iteration",
    y = "Simulated beta 2"
  ) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,3], color = "Class 1")) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,7], color = "Class 2")) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,11], color = "Class 3"))
ggsave(here("output", "beta2_MCMC.png"))

ggplot() +
  labs(
    x = "MCMC iteration",
    y = "Simulated beta 3"
  ) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,4], color = "Class 1")) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,8], color = "Class 2")) +
  geom_line(aes(x = 1:n_MCMC, y = betas_MCMC[,12], color = "Class 3"))
ggsave(here("output", "beta3_MCMC.png"))