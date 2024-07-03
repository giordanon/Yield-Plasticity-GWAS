data {
  // Training data
  int<lower=0> N;          // Number of observations
  int<lower=0> K;          // Number of covariates
  matrix[N, K] X;          // Matrix of covariates
  vector[N] y;             // Response variable
  cholesky_factor_cov[N] L;
  real lambda;            // Regularization hyper parameter 
}

parameters {
  //real beta0; // Intercept
  vector[K] B; // Covariates
  real<lower=0.001> alpha;  // dispersion of variance
  vector<lower=0.001>[N] z; // random effect primitive
  real<lower=1e-4> sigma_ran;
  
}

transformed parameters {
  vector[N] u;
  vector[N] mu;
  u = L * z;
  for (n in 1:N){
    mu[n] = exp(u[n] + X[n]*B);
  }
}

model {
  // Prior distributions
  alpha ~ gamma(7,1); // First moment of the Gamma distribution
  B ~ double_exponential(0, lambda); // modify lambda hyper parameter for shrinkage
  sigma_ran ~ inv_gamma(2, 1); // prior for the variance of the random effect
  z ~ normal(0, sigma_ran);
  
  for (n in 1:N){
    y[n] ~ gamma(alpha, alpha ./mu[n]);
  }
}

generated quantities {
  vector[N] Y_new;
  
  vector[N] indiv_squared_errors;
  real <lower = 0> sum_of_squares;
  real <lower = 0> root_mean_squared_error;
 
    for (n in 1:N)  {
      Y_new[n] = gamma_rng(alpha, alpha ./mu[n]);
      indiv_squared_errors[n] = (y[n] - Y_new[n])^2;
    }
  
  // sum of squares, data set size and scale specific
  sum_of_squares = sum(indiv_squared_errors);

  // divide by number of new / test points and sqrt for RMSE
  root_mean_squared_error = sqrt(sum_of_squares / N);
  
}