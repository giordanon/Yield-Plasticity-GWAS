functions {
  array[] int permutation_rng(int N) {
    array[N] int y;
    for (n in 1 : N) {
      y[n] = n;
    }
    vector[N] theta = rep_vector(1.0 / N, N);
    for (n in 1 : size(y)) {
      int i = categorical_rng(theta);
      int temp = y[n];
      y[n] = y[i];
      y[i] = temp;
    }
    return y;
  }
}

data {
  // Training data
  int<lower=0> N;          // Number of observations on training data
  int<lower=0> K;          // Number of covariates
  matrix[N, K] X;          // Matrix of covariates
  vector[N] y;             // Response variable
  int<lower = 0, upper = N> N_test; // Number of observations on testing data
  real lambda;            // Regularization hyper parameter 
}

transformed data {
  int N_train = N - N_test;
  int permutation[N] = permutation_rng(N);
  matrix[N_train,K] X_train = X[permutation[1 : N_train]];
  vector[N_train] y_train = y[permutation[1 : N_train]];
  matrix[N_test,K] X_test = X[permutation[N_train + 1 : N]];
  vector[N_test] y_test = y[permutation[N_train + 1 : N]];
}

parameters {
  real beta0; // Intercept
  vector[K] B; // Covariates
  real<lower=0.001> alpha;  // dispersion of variance
}
transformed parameters {
  vector[N_train] mu;
  
  for (n in 1:N_train){
    mu[n] = exp(beta0 + X_train[n]*B);
  }
}

model {
  // Prior distributions
  alpha ~ gamma(7,1); // First moment of the Gamma distribution
  B ~ double_exponential(0, lambda); // modify lambda hyper parameter for shrinkage
  
  for (n in 1:N_train){
    y_train[n] ~ gamma(alpha, alpha ./mu[n]);
  }
}

generated quantities {
  vector[N_test] y_hat;
  vector[N_test] mu_new;

  vector[N_test] indiv_squared_errors;
  real <lower = 0> sum_of_squares;
  real <lower = 0> root_mean_squared_error;
  
  for (n in 1:N_test){
    mu_new[n] = exp(beta0 + X_test[n]*B);
    }
    
  for (n in 1:N_test)  {
    y_hat[n] = gamma_rng(alpha, alpha ./mu_new[n]);
    indiv_squared_errors[n] = (y_test[n] - y_hat[n])^2;
    }
  
  // sum of squares, data set size and scale specific
  sum_of_squares = sum(indiv_squared_errors);

  // divide by number of new / test points and sqrt for RMSE
  root_mean_squared_error = sqrt(sum_of_squares / N_test);
  
}