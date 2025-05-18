//
  // Learn more about model development with Stan at:
  //
  //    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

// The input data
data {
  int<lower=2> n;
  int<lower=2> p;
  matrix[n,p] logZ;
  cov_matrix[n] C;
  real mean_prior_logY;
  real<lower=0> sd_prior_logY;
}

transformed data {
  vector[n] mu_logY;
  real var_Y;
  matrix[p,p] Id_p;
  matrix[n,n] Id_n;
  matrix[n,n] C0;
  mu_logY = rep_vector(mean_prior_logY, n);
  var_Y = sd_prior_logY^2;
  Id_p = diag_matrix(rep_vector(1.0, p));
  Id_n = diag_matrix(rep_vector(1.0, n));
  C0 = C;
  for (i in 1:n)
  {
    C0[i,i] = 0;
  }
}

// The parameters accepted by the model.
parameters {
  simplex[p] Z0;
  vector[n] logY_trans;
  cov_matrix[p] R;
  real<lower=0, upper=1> lambda;
}

transformed parameters {
  cov_matrix[n] Clambda;
  cholesky_factor_cov[n] cholfact_logY;
  vector[n] logY;
  Clambda = C0*lambda + Id_n;
  cholfact_logY = cholesky_decompose(var_Y*Clambda);
  logY = mu_logY + cholfact_logY*logY_trans;
}


// The model to be estimated.
model {
  matrix[n,p] logXdiff;
  for (i in 1:n)
  {
    logXdiff[i] = logZ[i] + logY[i] - log(Z0)';
  }
  target += -0.5*n*log_determinant(R) - 0.5*p*log_determinant(Clambda) - 0.5*trace(mdivide_left_spd(R, logXdiff') * mdivide_left_spd(Clambda, logXdiff));
  logY_trans ~ std_normal(); 
  R ~ inv_wishart(p, Id_p);
}
