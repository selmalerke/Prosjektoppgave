functions {
  // mean function
  // vector mu_fn(real mu_wk2, 
  //              real R, 
  //              int n_theta, 
  //              int nI){
  //                real R01 = R*2.5 + 0.5;
  //                vector[n_theta]  mP = rep_vector(mu_wk2,n_theta); 
  //                vector[nI]  mI = rep_vector((1/R01)*mu_wk2,nI);
  //                vector[n_theta + nI] mu;
  //                mu = append_row(mP, mI);
  //                return(mu);
  // }
  // Physics informed prior kernel of the WK2 model
  // Here are two funcs, I only have one
  matrix K_pendulum(matrix t_theta, 
               matrix t_E,
               real l,
               real sigma,
               real sigma_theta,
               real sigma_E,
               real R) {
    int n_theta = rows(t_theta);
    int nI = rows(t_E);
    matrix[n_theta + nI, n_theta + nI] K;
    // real R01 = R*2.5 + 0.5;
    //real C01 = C*2.5 + 0.5;
    // K_uu
    for (i in 1:(n_theta-1)){
      K[i,i] =   pow(sigma, 0.2e1);
      for (j in (i+1):n_theta){
        K[i,j] = exp(-0.5e0 * pow(t_theta[i,1] - t_theta[j,1], 0.2e1) * pow(l, -0.2e1));
        K[i,j] = pow(sigma, 0.2e1) * K[i,j];
        K[j,i] = K[i,j];
      }
      K[n_theta,n_theta] = pow(sigma, 0.2e1);
    }
    K[1:n_theta, 1:n_theta] = K[1:n_theta, 1:n_theta] + diag_matrix(rep_vector(pow(sigma_theta, 0.2e1), n_theta)); // observasjonsstøy
    
    // K_uf
    for (i in 1:n_theta){
      for (j in 1:nI){
        K[i, n_theta + j] =-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) + sigma * sigma * exp(-0.5e0 * pow((t_theta[i, 1] - t_E[j, 1]), 0.2e1) * pow(l, -0.2e1)) / R;

      }
    }
    
    
    // K_fu 
    
    for (i in 1:nI){
      for (j in 1:n_theta){
       K[n_theta + i, j]  =-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) + sigma * sigma * exp(-0.5e0 * pow((t_E[i, 1] - t_theta[j, 1]), 0.2e1) * pow(l, -0.2e1)) / R;

       // K[n_theta + i, j] = pow(sigma, 0.2e1) * K[n_theta + i, j];
     }
   }

   
   // K_ff
    for (i in 1:(nI-1)){
     K[n_theta + i, n_theta +i] =0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.10000e1 * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 / R * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[i, 1] - t_E[i, 1], 0.2e1) * pow(l, -0.2e1)));

     
     // K[n_theta + i, n_theta +i] = pow(sigma, 0.2e1) * K[n_theta + i, n_theta +i];
     // if(yI[i]!=yDia) K[n_theta + i, n_theta +i] = K[n_theta + i, n_theta +i] + pow(sigma_E,0.2e1);
     // if(yI[i]==yDia) K[n_theta + i, n_theta +i] = K[n_theta + i, n_theta +i] + 1e-12;
     for (j in (i+1):nI){
      K[n_theta + i, n_theta +j] =0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[i,1] - t_E[j,1], 0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.10000e1 * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 / R * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[i,1] - t_E[j,1], 0.2e1) * pow(l, -0.2e1)));

      // K[n_theta + i, n_theta +j] = pow(sigma, 0.2e1) * K[n_theta + i, n_theta +j];
      K[n_theta + j, n_theta +i] = K[n_theta + i, n_theta +j];
     }
     K[n_theta + nI, n_theta +nI] = 0.300e1 * sigma * sigma * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.2e1)) - 0.6000e1 * sigma * sigma * pow(l, -0.6e1) * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * exp(-0.5e0 * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.2e1)) + 0.10000e1 * sigma * sigma * pow(t_E[nI,1] - t_E[nI,1], 0.4e1) * pow(l, -0.8e1) * exp(-0.5e0 * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.2e1)) - 0.10e1 / R * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 / R * sigma * sigma * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * (-0.10e1 * sigma * sigma * pow(l, -0.2e1) * exp(-0.5e0 * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.2e1)) + 0.100e1 * sigma * sigma * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.4e1) * exp(-0.5e0 * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.2e1)) + 0.1e1 / R * sigma * sigma * exp(-0.5e0 * pow(t_E[nI,1] - t_E[nI,1], 0.2e1) * pow(l, -0.2e1)));

     //K[n_theta + nI, n_theta +nI] = pow(sigma, 0.2e1) * K[n_theta + nI, n_theta +nI];
     // if(yI[nI]!=yDia) K[n_theta + nI, n_theta +nI] = K[n_theta + nI, n_theta +nI] + pow(sigma_E,0.2e1);
     // if(yI[nI]==yDia) K[n_theta + nI, n_theta +nI] = K[n_theta + nI, n_theta +nI] + 1e-12;
     }
     K[(n_theta + 1):(n_theta + nI), (n_theta + 1):(n_theta + nI)] = K[(n_theta + 1):(n_theta + nI), (n_theta + 1):(n_theta + nI)]
            + diag_matrix(rep_vector(pow(sigma_E, 0.2e1), nI));
     return cholesky_decompose(K);
  }
}


data {
  int<lower=1> n_theta;
  int<lower=1> nI;
  // int<lower=1> n_theta_pred;
  // int<lower=1> nI_pred;
  matrix[n_theta,1] t_theta;
  matrix[nI,1] t_E;
  // matrix[n_theta_pred,1] t_theta_pred;
  // matrix[nI_pred,1] t_E_pred;
  vector[n_theta] yP;
  vector[nI] yI;
  real <lower=0> prior_sigma_theta_mu;
  real<lower=0>prior_sigma_theta_var;
  real<lower=0>range_sigma_theta_upper;
}

transformed data {
  vector[n_theta + nI] y = append_row(yP, yI);
  real sigma_E= 1e-3;
  // real sigma_theta= 1e-2;
}

parameters {
  // hyper-parameters
  real<lower=0.001> l;
  real<lower=0.001> sigma;
  // real<lower=0> mu_wk2;
  real<lower=0, upper= range_sigma_theta_upper> sigma_theta;
  //real<lower=0> sigma_E;
  // physical parameters
  real<lower=0> R;
  // real<lower=0, upper=1> C;
}

model {
  // Chol. of PI kernel
  matrix[n_theta + nI, n_theta + nI] L_K = K_pendulum(t_theta, t_E, l, sigma, sigma_theta, sigma_E, R);
  // mean vector
  vector[n_theta + nI] mu = rep_vector(0, n_theta+nI); // Put in zero-vector + støy? sigma_E
  // Basis priors
  l ~ normal(1,1);
  sigma ~ normal(1,0.5);
  sigma_theta ~ normal(prior_sigma_theta_mu, prior_sigma_theta_var);
  
  // Basis priors
  R ~ normal(0.25, 0.2);
  
  // Informative priors
  // R ~ normal(0.2, 0.001);
  
  // Uninformative priors
  // R ~ normal(0, 0.5);
  
  
  

  y~ multi_normal_cholesky(mu, L_K);
}




// generated quant_Et_Ees {
//   vector[n_theta_pred] f_P;
//   vector[n_theta_pred] y_P;
//   vector[nI_pred] f_I;
//   vector[nI_pred] y_I;
//   
//   f_P = Pgp_pred_rng(t_theta, t_E, yP, yI, t_theta_pred, sigma, l, R, C, sigma_theta, sigma_E);
//   for (n1 in 1:n_theta_pred)
//     y_P[n1] = normal_rng(f_P[n1], sigma_theta);
//   
//   f_I = Igp_pred_rng(t_theta, t_E, yP, yI, t_E_pred, sigma, l, R, C, sigma_theta, sigma_E);
//   for (n2 in 1:nI_pred)
//     y_I[n2] = normal_rng(f_I[n2], sigma_E);
// }




