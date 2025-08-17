////// fit_cholera_dpm_vax.stan
//
// This script defines the information necessary for conducting Bayesian MCMC methods of 
// parameter estimation using rstan in the format that is recognizable by the rstan function:
// stan(). This information is necessary to estimate parameters values specific to cholera
// transmission in Ethiopia during the 2023 and 2024 outbreaks.
// Key inputs include: the observed case data, the model (cholera_dpm_vax.R function), priors and a likelihood funtion. 
//
// This script doesn't directly generate outputs but is used by the stan() function
// in the 'fit_cholera_dpm_vax.R' script. 
//
//////


// define the transmission model fucntion
functions {
  real[] cholera_dpm_vax(real t, real[] y, real[] theta,
         real[] x_r, int[] x_i) {
       
      // declare initial values
      real S = y[1];
      real V1 = y[2];
      real V2 = y[3];
      real E = y[4];
      real Ev = y[5];
      real I1 = y[6]; 
      real I2 = y[7];
      real I3 = y[8];
      real I4 = y[9];
      real Ia = y[10];
      real Rc = y[11];
      real Rcp = y[12];
      real Ra = y[13];
      real BH = y[14];
      real BL = y[15];
      real CRInc = y[16];
      //
      
       
      // parameters to be fit to 
      real beta_e0 = theta[1]; 
      real beta_h = theta[2];
      real alpha_b = theta[3];
      real alpha_n = theta[4];
      
      
      // define fixed parameters
      real mu = x_r[1];        
      real gamma = x_r[2];	      
      real j = x_r[3];
      real f = x_r[4];
      real g = x_r[5];
      real omega_cf = x_r[6];
      real omega_cp = x_r[7];
      real omega_a = x_r[8]; 
      real theta_12 = x_r[9];
      real theta_23 = x_r[10];
      real theta_34 = x_r[11];
      real tau_c = x_r[12];
      real tau_a = x_r[13];
      real tau_t = x_r[14];
      real delta_1 = x_r[15];
      real delta_2 = x_r[16];
      real zeta = x_r[17];
      real sigma_2 = x_r[18];
      real sigma_3 = x_r[19];
      real sigma_4 = x_r[20];
      real rr_cd = x_r[21];
      real rr_t = x_r[22]; 
      real v_H = x_r[23];
      real b = x_r[24]; 
      real n_0 = x_r[25];
      real v_L = x_r[26];
      real KL = x_r[27];
      real K = x_r[28];
      real temp_mean = x_r[29];
      real temp_max = x_r[30];
      real rain_max = x_r[31];
      real A_temp = x_r[32];
      real A_rain = x_r[33];
      real phi_temp = x_r[34];
      real phi_rain = x_r[35];
      real rain_mean = x_r[36];
      real lag_tr = x_r[37]; 
      real omega_v1 = x_r[38];
      real omega_v2 = x_r[39];
      real a = x_r[40]; 
      real l1 = x_r[41];
      real l2 = x_r[42];
      real q = x_r[43];
      real r = x_r[44];
      real hyg_0 = x_r[45]; 
      real san_0 = x_r[46]; 
      real cwa_0 = x_r[47];
      real rho1 = x_r[48]; // will be zero unless a timestep changes it 
      real rho2 = x_r[49]; // will be zero unless a timestep changes it 
      
      
      // time-dependent parameters; as seen in the main function
      real P = S+V1+V2+E+Ev+I1+I2+I3+I4+Ia+Rc+Rcp+Ra;
      
      
        // interventional forcing terms 
      if (t >= 1930 && t < 2220) { //month 64 to 74
        hyg_0 = 0.5; // hyg_0 is replacing hyg_efct, defined in the model code, to prevent coding in a function; no linear increases during fitting time frame
      }
      if (t >= 1930 && t < 2220) { //month 64 to 74
        san_0 = 0.5; // san_0 is replacing san_efct, defined in the model code, to prevent coding in a function; no linear increases during fitting time frame
      }
      if (t >= 1930 && t < 2220) { //month 64 to 74
        cwa_0 = 0.8;  // cwa_0 is replacing cwa_efct, defined in the model code, to prevent coding in a function; no linear increases during fitting time frame
      }
      
      real SAN = 1 - san_0;
      real HYG = 1 - hyg_0;
      real CWA = 1 - cwa_0;
      
        // seasonality 
      real temp = A_temp * sin((pi()/6) * ((t/30+lag_tr) - phi_temp)) + temp_mean;
      real rain = A_rain * sin((pi()/6) * ((t/30+lag_tr) - phi_rain)) + rain_mean;
      
      
      real beta_e = beta_e0 * (1 + alpha_b * ((rain - rain_mean)/(rain_max - rain_mean)));
      
        // transmission forcing terms
      if (t >= 1790 && t < 1930) { 
        beta_e = beta_e + beta_e0 * 1.7;
        beta_h = beta_h * 5;
      }
      
      real KH = KL/700;
    
      real lam_e = beta_e * (BL/(KL + BL)) * CWA;
      real lam_h = beta_h * (BH/(KH + BH)) * HYG;
    
      real u = zeta * (I1 + sigma_2*I2 + sigma_3*I3 + sigma_4*I4 + sigma_4*Ia) * SAN;
      real n = n_0 * (1 + alpha_n * ((temp - temp_mean)/(temp_max - temp_mean)));
      
      
    // manually chnage the values of rho1 and rho2 based on time-dependent outputs from the base model; prevent coding in a function that slows run time and complicates math
      // rho1 triggering
      if (t >= 164 && t <= 167) {
        rho1 = 0.001368152;
      } else if (t >= 181 && t <= 185) {
        rho1 = 0.005899097;
      } else if (t >= 711 && t <= 719) {
        rho1 = 0.01141273;
      } else if (t >= 860 && t <= 866) {
        rho1 = 0.011377176;
      } else if (t >= 880 && t <= 887) {
        rho1 = 0.006952371;
      } else if (t >= 1038 && t <= 1042) {
        rho1 = 0.001434120;
      } else if (t >= 1073 && t <= 1082) {
        rho1 = 0.006707343;
      } else if (t >= 1193 && t <= 1200) {
        rho1 = 0.001300005;
      } else if (t >= 1223 && t <= 1231) {
        rho1 = 0.001886968;
      } else if (t >= 1453 && t <= 1460) {
        rho1 = 0.0008151785;
      } else if (t >= 1575 && t <= 1584) {
        rho1 = 0.01268335;
      } else if (t >= 1660 && t <= 1670) {
        rho1 = 0.01346582;
      } else if (t >= 1696 && t <= 1703) {
        rho1 = 0.01584711;
      } else if (t >= 1751 && t <= 1761) {
        rho1 = 0.008993750;
      } else if (t >= 1769 && t <= 1775) {
        rho1 = 0.008324294;
      } else if (t >= 1900 && t <= 1906) {
        rho1 = 0.008324294;
      }
     
      // rho2 triggering
       if (t == 181) {
        rho2 = 0.8181415;
      } else if (t == 182) {
        rho2 = 0.3133844;
      } else if (t == 183) {
        rho2 = 0.1835815;
      } else if (t == 184) {
        rho2 = 0.1285779;
      } else if (t == 353) {
        rho2 = 0.4651019;
      } else if (t >= 354 && t <= 363) {
        rho2 = 0.4605170;
      } else if (t >= 782 && t <= 788) {
        rho2 = 0.7675284;
      } else if (t == 950) {
        rho2 = 0.2037138;
      } else if (t >= 951 && t <= 957) { 
        rho2 = 0.6578815;
      } else if (t >= 1056 && t <= 1060) {
        rho2 = 1.15192925;
      } else if (t >= 1173 && t <= 1179) {
        rho2 = 0.7675284;
      }
      

      // ODEs
      real dS_dt = mu*P + omega_a*Ra + omega_cp*Rcp + omega_v1*V1 + omega_v2*V2 - rho1*S - (lam_e+lam_h)*S - mu*S;
    
      real dV1_dt = rho1*S - omega_v1*V1 - rho2*V1 - (1-l1)^a*(lam_e+lam_h)*V1 - mu*V1;
      real dV2_dt = rho2*V1 - omega_v2*V2 - (1-l2)^a*(lam_e+lam_h)*V2 - mu*V2;
    
      real dE_dt = (lam_e+lam_h)*S - (j+f+g)*gamma*E - (1-j-f-g)*gamma*E - mu*E;
      real dEv_dt = (1-l1)^a*(lam_e+lam_h)*V1 + (1-l2)^a*(lam_e+lam_h)*V2 - ((q*j)+(q*f)+(r*g))*gamma*Ev - (1-(q*j)-(q*f)-(r*g))*gamma*Ev - mu*Ev;
    
      real dI1_dt = j*gamma*(E + q*Ev) - (1-tau_t)*delta_1*theta_12*I1 - tau_t*theta_12*I1 - (1-tau_t)*(1-delta_1)*theta_12*I1 - mu*I1;
      real dI2_dt = f*gamma*(E + q*Ev) + (1-tau_t)*(1-delta_1)*theta_12*I1 - (1-tau_t)*delta_2*theta_23*I2 - tau_t*theta_23*I2 - (1-tau_t)*(1-delta_2)*theta_23*I2 - mu*I2;
      real dI3_dt = g*gamma*(E + r*Ev) + (1-tau_t)*(1-delta_2)*theta_23*I2 - theta_34*I3 - mu*I3;
      real dI4_dt = theta_34*I3 - tau_c*I4 - mu*I4;
    
      real dIa_dt = (1-j-f-g)*gamma*E + (1-(q*j)-(q*f)-(r*g))*gamma*Ev - tau_a*Ia - mu*Ia;
    
      real dRc_dt = tau_t*theta_12*I1 + tau_t*theta_23*I2 + tau_c*I4 - omega_cf*Rc - mu*Rc;
      real dRcp_dt = omega_cf*Rc - omega_cp*Rcp - mu*Rcp;
    
      real dRa_dt = tau_a*Ia - omega_a*Ra - mu*Ra;
    
      real dBH_dt = u - b*v_H*BH - (1-b)*v_H*BH;
      real dBL_dt = n*BL*(1-(BL/K)) + b*v_H*BH - v_L*BL;
    
      real dCDc_dt = (1-tau_t)*delta_1*theta_12*I1 + (1-tau_t)*delta_2*theta_23*I2;
      real dCInc_dt = (lam_e+lam_h)*S + (1-l1)^a*(lam_e+lam_h)*V1 + (1-l2)^a*(lam_e+lam_h)*V2;
      
      real dCRD_dt = ((rr_cd * (1-tau_t)*delta_1*theta_12*I1) + (rr_cd * (1-tau_t)*delta_2*theta_23*I2)); 
      
      real dCRInc_dt = ((rr_cd * (1-tau_t)*delta_1*theta_12*I1) + 
                 (rr_t * tau_t*theta_12*I1) +                  
                 (rr_cd * (1-tau_t)*delta_2*theta_23*I2) +      
                 (rr_t * tau_t*theta_23*I2)); 
                 
              
      return {dS_dt, 
              dV1_dt, dV2_dt, 
              dE_dt, dEv_dt, 
              dI1_dt, dI2_dt, dI3_dt, dI4_dt, dIa_dt,
              dRc_dt, dRcp_dt, dRa_dt,
              dBH_dt, dBL_dt, dCRInc_dt};
  }
  // solve the model at each timestep
  matrix solve_ode(real[] ts, real[] y0, real[] theta, data real[] x_r, data int[] x_i) { 
    return(to_matrix(integrate_ode_rk45(cholera_dpm_vax, y0, 0.0, ts, theta, x_r, x_i)));
  }
}


// define all data to be used by stan
data {
  // time data
  int<lower=1> n_days_sim;
  int<lower=1> n_months_sim;
  int<lower=0> burnin_months;
  int<lower=0> n_months_data;
  
  // fixed parameters (EVERY PARM in x_r)
  real<lower=0> mu;        
  real<lower=0> gamma;      
  real<lower=0> j;
  real<lower=0> f;
  real<lower=0> g;
  real<lower=0> omega_cf;
  real<lower=0> omega_cp;
  real<lower=0> omega_a;
  real<lower=0> theta_12;
  real<lower=0> theta_23;
  real<lower=0> theta_34;
  real<lower=0> tau_c;
  real<lower=0> tau_a;
  real<lower=0> tau_t;
  real<lower=0> delta_1;
  real<lower=0> delta_2;
  real<lower=0> zeta;
  real<lower=0> sigma_2;
  real<lower=0> sigma_3;
  real<lower=0> sigma_4;
  real<lower=0> rr_cd;
  real<lower=0> rr_t;
  real<lower=0> v_H;
  real<lower=0> b;
  real<lower=0> n_0;
  real<lower=0> v_L;
  real<lower=0> KL;
  real<lower=0> K;
  real<lower=0> temp_mean;
  real<lower=0> temp_max;
  real<lower=0> rain_max;
  real<lower=0> A_temp;
  real<lower=0> A_rain;
  real<lower=0> phi_temp;
  real<lower=0> phi_rain;
  real<lower=0> rain_mean;
  real<lower=0> lag_tr;
  real<lower=0> omega_v1;
  real<lower=0> omega_v2;
  real<lower=0> a;
  real<lower=0> l1;
  real<lower=0> l2;
  real<lower=0> q;
  real<lower=0> r;
  real<lower=0> hyg_0;
  real<lower=0> san_0;
  real<lower=0> cwa_0;
  real<lower=0> rho1;
  real<lower=0> rho2;
  
  array[16] real y0; // ode intial conditions

  array[n_months_data] int rinc_dat;   // number of obserevd data points
  
  int compute_likelihood; // for prior predictive checks 
}


// change format of data so recognizable by the model function 
transformed data {
  array[49] real x_r; // real numbers
  array[0] int x_i; // integars
  
  // RELIST EVERY FIXED PARAMETER
  x_r[1] = mu;        
  x_r[2] = gamma;      
  x_r[3] = j;
  x_r[4] = f;
  x_r[5] = g;
  x_r[6] = omega_cf;
  x_r[7] = omega_cp;
  x_r[8] = omega_a;
  x_r[9] = theta_12;
  x_r[10] = theta_23;
  x_r[11] = theta_34;
  x_r[12] = tau_c;
  x_r[13] = tau_a;
  x_r[14] = tau_t;
  x_r[15] = delta_1;
  x_r[16] = delta_2;
  x_r[17] = zeta;
  x_r[18] = sigma_2;
  x_r[19] = sigma_3;
  x_r[20] = sigma_4;
  x_r[21] = rr_cd;
  x_r[22] = rr_t;
  x_r[23] = v_H;
  x_r[24] = b;
  x_r[25] = n_0;
  x_r[26] = v_L;
  x_r[27] = KL;
  x_r[28] = K;
  x_r[29] = temp_mean;
  x_r[30] = temp_max;
  x_r[31] = rain_max;
  x_r[32] = A_temp;
  x_r[33] = A_rain;
  x_r[34] = phi_temp;
  x_r[35] = phi_rain;
  x_r[36] = rain_mean;
  x_r[37] = lag_tr;
  x_r[38] = omega_v1;
  x_r[39] = omega_v2;
  x_r[40] = a;
  x_r[41] = l1;
  x_r[42] = l2;
  x_r[43] = q;
  x_r[44] = r;
  x_r[45] = hyg_0;
  x_r[46] = san_0;
  x_r[47] = cwa_0;
  x_r[48] = rho1;
  x_r[49] = rho2;
  
  array[n_days_sim] real ts; // define timestep
  for (ii in 1:n_days_sim) {
    ts[ii] = ii;  
  }
}


parameters { // parameters I want to provide priors for and fit to (if estimating transformed parms, this should be the value used to calc)
  // define 'modeled parameters' for which priors will be provided :: parameter space constructed from these
  real<lower=0, upper=1> beta_e0;               
  real<lower=0, upper=1> beta_h; 
  real<lower=0> alpha_b;
  real<lower=0, upper=1.9> alpha_n;
}


transformed parameters {
  // define 'transformed paramteres' or derived quantities using modeled parameters
  // make parameter vector for simulation
  array[4] real theta; // increase by one for new parameter
  theta[1] = beta_e0;       
  theta[2] = beta_h; 
  theta[3] = alpha_b;
  theta[4] = alpha_n;
  
  
  // simulate transmission model
  //array[n_days_sim, 19] real<lower=0> y; // make an array to hold the output of each timestep, 5 for each equation 
  array[n_days_sim, 16] real<lower=0> y;
  
  //int n_months_sim = floor(n_days_sim/30);
  //int n_months_sim = n_days_sim/30;
  array[n_months_sim] real<lower=0> mreported_inc;
  
  y = integrate_ode_rk45(cholera_dpm_vax, y0, 0.0, ts, theta, x_r, x_i); // WHAT DOES 0.0 REFER TO
  
  //16s used to be 19s before removing equations
  mreported_inc[1] = y[30, 16]; // from model output, calculate mean reported weekly incidence
  for (ii in 2:n_months_sim) {
    mreported_inc[ii] = (y[30*ii, 16] - y[30*(ii - 1), 16]) ;
  }
}


model { 
  // modeled parameter priors
  beta_e0 ~  beta(2,8); //is beta(2,8) generic enough for it not be infomred by literature
  beta_h ~  beta(2,8);
  alpha_b ~ normal(0,1);
  alpha_n ~ normal(0,1);
 
  // likelihood
  if (compute_likelihood == 1) {
    for (ii in 1:n_months_data){ // match the data with the right simulated months here 
      rinc_dat[ii] ~ poisson(mreported_inc[ii + burnin_months]); 
    }
  }
}


generated quantities {
  // derive relevant quantities after each simulation
  array[n_months_sim] real<lower=0> pred_reported_inc;
 
  for (ii in 1:n_months_sim) {
    pred_reported_inc[ii] = poisson_rng(mreported_inc[ii]); // given mean cases generate pred_cases using poisson distribution
  }

  // log-likelihood for model comparison :: probabiltiy of observing the data given the parameters
  vector[n_months_sim] log_lik_cases;
  for (ii in 1:n_months_data) {
    log_lik_cases[ii] = poisson_lpmf(rinc_dat[ii] | mreported_inc[ii + burnin_months]); // log-likelihood for case data
  }
}

