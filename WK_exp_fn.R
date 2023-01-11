
# 
# # create data in the [0,1] scale
# create_data = function(flow, t_Eme, Rtrue, Ctrue, Ztrue=NULL, nc=3, n_theta=20, nI=17, n_pred=50, theta_noise=3, E_noise=4, seed = 0){
#   # 1. simulate data from WK2 or WK3
#   if(is.null(Ztrue)){
#     Psim = WK2_simulate(flow = flow, t_Eme = t_Eme, R = Rtrue, C = Ctrue)
#   }else{
#     Psim = WK3_simulate(flow = flow, t_Eme = t_Eme, R = Rtrue, C = Ctrue, Z=Ztrue)
#   }
#   
#   # 2. selected indices according to n_theta and nI
#   nflow = length(flow)
#   indP = round(seq(1, nflow, length.out = n_theta)); indI = round(seq(1, nflow, length.out = nI))
#   t_theta = matrix(d$t_Eme[indP], ncol = 1); t_E = matrix(d$t_Eme[indI], ncol = 1)
#   yP_real = Psim[indP]; yI_real = flow[indI]
#   # 3. Add noise
#   set.seed(seed)
#   theta_noise = rnorm(n_theta*nc, 0, theta_noise)
#   E_noise =rnorm(nI*nc, 0, E_noise)
#   yP_real = rep(yP_real,each=nc)
#   yI_real = rep(yI_real,each=nc)
#   
#   yP_temp = yP_real + theta_noise
#   yI_temp = yI_real + E_noise
#   y_temp = c(yP_temp, yI_temp)
#   # tranform y in [0,1]
#   rPI = range(c(yP_temp, yI_temp))
#   yP = (yP_temp - rPI[1])/diff(rPI)
#   yI = (yI_temp - rPI[1])/diff(rPI)
#   
#   ind_pred = round(seq(1,101, length.out = n_pred))
#   data_noisy_pred = list(n_theta = nc*n_theta, nI = nc*nI, t_theta = matrix(rep(t_theta,each=nc),ncol=1)
#                          , t_E = matrix(rep(t_E,each=nc),ncol=1), yP = yP, yI = yI
#                          , t_theta_pred = matrix(d$t_Eme[ind_pred],ncol = 1), n_theta_pred = n_pred
#                          , t_E_pred = matrix(d$t_Eme[ind_pred],ncol = 1), nI_pred = n_pred
#                          )
#   
#   # true data for plott_Eng
#   data_mod_true = data.frame(t_Eme = t_Eme, I = flow, P = Psim)
#   return(list(data_noisy_pred = data_noisy_pred, data_mod_true =data_mod_true, y=y_temp))
# }
# 
# pendulum <- function(theta0,L,runt_Eme,step = 0.001){
#   g <- 9.81
#   period <- 2*pi*sqrt(L/g)
#   
#   t_Emes <- seq(0.0,runt_Eme, by = step)
#   
#   angle=theta(t_Emes, period = period,theta0 = theta0)
#   
#   pen.x <- L*sin(theta(t_Emes, period = period,theta0 = theta0)) 
#   pen.y <- -L*cos(theta(t_Emes,period = period,theta0 = theta0))
#   
#   out_thetaut <- matrix(c(t_Emes, angle),ncol = 2)
#   return(out_thetaut)
# }
# 
# Psim =pendulum(theta0, L, maxt_Eme)
# plot(Psim[,1], Psim[,2])

# create data in the [0,1] scale
create_data = function(angle, time, Rtrue, Ctrue, theta0, L, Ztrue=NULL, nc=1, n_theta=20, nI=20, n_pred=50, theta_noise, E_noise, prior_sigma_theta_mu =0.01, prior_sigma_theta_var =0.05, range_sigma_theta_upper =0.02, seed = 0){

  # Define pendulum and the times the pendulum is observed
  t_theta=time
  t_E=time+0.05
  # # 2. selected indices according to n_theta and nI
  # nangle = length(angle)
  # indP = round(seq(1, nangle, length.out = n_theta)); indI = round(seq(0, nangle, length.out = nI))
  # print(indP)
  # print(indI)
  # t_theta = matrix(t_Eme[indP], ncol = 1); t_E = matrix(t_Eme[indP], ncol = 1)
  # yP_real = angle[indP]; yI_real = rep(0, length = length(yP_real))
  
  yP_real = angle; yI_real = rep(0, length = length(yP_real))
  
  # 3. Add noise
  set.seed(seed)
  theta_noise = rnorm(n_theta*nc, 0, theta_noise)
  E_noise =rnorm(nI*nc, 0, E_noise)
  yP_real = rep(yP_real,each=nc)
  yI_real = rep(yI_real,each=nc)

  
  yP_temp = yP_real + theta_noise
  yI_temp = yI_real + E_noise
  y_temp = c(yP_temp, yI_temp)
  
  # tranform y in [0,1] KAN KANSKJE DROPPE
  # rPI = range(c(yP_temp, yI_temp))
  # yP = (yP_temp - rPI[1])/diff(rPI)
  # yI = (yI_temp - rPI[1])/diff(rPI)
  
  #HVIS DROPPE OVER
  yP = yP_temp
  yI = yI_temp
  
  # ind_pred = round(seq(1,101, length.out = n_pred))
  # data_noisy_pred = list(n_theta = nc*n_theta, nI = nc*nI, t_theta = matrix(rep(t_theta,each=nc),ncol=1)
  #                        , t_E = matrix(rep(t_E,each=nc),ncol=1), yP = yP, yI = yI
  #                        , t_theta_pred = matrix(t_Eme[ind_pred]+0.005,ncol = 1), n_theta_pred = n_pred
  #                        , t_E_pred = matrix(t_Eme[ind_pred]+0.005,ncol = 1), nI_pred = n_pred
  # )
  
  data_noisy_pred = list(n_theta = nc*n_theta, nI = nc*nI, t_theta = matrix(rep(t_theta,each=nc),ncol=1)
                         , t_E = matrix(rep(t_E,each=nc),ncol=1), yP = yP, yI = yI
                         , prior_sigma_theta_mu = prior_sigma_theta_mu, prior_sigma_theta_var = prior_sigma_theta_var
                         , range_sigma_theta_upper =range_sigma_theta_upper
  )
  
  # true data for plott_Eng
  data_mod_true = data.frame(time =time, I = 0, P = angle)
  return(list(data_noisy_pred = data_noisy_pred, data_mod_true =data_mod_true, y=y_temp))
}

# transform to the real scale
transform_post = function(y, fit){
  post_df = as.data.frame(extract(fit))
  ry = diff(range(y))
  my = min(y)
  cn = colnames(post_df)
  
  # grep("_P", cn)
  # 
  # #post_df[,c("R")] = post_df[,c("R")]*2.5+0.5 #Why *2.5+0.5 here
  # indPI = c(grep("_P", cn), grep("_I", cn)) #This returns integer(0), dont understand
  # post_df[,indPI] = post_df[,indPI]*ry + my #data frame with 0 columns and 1500 rows = data frame with 0 columns and 0 rows
  # # ind_delta = grep("delta",cn)
  # ind_delta = c(grep("sigma",cn), grep("mu_wk2",cn))
  # md = apply(post_df[,ind_delta],2,min)
  # post_df[,ind_delta] = post_df[,ind_delta]*ry + md
  # post_df[,y]
  return(post_df)
}



