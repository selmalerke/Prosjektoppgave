#Load packages
library(foreign)
library(xtable)
library(stargazer)
library(rstan)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gmodels)
library(pracma)
library(BayesFactor)
theme_set(theme_classic()) # set ggplot theme
rstan_options(auto_write = TRUE)
options(mc.cores = 3) # allocate 3 cores (for each model we run 3 chains in parallel)

#Load fucntion from other file. 
source("WK_exp_fn.R") 
stanfile = 'Pendulum_model.stan'
# source("/Users/selmalerke/BC-with-PI-priors/functions/WK_exp_fn.R")
# stanfile = '/Users/selmalerke/Pendulum_model.stan'

#pi10
angledamped <- scan("angledamped.txt", what="", sep="\n")
timesdamped <- scan("timesdamped.txt", what="", sep="\n")

#pi5
# angledamped <- scan("angledampedpi5.txt", what="", sep="\n")
# timesdamped <- scan("timesdampedpi5.txt", what="", sep="\n")

#Computer
#angledamped <- scan("/Users/selmalerke/A Prosjektoppgave/file.txt", what="", sep="\n")
#timesdamped <- scan("/Users/selmalerke/A Prosjektoppgave/file2.txt", what="", sep="\n")
#angledamped <- scan("/Users/selmalerke/PycharmProjects/fysmat2/angledampedpi5.txt", what="", sep="\n")
#timesdamped <- scan("/Users/selmalerke/PycharmProjects/fysmat2/timesdampedpi5.txt", what="", sep="\n")


#Values used
L=2 #Length of pendulum rod. BASIS = 2
g=10 #Gravitational constant. BASIS = 10
R_real = L/g #The real value, which we want to estimate
theta0=pi/5 #Starting angle in radians BASIS = pi/5
sigma_theta=0.03 #Angle measurement noise. BASIS = 0.03
prior_sigma_theta_mu=0.02 #Goes into STAN code. Mean of sigma_theta. BASIS = 0.02
prior_sigma_theta_var = 0.04 #Goes into STAN code. Variance of sigma_theta. BASIS = 0.04
range_sigma_theta_upper = 0.06 #Goes into STAN code. Upper value of sigma_theta. BASIS = 0.06
realval <- data.frame(variable = c("sigma_theta", "R"), real = c(sigma_theta, R_real))

runtime = 3.5 #How long the pendulum will run for
nummeasurements = 12
steplength = runtime/nummeasurements #Length betweeen each measurements

#Values only used in true model
m = 3.0
I = 12.0
h = 0.001

#When running for credible intervals
sigma_theta_list = seq(0.01, 0.05, length.out = 100)

#Titles
plottitle = " "
plotsubtitle = " "
histtitle = "Posterior means"
savehisttitle = "Postmeans6.pdf"
devplottitle = "Deviation from true R-value"
savedevplot = "Dev.pdf"
ciplottitle = "Width of credible interval"
saveciplot = "Ci.pdf"
savehisttitle_free = "Histogram_freescales.pdf"
comparehisttitleR = " "
comparehisttitlesigmatheta = " "
savecomparehisttitleR = "CompareR1.pdf"
savecomparehisttitlesigmatheta = "Comparesigmatheta1.pdf"



#seeds=as.integer(seq(10, 1800, 350)) #For test runs
seeds=as.integer(seq(10, 1800, length.out = 100)) #For run on NTNU computer

#PENDULUM SIMULATION

# 1. Linearized model

pendulum <- function(R, theta0, runtime, nummeasurements){
  times = seq(0.0, runtime,length.out = nummeasurements) 
  angle = theta0*cos(sqrt(1/R)*times)
  df = data.frame(Time = times, Angle = angle)
}


#Function to plot pendulum
makeplot <- function(data, xaxis, yaxis){
  ggplot(data = data, aes(xaxis, yaxis)) + 
    theme_bw() + 
    geom_point(col = "tomato1", size = 2) + 
    labs(title = plottitle, subtitle = plotsubtitle) + 
    labs(x = "Time(s)", y = expression(paste("Angle ", theta, " (rad)"))) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5))
}

# makeplot_with_and_without_noise = function(data, xaxis, yaxis, wnoise, xwaxis, ywaxis){
#   ggplot(data = data, aes(xaxis, yaxis)) + 
#     theme_bw() + geom_point(col = "blue", size = 0.5) + geom_point(data = wnoise, aes(x = xwaxis, y = ywaxis), col = "tomato1", shape = 3, size = 3, stroke = 1.5)+
#     labs(title = plottitle, subtitle = plotsubtitle) + 
#     labs(x = "Time(s)", y = expression(paste("Angle ", theta, " (rad)"))) +
#     theme(plot.title = element_text(hjust = 0.5, size = 18),
#           plot.subtitle = element_text(hjust = 0.5))+
#     scale_shape_manual(name = "Data", values = c(Noisy = 3, Original = 16)) +  scale_color_manual(name = "Data", values = c(Noisy = "tomato1", Original = "blue"))
#   
# }
makeplot_with_and_without_noise = function(data, xaxis, yaxis, wnoise, xwaxis, ywaxis){
  ggplot(data = data, aes(xaxis, yaxis)) + 
    theme_bw() +  
    geom_point(data = wnoise, aes(x = xwaxis, y = ywaxis, shape = "Noisy", color = "Noisy"), size = 3, stroke = 1.5) +
    geom_point(aes(shape = "True", color = "True"), col = "blue", size = 0.3) + 
    labs(title = plottitle, subtitle = plotsubtitle) + 
    labs(x = "Time(s)", y = expression(paste("Angle ", theta, " (rad)"))) +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5))+
    scale_shape_manual(name = "Data", values = c(Noisy = 3, True = 16)) +
    scale_color_manual(name = "Data", values = c(Noisy = "tomato1", True = "blue"))
}

# df_true=pendulum_true(theta0, runtime, h, m, g, L, I, nummeasurements)
# df_t_cont =pendulum_true(theta0, runtime, h, m, g, L, I, 10000)
# makeplot(df_true, df_true$Time, df_true$Angle)
# 
# 
#Plot with noise
# plot_it = function(df_big, df_small){
#   wnoise = create_data(df_small$Angle, df_small$Time, Rtrue = L/g, Ctrue =1, theta0, L, Ztrue=0, nc=1, n_theta=as.integer(length(df$Time)), nI=as.integer(length(df$Time)), n_pred=50, sigma_theta, E_noise=1e-3, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, seed = 18)
#   df_t =data.frame(Time = wnoise$data_noisy_pred$t_theta, Angle = wnoise$data_noisy_pred$yP)
#   makeplot_with_and_without_noise(df_big, df_big$Time, df_big$Angle, df_t, df_t$Time, df_t$Angle)
# }
# 
# df  =pendulum(R_real, theta0, runtime, nummeasurements)
# df_cont =pendulum(R_real, theta0, runtime, 10000)
# 
# df_damped = pendulum_damped(angledamped, timesdamped, nummeasurements)
# 
# set.seed(16)
# plot_it(df_cont, df_damped)


compare_with_noise <- function(df, xaxis, yaxis, sigma_theta){
  df2 = create_data(df$Angle, df$Time, Rtrue = L/g, Ctrue =1, theta0, L, Ztrue=0, nc=1, n_theta=as.integer(length(df$Time)), nI=as.integer(length(df$Time)), n_pred=50, sigma_theta, E_noise=1e-3, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, seed = 12)
  df3 = df2$data_noisy_pred
  xaxis = df3$t_theta
  yaxis = df3$yP
  
  dt <- 0.01  # Time step
  t <- seq(0, runtime, by = dt)  # Time values
  
  # Calculate the pendulum position at each time step
  position <- cos(sqrt(R_real) * t)
  withoutnoise = data.frame(t, position)
  # Plot the pendulum position over time
  plot(t, position, type = "l", xlab = "Time (s)", ylab = "Position (m)")
  
  ggplot(data = df2, aes(xaxis, yaxis)) + 
    theme_bw() + 
    geom_point(col = "tomato1", size = 2) + 
    labs(title = plottitle, subtitle = plotsubtitle) + 
    labs(x = "Time(s)", y = expression(paste("Angle ", theta, " (rad)"))) + geom_line(data = withoutnoise, aes(t, position))+
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5))
}


# Define pendulum 
df=pendulum(R_real, theta0, runtime, nummeasurements)

#Make plot of pendulum motion
makeplot(df, df$Time, df$Angle)


#2. True model

pendulum_true <- function(theta0, runtime, h, m, g, L, I, nummeasurements) {
  t = c( seq(0, runtime, by=h))
  n = length(t)
  
  y=c(rep(0, n))
  v=c(rep(0, n))
  
  accel<-function(theta){
    return(-(m*g*L/I)*sin(theta))
  }
  
  y[1] = theta0
  v[1] = 0
  
  for (i in 1:n){
    k1y = h*v[i]
    k1v = h*accel(y[i])
    
    k2y = h*(v[i]+0.5*k1v)
    k2v = h*accel(y[i]+0.5*k1y)
    
    k3y = h*(v[i]+0.5*k2v)
    k3v = h*accel(y[i]+0.5*k2y)
    
    k4y = h*(v[i]+k3v)
    k4v = h*accel(y[i]+k3y)
    
    # Update next value of y
    y[i+1] = y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0
    v[i+1] = v[i] + (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0
  }
  y=y[1:n]
  
  #extracting every 70th element for shorter run time of MCMC
  t <- t[seq(1, length(t), length.out = nummeasurements)]
  y <- y[seq(1, length(y), length.out = nummeasurements)]
  
  df_true <- data.frame(Time=t, Angle=y)
  return(df_true)
}

df_true=pendulum_true(theta0, runtime, h, m, g, L, I, nummeasurements)

# 3. Damped pendulum
#This is coded in python, jupiter
# Read in the data

pendulum_damped <- function(angledamped, times, nummeasurements){
  angledamped=as.double(angledamped)
  times = as.double(times)
  
  times <- times[seq(1, length(times), length.out = nummeasurements)]
  angledamped <- angledamped[seq(1, length(angledamped), length.out = nummeasurements)]
  df = data.frame(Time= times, Angle=angledamped)
  return(df)
}

df_damped = pendulum_damped(angledamped, timesdamped, nummeasurements)
makeplot(df_damped, df_damped$Time, df_damped$Angle)


#Fit model to the data
fit_model_once <- function(process, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper) {
  #meanlist = data.frame(l=rep(0, length(seeds)), sigma=rep(0, length(seeds)), sigma_theta=rep(0, length(seeds)), R=rep(0, length(seeds)))
  df2 = create_data(df$Angle, df$Time, Rtrue = L/g, Ctrue =1, theta0, L, Ztrue=0, nc=1, n_theta=as.integer(length(df$Time)), nI=as.integer(length(df$Time)), n_pred=50, sigma_theta, E_noise=1e-3, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, seed = 12)
  
  #Plot the angle with noise
  plot(df2$data_noisy_pred$t_theta, df2$data_noisy_pred$yP)
  
  #Fit the model
  result = stan(file=stanfile,
              data=df2$data_noisy_pred,
              chains=3,
              iter=1000, 
              seed=150
  )
 
  return(result)
}

# #PLOT PRIOR AND POSTERIOR
# df=pendulum(R_real, theta0, runtime, nummeasurements)
# df_true=pendulum_true(theta0, runtime, h, m, g, L, I, nummeasurements)
# 
# makeplot(df, df$Time, df$Angle)
# makeplot(df_true, df_true$Time, df_true$Angle)
# df2 = create_data(df_true$Angle, df_true$Time, Rtrue = L/g, Ctrue =1, theta0, L, Ztrue=0, nc=1, n_theta=as.integer(length(df$Time)), nI=as.integer(length(df$Time)), n_pred=50, sigma_theta, E_noise=1e-3, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, seed = 15)
# #Plot the angle with noise
# plot(df2$data_noisy_pred$t_theta, df2$data_noisy_pred$yP)
# 
# set.seed(16)
# testing_true = fit_model_once(df, theta0, L, sigma_theta, prior_sigma_theta_mu=sigma_theta, prior_sigma_theta_var, range_sigma_theta_upper)
#   # # Make histogram ready list
# y = df2$y # the original observed data
# pp_se= transform_post(y, testing_true)
# pp_se$R
# pl_df=pp_se[,names(testing_true)[3:4]]
# pl_df$sample=1:nrow(pl_df)
# m_pl_df = melt(pl_df, id="sample")
# m_pl_df
# 
# histogram = ggplot(data=m_pl_df)+scale_color_manual(values=c("thistle3", "thistle3", "thistle3", "thistle3"))+ggtitle(histtitle)+
#     geom_histogram(aes(x=value, y=..density.., color=variable), inherit.aes = FALSE, fill="lavender")+theme_classic()+facet_wrap(~variable,nrow = 2, scales = "free")+
#     theme(legend.position="none")+geom_vline(data=realval, aes(xintercept=real, color=variable), color = "springgreen4", size=1)+theme(text = element_text(size = 15))
# 
# 
# plot_prior_posterior <- function(posterior, mean, sd) {
#   ggplot() +
#     geom_density(aes(x = posterior), color = "blue") +
#     stat_function(fun = dnorm, args = list(mean = mean, sd = sd),
#                   color = "red", alpha = 0.5) +scale_color_manual(name = "Line Type", values = c("blue" = "Posterior", "red" = "Prior"))+
#     xlim(min(min(posterior), mean - 2 * sd), max(max(posterior), mean + 2 * sd))+xlab("Value") +
#     ylab("Density")+ geom_vline(xintercept = R_real, color = "springgreen4") + ggtitle("Prior and posterior")
# }
# histogram
# plot_prior_posterior(pp_se$R, 0.25, 0.2)
# plot_prior_posterior(pp_se$R, 0.2, 0.04)

  
fit_model <- function(df, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper) {
  
  meanlist = data.frame(l=rep(0, length(seeds)), sigma=rep(0, length(seeds)), sigma_theta=rep(0, length(seeds)), R=rep(0, length(seeds)))
  
  for (i in 1:length(seeds)){
    print("Round")
    print(i)
    #Sort data and add noise
    df2 = create_data(df$Angle, df$Time, Rtrue = L/g, Ctrue =1, theta0, L, Ztrue=0, nc=1, n_theta=as.integer(length(df$Time)), nI=as.integer(length(df$Time)), n_pred=50, sigma_theta, E_noise=1e-3, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper, seed = seeds[i])
    
    #Plot the angle with noise
    plot(df2$data_noisy_pred$t_theta, df2$data_noisy_pred$yP)
    
    #Plot the energy balance with noise
    plot(df2$data_noisy_pred$t_E, df2$data_noisy_pred$yI)
    
    #Fit the model
    temp = stan(file=stanfile,
                data=df2$data_noisy_pred,
                chains=3,
                iter=1000, 
                seed=seeds[i]
    )
    temp = as.data.frame(extract(temp))
    
    #Summary of the model
    meanlist$sigma_theta[i] = mean(temp$sigma_theta)
    meanlist$R[i] = mean(temp$R)
  }
  
  return(meanlist)
}

#Fit model to the data
fit_model_extract_ci_dev <- function(df, seeds, theta0, L, sigma_theta_list, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper) {
  
  prior_sigma_theta_mu=sigma_theta_list+0.01
  prior_sigma_theta_var = prior_sigma_theta_var
  range_sigma_theta_upper = sigma_theta_list+prior_sigma_theta_var
  meanlist_R = rep(0, length(sigma_theta_list))
  ci_R=rep(0, length(sigma_theta_list))
  
  meancred = data.frame(dev = meanlist_R, ci_R=ci_R)
  meanlist = data.frame(l=rep(0, length(seeds)), sigma=rep(0, length(seeds)), sigma_theta=rep(0, length(seeds)), R=rep(0, length(seeds)))
  print(meancred)
  for (i in 1:length(sigma_theta_list)){
    print("Round")
    print(i)
    #Sort data and add noise
    R_true = L/g
    df2 = create_data(df$Angle, df$Time, Rtrue = L/g, Ctrue =1, theta0, L, Ztrue=0, nc=1, n_theta=as.integer(length(df$Time)), nI=as.integer(length(df$Time)), n_pred=50, sigma_theta_list[i], E_noise=1e-3, prior_sigma_theta_mu[i], prior_sigma_theta_var, range_sigma_theta_upper[i], seed = seeds[i])
    
    #Fit the model
    temp = stan(file=stanfile,
                data=df2$data_noisy_pred,
                chains=3,
                iter=1000, 
                seed=seeds[i]
    )
    temp = as.data.frame(extract(temp))
    
    #Summary of the model
    meanlist$sigma_theta[i] = mean(temp$sigma_theta)
    meanlist$R[i] = mean(temp$R)
    
    # Calculate the credible interval
    ci <-quantile(unlist(temp$R), c(0.025, 0.975))
    
    # Extract the lower and upper bounds of the interval
    lower_bound <- ci[1]
    upper_bound <- ci[2]
    
    # Calculate the length of the interval
    length <- upper_bound - lower_bound
    
    print(unname(length))
    print(mean(temp$R))
    
    # Add interval lenght to the 
    meancred$ci_R[i] = unname(length)
    meancred$dev[i] = abs(mean(temp$R) - R_true)
  }
  
  return(meancred)
}


#Find means of the meanlists from the 100 fitted models
maketot_means = function(meanlist){
  #Summarized mean
  means = data.frame(matrix(nrow=2, ncol=0))
  means$variable = c("sigma_theta", "R")
  means$grp.mean = c( mean(meanlist$sigma_theta), mean(meanlist$R))
  return(means)
}

#With free scale
plot_histogram_free_scale <- function(data, tot_means, savehisttitle_free){
  
  histogram = ggplot(data=data)+scale_color_manual(values=c("thistle3", "thistle3", "thistle3", "thistle3"))+ggtitle(histtitle)+
    geom_histogram(aes(x=value, y=..density.., color=variable), inherit.aes = FALSE, fill="lavender")+theme_classic()+facet_wrap(~variable,nrow = 2, scales = "free")+
    geom_vline(data=tot_means, aes(xintercept=grp.mean, color=variable), linetype="dashed", size=1, col = "thistle4") + theme(legend.position="none")+geom_vline(data=realval, aes(xintercept=real, color=variable), color = "springgreen4", size=1)+theme(text = element_text(size = 15))

  ggsave(savehisttitle_free, plot = histogram)
  return(histogram)
}



plot_dev = function(sigma_theta_list, devlist, savedevplot){
  devlist= as.data.frame(devlist)
  
  plot =ggplot(data = devlist, aes(sigma_theta_list, devlist)) + theme_bw()+ geom_point(col="purple4", size=3)+labs(title="Deviation from mean of R", subtitle= "For different angle measurement noises")+labs(x= expression(paste(sigma, "(", theta, ")")), y=expression(paste("Deviation")))+theme(
    plot.title = element_text(hjust = 0.5, size=18),
    plot.subtitle = element_text(hjust = 0.5))+theme(text = element_text(size = 20))  + ylim(0,0.03)
  
  ggsave(savedevplot, plot = plot)
  return(plot)
}

plot_ci = function(sigma_theta_list, cilist, saveciplot){
  cilist= as.data.frame(cilist)
  
  plot =ggplot(data = cilist, aes(sigma_theta_list, cilist)) + theme_bw()+ geom_point(col="springgreen4", size=3)+labs(title="Width of 95% credible interval for R", subtitle= "For different angle measurement noises")+labs(x= expression(paste(sigma, "(", theta, ")")), y=expression(paste("Width")))+theme(
    plot.title = element_text(hjust = 0.5, size=18),
    plot.subtitle = element_text(hjust = 0.5))+theme(text = element_text(size = 20))  + ylim(0,0.05)
  
  ggsave(saveciplot, plot = plot)
  return(plot)
}

plot_multi_histogram <- function(df, feature, label_column, title, trueval) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="white") +
    geom_density(alpha=0.7) +
    labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=label_column))+geom_vline(xintercept = trueval, color = "springgreen4", size=1)+theme(text = element_text(size = 15))+ggtitle(title)+theme(text = element_text(size = 20))  
  
}

comparemodshist <- function(meanlist, meanlist2, firstmodname, secondmodname){
  a <-data.frame(R=meanlist$R, Model=rep(firstmodname, length(meanlist$R)))
  b <-data.frame(R=meanlist2$R, Model=rep(secondmodname, length(meanlist2$R)))
  many_distros <- do.call('rbind', list(a,b))
  
  Rcomparison = plot_multi_histogram(many_distros, 'R', 'Model', comparehisttitleR, R_real)
  ggsave(savecomparehisttitleR, plot = Rcomparison)
  
  a <-data.frame(sigma_theta=meanlist$sigma_theta, Model=rep(firstmodname, length(meanlist$sigma_theta)))
  b <-data.frame(sigma_theta=meanlist2$sigma_theta, Model=rep(secondmodname, length(meanlist2$sigma_theta)))
  many_distros <- do.call('rbind', list(a,b))
  
  sigmathetacomparison = plot_multi_histogram(many_distros, 'sigma_theta', 'Model', comparehisttitlesigmatheta, sigma_theta)
  ggsave(savecomparehisttitlesigmatheta, plot = sigmathetacomparison)
}

#TIRSDAG: ci og dev
# ci_dev = fit_model_extract_ci_dev(df, seeds, theta0, L, sigma_theta_list, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# ci_dev_true = fit_model_extract_ci_dev(df_true, seeds, theta0, L, sigma_theta_list, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# ci_dev_damped = fit_model_extract_ci_dev(df_damped, seeds, theta0, L, sigma_theta_list, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# 
# plotdev = plot_dev(sigma_theta_list, ci_dev$dev, "Dev_lin.pdf")
# plotdev
# plotci = plot_ci(sigma_theta_list, ci_dev$ci_R, "Ci_lin.pdf")
# plotci
# 
# plotdev = plot_dev(sigma_theta_list, ci_dev_true$dev, "Dev_true.pdf")
# plotdev
# plotci = plot_ci(sigma_theta_list, ci_dev_true$ci_R, "Ci_true.pdf")
# plotci
# 
# plotdev = plot_dev(sigma_theta_list, ci_dev_damped$dev, "Dev_damped.pdf")
# plotdev
# plotci = plot_ci(sigma_theta_list, ci_dev_damped$ci_R, "Ci_damped.pdf")
# plotci

# #TIRSDAG: Comparison plots:
# 
# #TRUE PROCESS
# meanlist_true = fit_model(df_true, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_true = meanlist_true[3:4]
# tot_means_true = maketot_means(meanlist_true)
# #
# # # Make histogram ready list
# meanlist_true$sample=1:nrow(meanlist_true)
# meanlist_freq_true <- melt(meanlist_true, id.vars = "sample")
# 
# #LINEAR PROCESS
# meanlist_lin = fit_model(df, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_lin = meanlist_lin[3:4]
# tot_means_lin = maketot_means(meanlist_lin)
# 
# # # Make histogram ready list
# meanlist_lin$sample=1:nrow(meanlist_lin)
# meanlist_freq_lin <- melt(meanlist_lin, id.vars = "sample")
# 
# #DAMPED PROCESS
# meanlist_damped = fit_model(df_damped, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_damped = meanlist_damped[3:4]
# tot_means_damped = maketot_means(meanlist_damped)
# 
# # # Make histogram ready list
# meanlist_damped$sample=1:nrow(meanlist_damped)
# meanlist_freq_damped <- melt(meanlist_damped, id.vars = "sample")
# 
# #SUMMARIES
# xtable(summary(meanlist_true))
# xtable(summary(meanlist_lin))
# xtable(summary(meanlist_damped))
# 
# #MAKE PLOTS
# plottrue = plot_histogram_free_scale(meanlist_freq_true, tot_means_true, "TP_histogram2_pi5_003.pdf")
# plotlin = plot_histogram_free_scale(meanlist_freq_lin, tot_means_lin, "LP_histogram2_pi5_003.pdf")
# plotdamped = plot_histogram_free_scale(meanlist_freq_damped, tot_means_damped, "DP_histogram2_pi5_003.pdf")
# 
# savecomparehisttitlesigmatheta = "Comp_TP_LP2_pi5_003.pdf"
# savecomparehisttitleR= "Comp_TP_LP_R2_pi5_003.pdf"
# comparemodshist(meanlist_true, meanlist_lin, "TP", "LP")
# savecomparehisttitlesigmatheta = "Comp_DP_LP2_pi5_003.pdf"
# savecomparehisttitleR= "Comp_DP_LP_R2_pi5_003.pdf"
# comparemodshist(meanlist_damped, meanlist_lin, "DP", "LP")
# 


# # #TIRSDAG: Comparison plots with different number of measurement points:
# #6 data points
# nummeasurements = 6
# df_lin_6 = pendulum(R_real, theta0, runtime, nummeasurements)
# df_true_6 = pendulum_true(theta0, runtime, h, m, g, L, I, nummeasurements)
# df_damped_6 = pendulum_damped(angledamped, timesdamped, nummeasurements)
# 
# #12 data points
# nummeasurements = 12
# df_lin_12 = pendulum(R_real, theta0, runtime, nummeasurements)
# df_true_12 = pendulum_true(theta0, runtime, h, m, g, L, I, nummeasurements)
# df_damped_12 = pendulum_damped(angledamped, timesdamped, nummeasurements)
# 
# 
# #20 data points
# nummeasurements = 20
# df_lin_20 = pendulum(R_real, theta0, runtime, nummeasurements)
# df_true_20 = pendulum_true(theta0, runtime, h, m, g, L, I, nummeasurements)
# df_damped_20 = pendulum_damped(angledamped, timesdamped, nummeasurements)
# 
# #Change to which you want to check
# df_6 = df_true_6
# df_12 = df_true_12
# df_20 = df_true_20
# 
# #6dp
# meanlist_6 = fit_model(df_6, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_6 = meanlist_6[3:4]
# tot_means_6 = maketot_means(meanlist_6)
# 
# # # Make histogram ready list
# meanlist_6$sample=1:nrow(meanlist_6)
# meanlist_freq_6 <- melt(meanlist_6, id.vars = "sample")
# 
# #12dp
# meanlist_12 = fit_model(df_12, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_12 = meanlist_12[3:4]
# tot_means_12 = maketot_means(meanlist_12)
# 
# # # Make histogram ready list
# meanlist_12$sample=1:nrow(meanlist_12)
# meanlist_freq_12 <- melt(meanlist_12, id.vars = "sample")
# 
# #20dp
# meanlist_20 = fit_model(df_20, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_20 = meanlist_20[3:4]
# tot_means_20 = maketot_means(meanlist_20)
# 
# # # Make histogram ready list
# meanlist_20$sample=1:nrow(meanlist_20)
# meanlist_freq_20 <- melt(meanlist_20, id.vars = "sample")
# 
# #SUMMARIES
# xtable(summary(meanlist_6))
# xtable(summary(meanlist_12))
# xtable(summary(meanlist_20))
# 
# #MAKE PLOTS
# plot_6 = plot_histogram_free_scale(meanlist_freq_6, tot_means_6, "pi10TP_histogram_6.pdf")
# plot_12 = plot_histogram_free_scale(meanlist_freq_12, tot_means_12, "pi10TP_histogram_12.pdf")
# plot_20 = plot_histogram_free_scale(meanlist_freq_20, tot_means_20, "pi10TP_histogram_20.pdf")
# 
# savecomparehisttitlesigmatheta = "pi10TP_Compsigmatheta_6_12.pdf"
# savecomparehisttitleR= "pi10TP_Comp_R_6_12.pdf"
# comparemodshist(meanlist_6, meanlist_12, "6 Points", "12 Points")
# 
# savecomparehisttitlesigmatheta = "pi10TP_Compsigmatheta_6_20.pdf"
# savecomparehisttitleR= "pi10TP_Comp_R_6_20.pdf"
# comparemodshist(meanlist_6, meanlist_20, "6 Points", "20 Points")
# 
# savecomparehisttitlesigmatheta = "pi10TP_Compsigmatheta_12_20.pdf"
# savecomparehisttitleR= "pi10TP_Comp_R_12_20.pdf"
# comparemodshist(meanlist_12, meanlist_20, "12 Points", "20 Points")


#ONSDAG: Priors:
#Choose which function
sigma_theta=0.03 #Angle measurement noise. BASIS = 0.03
prior_sigma_theta_mu=0.02 #Goes into STAN code. Mean of sigma_theta. BASIS = 0.02
prior_sigma_theta_var = 0.04 #Goes into STAN code. Variance of sigma_theta. BASIS = 0.04
range_sigma_theta_upper = 0.06 #Goes into STAN code. Upper value of sigma_theta. BASIS = 0.06

df = pendulum(R_real, theta0, runtime, 12)
df_true = pendulum_true(theta0, runtime, h, m, g, L, I, 12)
df_damped = pendulum_damped(angledamped, timesdamped, 12)
df = df_damped

#BASIS PRIORS
stanfile = 'Pendulum_model.stan'
meanlist = fit_model(df, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
meanlist = meanlist[3:4]
tot_means = maketot_means(meanlist)

# # Make histogram ready list
meanlist$sample=1:nrow(meanlist)
meanlist_freq <- melt(meanlist, id.vars = "sample")

#WEAK PRIORS
stanfile = 'Pendulum_model_weak.stan'
meanlist_weak = fit_model(df, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
meanlist_weak = meanlist_weak[3:4]
tot_means_weak = maketot_means(meanlist_weak)

# # Make histogram ready list
meanlist_weak$sample=1:nrow(meanlist_weak)
meanlist_freq_weak <- melt(meanlist_weak, id.vars = "sample")

#INFORMATIVE PRIORS
stanfile = 'Pendulum_model_informative.stan'
meanlist_inf = fit_model(df, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
meanlist_inf = meanlist_inf[3:4]
tot_means_inf = maketot_means(meanlist_inf)

# # Make histogram ready list
meanlist_inf$sample=1:nrow(meanlist_inf)
meanlist_freq_inf <- melt(meanlist_inf, id.vars = "sample")

#MAKE PLOTS CHANGE NAMES
plot_histogram_free_scale(meanlist_freq_weak, tot_means_weak, "DP_weak_histogram.pdf")
plot_histogram_free_scale(meanlist_freq_inf, tot_means_inf, "DP_inf_histogram.pdf")
plot_histogram_free_scale(meanlist_freq, tot_means, "DP_basis_histogram.pdf")

savecomparehisttitleR = "DP_compweakbasisR.pdf"
savecomparehisttitlesigmatheta = "DP_compweakbasissigmatheta.pdf"
comparemodshist(meanlist_weak, meanlist, "Weak", "Basis")

savecomparehisttitleR = "DP_compinfbasisR.pdf"
savecomparehisttitlesigmatheta = "DP_compinfbasissigmatheta.pdf"
comparemodshist(meanlist_inf, meanlist, "Informative", "Basis")

savecomparehisttitleR = "DP_compinfweakR.pdf"
savecomparehisttitlesigmatheta = "DP_compinfweaksigmatheta.pdf"
comparemodshist(meanlist_inf, meanlist_weak, "Informative", "Weak")


# xtable(summary(meanlist_weak))
# xtable(summary(meanlist_inf))








#KJØRING 2 og 3 (2: use df, 3: use df_true)
#RUN FUNCTIONS
# 
# df = df_damped
# 
# stanfile = 'Pendulum_model.stan'
# meanlist_lin = fit_model(df, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_lin = meanlist_lin[3:4]
# tot_means_lin = maketot_means(meanlist_lin)
# 
# # # Make histogram ready list
# meanlist_lin$sample=1:nrow(meanlist_lin) 
# meanlist_freq_lin <- melt(meanlist_lin, id.vars = "sample")
# 
# #RUN FUNCTIONS
# stanfile = 'Pendulum_model_weak.stan'
# meanlist_weak = fit_model(df, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_weak = meanlist_weak[3:4]
# tot_means_weak = maketot_means(meanlist_weak)
# 
# # # Make histogram ready list
# meanlist_weak$sample=1:nrow(meanlist_weak) 
# meanlist_freq_weak <- melt(meanlist_weak, id.vars = "sample")
# 
# #RUN FUNCTIONS
# stanfile = 'Pendulum_model_informative.stan'
# meanlist_inf = fit_model(df, seeds, theta0, L, sigma_theta, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# meanlist_inf = meanlist_inf[3:4]
# tot_means_inf = maketot_means(meanlist_inf)
# 
# # # Make histogram ready list
# meanlist_inf$sample=1:nrow(meanlist_inf)
# meanlist_freq_inf <- melt(meanlist_inf, id.vars = "sample")
# 
# #MAKE PLOTS
# plot_histogram_free_scale(meanlist_freq_weak, tot_means_weak, "TP_weak_histogram2.pdf")
# plot_histogram_free_scale(meanlist_freq_inf, tot_means_inf, "TP_inf_histogram2.pdf")
# 
# savecomparehisttitleR = "TP_compweakbasisR2.pdf"
# savecomparehisttitlesigmatheta = "TP_compweakbasissigmatheta2.pdf"
# comparemodshist(meanlist_weak, meanlist_lin, "Weak", "Basis")
# 
# savecomparehisttitleR = "TP_compinfbasisR2.pdf"
# savecomparehisttitlesigmatheta = "TP_compinfbasissigmatheta2.pdf"
# comparemodshist(meanlist_inf, meanlist_lin, "Informative", "Basis")
# 
# 
# #TP dev
# ci_dev = fit_model_extract_ci_dev(df, seeds, theta0, L, sigma_theta_list, prior_sigma_theta_mu, prior_sigma_theta_var, range_sigma_theta_upper)
# plotdev = plot_dev(sigma_theta_list, ci_dev$dev, "Dev_true.pdf")
# plotdev
# plotci = plot_ci(sigma_theta_list, ci_dev$ci_R, "Ci_true.pdf")
# plotci
# 
# xtable(summary(meanlist_weak))
# xtable(summary(meanlist_inf))
# 
