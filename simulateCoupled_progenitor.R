#####
#playing with the natural frequency of the oscillators
#originally, natural frequency was a single value input
#making it a vector now





#RSansale 3/22/20
#simulateCoupled function

#need to return a few different objects:
#time
#x (matrix)
#theta (matrix)
#omega (matrix)
#A (matrix)
#B (matrix)
#C (matrix)
#D (matrix)

#inputs
#t0, tf, delT, numOscillators, numSwitches, kxx, kxtheta, kthetax, kthetatheta
#pa, pb, pc, pd, tau, naturalFrequency
#0, tMax,0.01, 100, 100, kxxOn,kOtherOn,kOn,kttOn,1,1,1,1,1,omegaHatOn)


#naturalFrequency_clique <- t(c(0,0,0,0,0, 1,1,1,1,1))

simulateCoupled <- function(t0, tf, timesteps, delT, numOscillators_t1, numSwitches_t1,
                            numOscillators_t2, numSwitches_t2,
                            pa, pb, pc, pd, tau,
                            naturalFrequency_t1,
                            naturalFrequency_t2,
                            working_directory){
  
  #initialize the variables
  time_t1  = seq(t0, tf[1], by = delT)
  time_length = length(time_t1)
  x     = matrix(0, nrow = numSwitches_t1, ncol = length(time_t1))
  theta = matrix(0, numOscillators_t1, length(time_t1))
  omega = matrix(0, numOscillators_t1, length(time_t1))
  
  #set the initial conditions
  #x[, 1]    = rnorm(n = numSwitches, mean = 0, sd = 1)
  x[,1]   = 0
  theta[,1] = 2*pi*runif(n = numOscillators_t1, min = 0, max = 1)
  omega[,1] = naturalFrequency_t1
  
  A = matrix(1, nrow = numSwitches_t1, ncol = numSwitches_t1)
  diag(A) <- 0
  
  B = matrix(1, nrow = numOscillators_t1, ncol = numOscillators_t1)
  #C is the matrix involved in the oscillator equation
  C = matrix(1, nrow = numOscillators_t1, ncol = numOscillators_t1)
  diag(C) <- 0 
  
  D = matrix(1, nrow = numOscillators_t1, ncol = numSwitches_t1)
  
  mod <- function(x,m){
    tl <- floor(x/m)
    return(x-tl*m)
  }
  
  # run the simulation
  #sprintf("%s/computeModelDerivative_progenitor.R", working_directory)
  source(sprintf("%s/computeModelDerivative_progenitor.R", working_directory))
#for(timestep in 1:timesteps){
  #if(timestep==1){naturalFrequency = naturalFrequency_t1}else{naturalFrequency = naturalFrequency_t2}
  for(iTime in 2:length(time_t1)){
    #computeModelDerivative <- function(x, theta, omega, natural_omega, tau, cutofftheta, kxx_parameter, kthetatheta_parameter, kthetax_parameter, 
    #                                   kxtheta_parameter, A, B, C, D, current_time, working_directory){
    tt = iTime
    output = computeModelDerivative(as.matrix(x[,iTime-1]), as.matrix(theta[,iTime-1]), as.matrix(omega[,iTime-1]),
                                    naturalFrequency_t1, 0, pi,
                                    A, B, C, D,
                                    working_directory,
                                    tt)

    x[, iTime]     = x[,iTime-1]    + output$dx%*%delT
    tchange = output$dtheta*delT
    ttchange = theta[,iTime-1] + tchange
    theta[,iTime] = mod((ttchange),2*pi)
    omega[,iTime] = omega[,iTime-1]  + output$domega*delT
    iTime = iTime + 1
    print(c(iTime))
    if(iTime==length(time_t1)){
    newlist <- list("x" = x, "theta" = theta, "omega" = omega, "raw_theta" = tchange,"raw_theta2" = ttchange)
    return(newlist)}
  }

}



#  x_t1 = x
#  theta_t1 = theta
#  omega_t1 = omega

#set the initial conditions for next time step
#empty switch matrix for just the new switches
'  x = matrix(0, nrow = numSwitches_t2, ncol = length(time_t2))
  x[,1]   = 0
  theta = matrix(0, numOscillators_t2, length(time_t2))
  theta[,1] = theta_t1[,iTime-2]
  omega = matrix(0, numOscillators_t2, length(time_t2))
  omega[,1] = naturalFrequency_t2
  iTime = NULL
  naturalFrequency = naturalFrequency_t2
  A = matrix(1, nrow = numSwitches_t2, ncol = numSwitches_t2)
  diag(A) <- 0
  
  B = matrix(1, nrow = numOscillators_t2, ncol = numOscillators_t2)
  #C is the matrix involved in the oscillator equation
  C = matrix(1, nrow = numOscillators_t2, ncol = numOscillators_t2)
  diag(C) <- 0 
  
  D = matrix(1, nrow = numOscillators_t2, ncol = numSwitches_t2)
  '
#}
#newlist <- list("x_t1" = x_t1, "x_t2" = x,
#  "omega_t1" = omega_t1, "omega_t2" = omega,
#  "theta_t1" = theta_t1, "theta_t2" = theta)  

