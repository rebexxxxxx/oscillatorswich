#Rsansale
#computeModelDerivative
#inputs: x, theta, omega,natural_omega, tau, cutofftheta, kxx, kthetatheta, kthetax, kxtheta, A, B, C, D
#outputs: dx, dtheta, domega 

computeModelDerivative <- function(x, theta, omega, natural_omega, tau, cutofftheta, 
                                    A, B, C, D, working_directory, total_time){
  
  numSwitch = nrow(x)
  numOscillator = nrow(theta)
  
  # make other variables the correct dimension
  if (length(tau) == 1){
    tau2 = data.frame(rep(0, times = numSwitch))
    colnames(tau2) <- c("values")
    tau2$values = rep((tau*1), times = numSwitch)
    tau = tau2$values
    tau2 = NULL
  }
  
  sprintf("%s/getBinaryValue_progenitor.R", working_directory)
  source(sprintf("%s/getBinaryValue_progenitor.R", working_directory))
  # get binary values for x and theta
  #lets make the cut off value different for each switch
  #c(off,on) for the last term
  binX     = matrix(getBinaryValue(x, 0, c(0,1)))
  binTheta = getBinaryValue(theta, cutofftheta, c(0,1))

  #switch to switch connection
  #pull in connection matrices 
  'kxx_matrix = as.matrix(read.csv(paste0(working_directory, "/switch-to-switch_t1.csv"), header = FALSE),nrow = numSwitch, ncol = numSwitch)
  kxx_matrix[upper.tri(kxx_matrix)] <- t(kxx_matrix)[upper.tri(kxx_matrix)]
  kxx_matrix[1,] <- -kxx_matrix[1,]
  kxx_matrix[,1] <- -kxx_matrix[,1]
  diag(kxx_matrix) <- 1
  kxt_matrix = as.matrix(read.csv(paste0(working_directory, "/switch-oscillator_t1.csv"), header = FALSE), nrow = numSwitch, ncol = numOscillator)
  kxt_matrix[upper.tri(kxt_matrix)] <- t(kxt_matrix)[upper.tri(kxt_matrix)]/10
  kxt_matrix[lower.tri(kxt_matrix)] <- kxt_matrix[lower.tri(kxt_matrix)]/5
  ktx_matrix = as.matrix(t(read.csv(paste0(working_directory, "/switch-oscillator_t1.csv"), header = FALSE)),nrow = numOscillator, ncol = numSwitch)
  ktx_matrix[upper.tri(ktx_matrix)] <- t(ktx_matrix)[upper.tri(ktx_matrix)]/10
  ktx_matrix[lower.tri(ktx_matrix)] <- ktx_matrix[lower.tri(ktx_matrix)]/5
  ktt_matrix = as.matrix(read.csv(paste0(working_directory, "/oo_t1.csv"), header = FALSE))
  ktt_matrix[upper.tri(ktt_matrix)] <- t(ktt_matrix)[upper.tri(ktt_matrix)]
'
  if(total_time < 1000){
  kxx_matrix = as.matrix(read.csv(paste0(working_directory, "/switch-to-switch_time1.csv"), header = FALSE),nrow = numSwitch, ncol = numSwitch)
  kxt_matrix = as.matrix(read.csv(paste0(working_directory, "/switch-oscillator_time1.csv"), header = FALSE), nrow = numSwitch, ncol = numOscillator)
  diag(kxt_matrix) <-1
  ktx_matrix = as.matrix(t(read.csv(paste0(working_directory, "/switch-oscillator_time1.csv"), header = FALSE)),nrow = numOscillator, ncol = numSwitch)
  diag(ktx_matrix) <-1
  ktt_matrix = as.matrix(read.csv(paste0(working_directory, "/oo_time1.csv"), header = FALSE))
  }else{
  kxx_matrix = as.matrix(read.csv(paste0(working_directory, "/switch-to-switch_time2.csv"), header = FALSE),nrow = numSwitch, ncol = numSwitch)
  kxt_matrix = as.matrix(read.csv(paste0(working_directory, "/switch-oscillator_time2.csv"), header = FALSE), nrow = numSwitch, ncol = numOscillator)
  diag(kxt_matrix) <-1
  ktx_matrix = as.matrix(t(read.csv(paste0(working_directory, "/switch-oscillator_time2.csv"), header = FALSE)),nrow = numOscillator, ncol = numSwitch)
  diag(ktx_matrix) <-1
  ktt_matrix = as.matrix(read.csv(paste0(working_directory, "/oo_time2.csv"), header = FALSE))
  }
  #change in switches
  dx      = -x + kxx_matrix%*%(binX) + rowSums(kxt_matrix%*%(binTheta))
  # Creates the square matrix sin(\theta_j - \theta_i)
  thetasq = matrix(rep(theta, times = numOscillator), nrow = numOscillator, ncol = numOscillator)
  delta = sin(thetasq-t(thetasq))
  # The differential equation for \theta
  dtheta  = omega + ktt_matrix%*%diag(C%*%delta)

  # The rectangular matrix \tilde{x}_j \hat{\omega}_\mu - \omega_\mu
  onevec  = matrix(rep(1, times = length(x)), nrow = numSwitch, ncol = 1)
  delom   = (binX %*% t(natural_omega)) - onevec %*% t(omega)
  # The differential equation for \omega
  domega  = rowSums(ktx_matrix%*%delom)
  
  newList <- list("dx" = dx, "dtheta" = dtheta, "domega" = domega)
  return(newList)
}


