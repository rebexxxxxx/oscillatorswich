#3 switches, 3 oscillators
#plots
#demonstrating slow promoter
install.packages("xlsx")
install.packages("beepr")
library(xlsx)
library(beepr)

#connection and frequency changes 
#3 genes (switches)
#3 transcription factors (oscillators)
#1 cliques

#pull in connection matrices 
wd = c("/Users/rebeccasansale/Desktop/progenitor_development")

#timing of simulations
time = c(20,30)

nat_frequency_clique_t1 = as.matrix(c(0.5,0.5,0.5,0.02,0.02,0.02,0.02))

numosc_t1 = 7
numswitch_t1 = 5


source(sprintf("%s/simulateCoupled_progenitor.R", wd))
#simulateCoupled <- function(t0, tf, timesteps, delT, numOscillators_t1, numSwitches_t1,
#numOscillators_t2, numSwitches_t2,kxx, kxtheta, kthetax, kthetatheta,
#pa, pb, pc, pd, tau,naturalFrequency_t1,naturalFrequency_t2,working_directory)
t <- simulateCoupled(0, tf = 20,timesteps = .01, delT=.01, numosc_t1, numswitch_t1,
                           numosc_t2, numswitch_t2,
                           1,1,1,1,1,
                           nat_frequency_clique_t1, nat_frequency_clique_t2, working_directory = wd)
switches = data.frame(s1 = t$x[1,],
                      s2 = t$x[2,],
                      s3 = t$x[3,],
                      s4 = t$x[4,],
                      s5 = t$x[5,])
#switches = paste0(wd, "/figures/switches_time1_time2.pdf")
#pdf(file = switches)
longswitch <- gather(switches, key = switch_number, value = switch_energy)
longswitch$switch_number <- as.factor(longswitch$switch_number)
longswitch$time <- rep(seq(1,2001),5)

g <- ggplot(data = longswitch,
         aes(x=time, y = switch_energy, color = switch_number)) +
        geom_line(aes(color = switch_number)) +

g.animation<-g +
  transition_reveal(time)
animate(g.animation, height = 500, width = 800, fps = 30, duration = 10,
        end_pause = 60, res = 100)
anim_save("/Users/rebeccasansale/Desktop/progenitor_development/2switches.gif")

'
plot(x = seq(1,21,0.01), y = t$x[1,], type = 'l', ylim = c(-3,11), col= "blue", main = "Gene Expression, Time 1 and Time 2", ylab = "mRNA concentration", xlab = "time")
lines(x = seq(1,21,0.01), y = t$x[2,], type = 'l', col = "blue2")
lines(x = seq(1,21,0.01), y = t$x[3,], type = 'l', col = "red")
lines(x = seq(1,21,0.01), y = t$x[4,], type = 'l', col = "red2")
lines(x = seq(1,21,0.01), y = t$x[5,], type = 'l', col = "red3")
dev.off()
'

osc_location = data.frame(osc_location1 = t$raw_theta2[1,],
                          osc_location2 = t$raw_theta2[2,],
                          osc_location3 = t$raw_theta2[3,],
                          osc_location4 = t$raw_theta2[4,],
                          osc_location5 = t$raw_theta2[5,])





#same for oscillators
oscillators = data.frame(o1 = t$theta[1,],
                      o2 = t$theta[2,],
                      o3 = t$theta[3,],
                      o4 = t$theta[4,],
                      o5 = t$theta[5,],
                      o6 = t$theta[6,],
                      o7 = t$theta[7,])
#switches = paste0(wd, "/figures/switches_time1_time2.pdf")
#pdf(file = switches)
longosc <- gather(oscillators, key = osc_number, value = osc_energy)
longosc$osc_number <- as.factor(longosc$osc_number)
longosc$time <- rep(seq(1,2001),7)

g <- ggplot(data = longosc,
            aes(x=time, y =osc_energy, color = osc_number)) +
  geom_line(aes(color = osc_number)) 
  
  g.animation<-g +
  transition_reveal(time)
animate(g.animation, height = 500, width = 1800, fps = 30, duration = 20,
        end_pause = 60, res = 100)
anim_save("/Users/rebeccasansale/Desktop/progenitor_development/7osc.gif")





osc = paste0(wd, "/figures/osc_time1_time2.pdf")
pdf(file = osc)
plot(x = seq(1,21,0.01), y = t$theta[1,], type = 'l', ylim = c(0,6.5), col = "blue", main = "Transcription Factor Networks, Time 1 and Time 2", xlab = "TF Concentrations", ylab = "time")
lines(x = seq(1,21,0.01), y = t$theta[2,], type = 'l', col = "blue2")
lines(x = seq(1,21,0.01), y = t$theta[3,], type = 'l', col = "red")
lines(x = seq(1,21,0.01), y = t$theta[4,], type = 'l', col = "red2")
lines(x = seq(1,21,0.01), y = t$theta[5,], type = 'l', col = "red3")
lines(x = seq(1,21,0.01), y = t$theta[6,], type = 'l', col = "red")
lines(x = seq(1,21,0.01), y = t$theta[7,], type = 'l', col = "red2")
dev.off()


csv_file_name = paste0(wd, "/simulation_output_data_time1time2.xlsx")
write.xlsx(outputs$x, file = csv_file_name, sheetName = "switches")
write.xlsx(outputs$theta, file = csv_file_name, sheetName = "oscillators", append = TRUE)
write.xlsx(outputs$omega, file = csv_file_name, sheetName = "frequencies", append = TRUE)
write.xlsx(nat_frequency_clique_t1, file = csv_file_name, sheetName = "nat_freq", append = TRUE)






nat_frequency_clique_t2 = as.matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,0.2))
numosc_t2 = 20
numswitch_t2 = 15

sprintf("%s/simulateCoupled_progenitort2.R", wd)
source(sprintf("%s/simulateCoupled_progenitort2.R", wd))
#go into that file and change initial conditions to the following:
init_x <- outputs$x[,(ncol(outputs$x)-5)]
init_theta <- outputs$theta[,(ncol(outputs$theta)-5)]
init_omega <- outputs$omega[,(ncol(outputs$omega)-5)]

outputs_t2 <- simulateCoupledt2(t0 = 0, tf = 200,timesteps = .01, delT=.01,
                           numosc_t2, numswitch_t2,
                           1,1,1,1,1,
                          nat_frequency_clique_t2,
                          working_directory = wd,
                          initial_x = init_x, 
                          initial_theta = init_theta,
                          initial_omega = init_omega)


beep(sound = 7)

osc_colors <- c("orchid", "orchid1", "orchid2", "orchid3", "orchid4",
                "palegreen", "palegreen1", "palegreen2", "palegreen3", "palegreen4")
switch_colors <- c("orchid", "orchid1", "orchid2", "orchid3", "orchid4",
                   "palegreen", "palegreen1", "palegreen2", "palegreen3", "palegreen4")

#changes in phase
all_file_name = paste0(wd, "/figures/osc_switch_time1_time2.pdf")
pdf(file = all_file_name)



par(mfrow = c(2,2))
plot(x = seq(1,19,by = 0.01), y = outputs$theta[1,1:length(seq(1,19,by = 0.01))], col = osc_colors[1], type = 'l', ylim = c(0,8), xlim = c(0,21), main = c("Epithelial Stem Cell Cycle, Cell Cycle 1"), xlab = "Time", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[2,1:length(seq(1,19,by = 0.01))], col = osc_colors[2], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[3,1:length(seq(1,19,by = 0.01))], col = osc_colors[3], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[4,1:length(seq(1,19,by = 0.01))], col = osc_colors[4], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[5,1:length(seq(1,19,by = 0.01))], col = osc_colors[5], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)


plot(x = seq(1,19,by = 0.01), y = outputs_t2$theta[1,1:length(seq(1,19,by = 0.01))], col = osc_colors[1], type = 'l', ylim = c(0,8), xlim = c(0,20), main = c("Epithelial Stem Cell Cycle, Cell Cycle 2"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[2,1:length(seq(1,19,by = 0.01))], col = osc_colors[2], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[3,1:length(seq(1,19,by = 0.01))], col = osc_colors[3], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[4,1:length(seq(1,19,by = 0.01))], col = osc_colors[4], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[5,1:length(seq(1,19,by = 0.01))], col = osc_colors[5], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[6,1:length(seq(1,19,by = 0.01))], col = osc_colors[1], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[7,1:length(seq(1,19,by = 0.01))], col = osc_colors[2], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[8,1:length(seq(1,19,by = 0.01))], col = osc_colors[3], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[9,1:length(seq(1,19,by = 0.01))], col = osc_colors[4], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[10,1:length(seq(1,19,by = 0.01))], col = osc_colors[5], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[11,1:length(seq(1,19,by = 0.01))], col = osc_colors[1], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[15,1:length(seq(1,19,by = 0.01))], col = osc_colors[6], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[16,1:length(seq(1,19,by = 0.01))], col = osc_colors[7], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[17,1:length(seq(1,19,by = 0.01))], col = osc_colors[8], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)
lines(x = seq(1,19,by = 0.01), y = outputs_t2$theta[18,1:length(seq(1,19,by = 0.01))], col = osc_colors[9], type = 'l', ylim = c(0,7), xlim = c(0,100), main = c("Change in Phase of each Oscillator"), sub = c("4 Oscillators 10 switches"), xlab = "Time", ylab = "Transcription Factor Concentration", lwd = 4)


#setting font size for legend
#op <- par(cex = 0.75)
#legend(x = 0, y = 8, c("Oscillator 1", "Oscillator 2","Oscillator 3"), lty = c(1,1,1), lwd = c(4,4,4), col = osc_colors, ncol = 2)
#dev.off()

#changes in switches
#switch_file_name = paste0(wd, "/figures/switches.pdf")
#pdf(file = switch_file_name)

#time 1
plot(x = seq(1,19,by = 0.01), y = outputs$x[1,1:length(seq(1,19,by = 0.01))], col = switch_colors[1], type = 'l', ylim = c(0,7), xlim = c(0,21), main = c("mRNA Concentrations of Differentiated Cells"), xlab = "Time", ylab = "mRNA Concentration", lwd = 1)
lines(x = seq(1,19,by = 0.01), y = outputs$x[2,1:length(seq(1,19,by = 0.01))], col = switch_colors[2], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,19,by = 0.01), y = outputs$x[3,1:length(seq(1,19,by = 0.01))], col = switch_colors[3], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,19,by = 0.01), y = outputs$x[4,1:length(seq(1,19,by = 0.01))], col = switch_colors[4], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,19,by = 0.01), y = outputs$x[5,1:length(seq(1,19,by = 0.01))], col = switch_colors[5], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 1)

#time 2
plot(x = seq(1,29,by = 0.01), y = outputs_t2$x[1,1:length(seq(1,29,by = 0.01))], col = switch_colors[1], type = 'l', ylim = c(0,19), xlim = c(0,21), main = c("mRNA Concentrations of Differentiated Cells"), xlab = "Time", ylab = "mRNA Concentration", lwd = 1)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[2,1:length(seq(1,29,by = 0.01))], col = switch_colors[2], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[3,1:length(seq(1,29,by = 0.01))], col = switch_colors[3], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[4,1:length(seq(1,29,by = 0.01))], col = switch_colors[4], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[5,1:length(seq(1,29,by = 0.01))], col = switch_colors[5], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 1)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[6,1:length(seq(1,29,by = 0.01))], col = switch_colors[2], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[7,1:length(seq(1,29,by = 0.01))], col = switch_colors[3], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[8,1:length(seq(1,29,by = 0.01))], col = switch_colors[4], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[9,1:length(seq(1,29,by = 0.01))], col = switch_colors[5], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 1)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[10,1:length(seq(1,29,by = 0.01))], col = switch_colors[6], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[11,1:length(seq(1,29,by = 0.01))], col = switch_colors[7], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[12,1:length(seq(1,29,by = 0.01))], col = switch_colors[8], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[13,1:length(seq(1,29,by = 0.01))], col = switch_colors[9], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 1)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[14,1:length(seq(1,29,by = 0.01))], col = switch_colors[10], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 3)
lines(x = seq(1,29,by = 0.01), y = outputs_t2$x[15,1:length(seq(1,29,by = 0.01))], col = switch_colors[5], type = 'l', ylim = c(0,.5), xlim = c(0,100), main = c("Change in Switch"), sub = c("4 Osc 10 switches"), lwd = 1)



#op <- par(cex = 0.85)
#legend(x = 0, y = 9, c("Switch 1","Switch 2", "Switch 3"), lty = c(1,1,1), lwd = c(2,2,2), col = c(switch_colors[1], switch_colors[8]), ncol = 2)
dev.off()


#csv_file_name = paste0(wd, "/simulation_output_data_t2.xlsx")
#write.xlsx(outputs$x, file = csv_file_name, sheetName = "switches")
#write.xlsx(outputs$theta, file = csv_file_name, sheetName = "oscillators", append = TRUE)
#write.xlsx(outputs$omega, file = csv_file_name, sheetName = "frequencies", append = TRUE)


#combine times 1 and 2



