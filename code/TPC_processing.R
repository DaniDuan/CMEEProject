setwd('/home/danica/Documents/CMEEProject/data')
rm(list=ls())
graphics.off()
library(minpack.lm)
library(ggplot2)
library(reshape2)
library(grid)

k = 8.61*10^(-5) 

logSchoolfield = function(Temp, lnB0, T_pk, Ea, E_D){
  return(lnB0-Ea*(1/Temp-1/273.15)/k-log(1+(Ea/(E_D-Ea))*exp(E_D*(1/T_pk-1/Temp)/k)))
}

Schoolfield = function(Temp, B0, T_pk, Ea, E_D){
  return(B0*exp(-Ea*(1/Temp-1/273.15)/k)/(1+(Ea/(E_D-Ea))*exp(E_D*(1/T_pk-1/Temp)/k)))
}

CUE = function(G, R){
  return(G/(G+R))
}


data = read.csv("./aerobic_tpc_data.csv")
data = data[-which(data$temps_before_peak_resp<=3 | data$temps_before_peak_growth<=3),]
res_data = data.frame(data$E_resp, data$B0_resp, data$E_D_resp, data$Tpk_resp, data$Ppk_resp, data$r_sq_resp)
names(res_data) = c('Ea', 'B0', 'E_D', 'T_pk', 'P_pk','r_sq')
grow_data = data.frame(data$E_growth, data$B0_growth, data$E_D_growth, data$Tpk_growth, data$Ppk_growth, data$r_sq_growth)
names(grow_data) = c('Ea', 'B0', 'E_D', 'T_pk', 'P_pk','r_sq')

###### Growth ######
plot(1, type="n", xlab = "Temperature(celsius)", ylab = 'growth rate', xlim=c(0, 50), ylim = c(0, 0.8))
for(i in 1:nrow(grow_data)){
  lines(Schoolfield(Temp = 273.15:323.15, B0 = grow_data$B0[i], T_pk = grow_data$T_pk[i], Ea = grow_data$Ea[i], E_D = grow_data$E_D[i]), typ = 'l')
}

plot(1, type="n", xlab = "Temperature(celsius)", ylab = 'ln(growth rate)', xlim=c(0, 50), ylim=c(-15, 0))
for(i in 1:nrow(grow_data)){
  lines(logSchoolfield(Temp = 273.15:323.15, lnB0 = log(grow_data$B0[i]), T_pk = grow_data$T_pk[i], Ea = grow_data$Ea[i], E_D = grow_data$E_D[i]), typ = 'l')
}


####### Respiration ######
plot(1, type="n", xlab = "Temperature(celsius)", ylab = 'respiration rate', xlim=c(0, 50), ylim = c(0, 6))
for(i in 1:nrow(res_data)){
  lines(Schoolfield(Temp = 273.15:323.15, B0 = res_data$B0[i], T_pk = res_data$T_pk[i], Ea = res_data$Ea[i], E_D = res_data$E_D[i]), typ = 'l')
}

plot(1, type="n", xlab = "Temperature(celsius)", ylab = 'ln(respiration rate)', xlim=c(0, 50), ylim=c(-7, 2))

for(i in 1:nrow(res_data)){
  lines(logSchoolfield(Temp = 273.15:323.15, lnB0 = log(res_data$B0[i]), T_pk = res_data$T_pk[i], Ea = res_data$Ea[i], E_D = res_data$E_D[i]), typ = 'l')
}


G = Schoolfield(Temp = 273.15:323.15, B0 = grow_data$B0[i], T_pk = grow_data$T_pk[i], Ea = grow_data$Ea[i], E_D = grow_data$E_D[i])
R = Schoolfield(Temp = 273.15:323.15, B0 = res_data$B0[i], T_pk = res_data$T_pk[i], Ea = res_data$Ea[i], E_D = res_data$E_D[i])

plot(1, type="n", xlab = "Temperature(celsius)", ylab = 'CUE', xlim=c(0, 50), ylim=c(0, 1))
for(i in 1:nrow(grow_data)){
  lines(CUE(G = Schoolfield(Temp = 273.15:323.15, B0 = grow_data$B0[i], T_pk = grow_data$T_pk[i], Ea = grow_data$Ea[i], E_D = grow_data$E_D[i]), R = Schoolfield(Temp = 273.15:323.15, B0 = res_data$B0[i], T_pk = res_data$T_pk[i], Ea = res_data$Ea[i], E_D = res_data$E_D[i])), typ = 'l')
}

CUE_Values= data.frame()
for(i in 1:nrow(grow_data)){
  CUE_Values = rbind(CUE_Values, CUE(G = Schoolfield(Temp = 273.15, B0 = grow_data$B0[i], T_pk = grow_data$T_pk[i], Ea = grow_data$Ea[i], E_D = grow_data$E_D[i]), R = Schoolfield(Temp = 273.15:323.15, B0 = res_data$B0[i], T_pk = res_data$T_pk[i], Ea = res_data$Ea[i], E_D = res_data$E_D[i])))
}
CUE(G = Schoolfield(Temp = 273.15:323.15, B0 = grow_data$B0[i], T_pk = grow_data$T_pk[i], Ea = grow_data$Ea[i], E_D = grow_data$E_D[i]), R = Schoolfield(Temp = 273.15:323.15, B0 = res_data$B0[i], T_pk = res_data$T_pk[i], Ea = res_data$Ea[i], E_D = res_data$E_D[i]))

CUE_0 = grow_data$B0/(grow_data$B0+res_data$B0)
hist(CUE_0)
mean(CUE_0)
mean(CUE_0) - sd(CUE_0)/sqrt(nrow(grow_data))
mean(CUE_0) + sd(CUE_0)/sqrt(nrow(grow_data))

###### Biotrait Growth ######

database = read.csv("./database.csv")
database = database[which(database$StandardisedTraitName == "Specific Growth Rate"),]
database[which(database$AmbientTemp > 150),]$AmbientTemp = database[which(database$AmbientTemp > 150),]$AmbientTemp - 273.15
database = database[which(database$ConKingdom == 'Bacteria'),]
scatter.smooth(database$AmbientTemp, database$StandardisedTraitValue, xlab ='Temperature', ylab = 'Growth rate')

summ = read.csv("./summary.csv")
summ = summ[-which(summ$Points_Before_Peak < 3),]
summ = summ[which(summ$Trait == "Specific Growth Rate"),]
summ = summ[-which(summ$AIC > 100),]
summ = summ[-which(summ$E_D < summ$E),]
summ = summ[-which(summ$T_pk > 370),]
summ = summ[-which(summ$E_D > 30),]
summ = summ[which(summ$ConKingdom == 'Bacteria'),]
#summ = summ[which(summ$B0>0),]
s_data = data.frame(summ$X, summ$Species, summ$B0, summ$E, summ$T_pk, summ$E_D, summ$R_Squared, summ$AIC, summ$BIC)
names(s_data) = c('ID', 'Species', 'B0', 'Ea', 'T_pk', 'E_D', 'r_sq', 'AIC', 'BIC')

plot(1, type="n", xlab = "Temperature (celsius)", ylab = 'Growth Rate', xlim=c(0,100), ylim = c(0, 90))
for(i in 1:nrow(s_data)){
  lines(x = 0:100, y = Schoolfield(Temp = 273.15:373.15, B0 = exp(s_data$B0[i]), T_pk = s_data$T_pk[i], Ea = s_data$Ea[i], E_D = s_data$E_D[i]), typ = 'l')
}

png('BioTrait_growth.png', width = 800, height = 600)
par(mar=c(7,7,7,1)+.1)
plot(1, type="n", xlab = "Temperature (celsius)", ylab = 'Growth Rate (1/Time)', xlim=c(0,60), ylim = c(0, 90), cex.lab=2.5, cex.axis=2.5)
for(i in 1:nrow(s_data)){
  lines(x = 0:60, y = Schoolfield(Temp = 273.15:333.15, B0 = exp(s_data$B0[i]), T_pk = s_data$T_pk[i], Ea = s_data$Ea[i], E_D = s_data$E_D[i]), typ = 'l', col = alpha('black',0.5), lwd = 2)
}
text(-8,105, 'A', cex = 3, xpd = NA)
dev.off()

plot(1, type="n", xlab = "Temperature (celsius)", ylab = 'ln(Growth Rate)', xlim=c(0, 100), ylim=c(-5, 3.5))
for(i in 1:nrow(s_data)){
  lines(logSchoolfield(Temp = 273.15:373.15, lnB0 = exp(s_data$B0[i]), T_pk = s_data$T_pk[i], Ea = s_data$Ea[i], E_D = s_data$E_D[i]), typ = 'l')
}

hist(exp(s_data$B0))
max(exp(s_data$B0))
mean(exp(s_data$B0))
