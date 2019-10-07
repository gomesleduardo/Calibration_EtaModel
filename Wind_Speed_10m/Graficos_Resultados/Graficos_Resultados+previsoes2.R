# --------------------------------------------------------- #
# |     UNIVERSIDADE FEDERAL DO RIO DE JANEIRO - UFRJ     | #
# |               INSTITUTO DE MATEMATICA - IM            | #
# |      DEPARTAMENTO DE METODOS ESTATISTICOS - DME       | #
# |-------------------------------------------------------| #
# |            ORIENTADORA: THAIS C. O. FONSECA           | #
# |           ORIENTADORA: KELLY C. M. GONCALVES          | #
# |              ALUNO: LUIZ EDUARDO S. GOMES             | #
# |-------------------------------------------------------| #
# |       CALIBRACAO DINAMICA DE PREVISOES NUMERICAS      | #
# |               CRITERIOS DE COMPARACAO                 | #
# |-------------------------------------------------------| #

# PACOTES ----
pacotes<-c("mvtnorm","forecast","coda","matrixStats",
           "corpcor","truncnorm","Matrix","akima",
           "Rcpp","Rcpp11","RcppArmadillo",
           "geoR","adaptMCMC","tmvtnorm",
           "sp","colorRamps","RColorBrewer")

for (i in 1:length(pacotes))
{
  if (length(names(installed.packages()[,1])[names(installed.packages()[,1])==pacotes[i]])==0)
  {install.packages(pacotes[i], repos="http://cran.fiocruz.br/")}
  library(pacotes[i],character.only = TRUE)
}
rm(i,pacotes)

# GRAFICO - CADEIAS (diario) ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
# load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
#             datas[data.c],"STEMOS2.Rda"))
# sampMH2 <- sampMH

load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"SEMOS.Rda"))

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))

#------ b0 ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,1]),
       sampMH$samples[20000:42970,1]),
     type="l",ylab=expression(beta[0]),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,1],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,1],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ b1 ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,2]),
       sampMH$samples[20000:42970,2]),
     type="l",ylab=expression(beta[1]),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,2],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,1],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ phi ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,3]),
       sampMH$samples[20000:42970,3]),
     type="l",ylab=expression(phi),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,3],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,3],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ lambda ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,4]),
       sampMH$samples[20000:42970,4]),
     type="l",ylab=expression(lambda),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,4],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,4],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)




par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))

#------ theta0 ------#
d <- density(x = Thetat.chain.final[1,1,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[0]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,1,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,1,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,1,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,1,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta1 ------#
d <- density(x = Thetat.chain.final[1,2,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[1]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,2,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,2,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,2,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,2,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta2 ------#
d <- density(x = Thetat.chain.final[1,3,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[2]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,3,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,3,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,3,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,3,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta3 ------#
d <- density(x = Thetat.chain.final[1,4,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[3]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,4,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,4,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,4,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,4,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta4 ------#
d <- density(x = Thetat.chain.final[1,5,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[4]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,5,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,5,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,5,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,5,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta5 ------#
d <- density(x = Thetat.chain.final[1,6,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[5]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,6,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,6,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,6,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,6,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)




sigma.b0.prior <-  10
sigma.b1.prior <-  10
sigma.lambda.prior <-  10

#------ b0 ------#
par(mar=c(4.5,4.5,.5,.5))
acf(PSI.chain$b0,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$b0, adjust = 5)
plot(d, xlab=expression(beta[0]), main="", ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.b0.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ b1 ------#
acf(PSI.chain$b1,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$b1, adjust = 5,from = 0)
plot(d, xlab=expression(beta[1]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.b1.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ phi ------#
acf(PSI.chain$phi,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$phi, adjust = 5)
plot(d, xlab=expression(phi), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ lambda ------#
acf(PSI.chain$lambda,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$lambda, adjust = 5)
plot(d, xlab=expression(lambda), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

rm(list=ls())








datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
# load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
#             datas[data.c],"STEMOS2.Rda"))
# sampMH2 <- sampMH

load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"DGOP.Rda"))

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))

#------ phi ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,1]),
       sampMH$samples[20000:54950,1]),
     type="l",ylab=expression(phi),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,1],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)

#------ lambda ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,2]),
       sampMH$samples[20000:54950,2]),
     type="l",ylab=expression(lambda),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,2],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)

#------ Thetat ------#
par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))
for (m in (1:r))
{
  plot(1,type="n",
       ylim=c(min(graf.Thetat[3,m,-(1:15)],na.rm=T),max(graf.Thetat[1,m,-(1:15)],na.rm=T)),
       xlim=c(1,T),ylab=bquote(theta[.(m-1)]),xlab="Tempo",cex.lab=1.5,cex.axis=1.5)
  polygon(c(rev(1:T),1:T),c(rev(graf.Thetat[3,m,-1]),graf.Thetat[1,m,-1]),col="grey80",border=NA)
  lines(graf.Thetat[2,m,-1],lwd=2)
  abline(h=0,lty=3,lwd=2)
}

sigma.lambda.prior <-  10

#------ phi ------#
acf(PSI.chain$phi,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$phi, adjust = 5)
plot(d, xlab=expression(phi), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ lambda ------#
acf(PSI.chain$lambda,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$lambda, adjust = 5)
plot(d, xlab=expression(lambda), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)


rm(list=ls())








datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
# load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
#             datas[data.c],"STEMOS2.Rda"))
# sampMH2 <- sampMH

load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"STEMOS.Rda"))

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))

#------ b0 ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,1]),
       sampMH$samples[20000:54950,1]),
     type="l",ylab=expression(beta[0]),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,1],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,1],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ b1 ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,2]),
       sampMH$samples[20000:54950,2]),
     type="l",ylab=expression(beta[1]),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,2],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,1],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ phi ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,3]),
       sampMH$samples[20000:54950,3]),
     type="l",ylab=expression(phi),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,3],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,3],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ lambda ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,4]),
       sampMH$samples[20000:54950,4]),
     type="l",ylab=expression(lambda),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,4],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,4],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ Thetat ------#
par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))
for (m in (1:r))
{
  plot(1,type="n",
       ylim=c(min(graf.Thetat[3,m,-(1:15)],na.rm=T),max(graf.Thetat[1,m,-(1:15)],na.rm=T)),
       xlim=c(1,T),ylab=bquote(theta[.(m-1)]),xlab="Tempo",cex.lab=1.5,cex.axis=1.5)
  polygon(c(rev(1:T),1:T),c(rev(graf.Thetat[3,m,-1]),graf.Thetat[1,m,-1]),col="grey80",border=NA)
  lines(graf.Thetat[2,m,-1],lwd=2)
  abline(h=0,lty=3,lwd=2)
}

sigma.b0.prior <-  10
sigma.b1.prior <-  10
sigma.lambda.prior <-  10

#------ b0 ------#
par(mar=c(4.5,4.5,.5,.5))
acf(PSI.chain$b0,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$b0, adjust = 5)
plot(d, xlab=expression(beta[0]), main="", ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.b0.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ b1 ------#
acf(PSI.chain$b1,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$b1, adjust = 5,from = 0)
plot(d, xlab=expression(beta[1]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.b1.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ phi ------#
acf(PSI.chain$phi,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$phi, adjust = 5)
plot(d, xlab=expression(phi), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ lambda ------#
acf(PSI.chain$lambda,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$lambda, adjust = 5)
plot(d, xlab=expression(lambda), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

BETA0 <- data.frame(NULL)
BETA1 <- data.frame(NULL)
PHI <- data.frame(NULL)
LAMBDA <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
              datas[data.c],"STEMOS.Rda"))
  BETA0[1:2000,data.c] <- PSI.chain$b0
  BETA1[1:2000,data.c] <- PSI.chain$b1
  PHI[1:2000,data.c] <- PSI.chain$phi
  LAMBDA[1:2000,data.c] <- PSI.chain$lambda
}

names(BETA0) <- c("Dez","Jan","Fev","Mar","Abr","Mai",
                  "Jun","Jul","Ago","Set","Out","Nov")
names(BETA1) <- c("Dez","Jan","Fev","Mar","Abr","Mai",
                  "Jun","Jul","Ago","Set","Out","Nov")
names(PHI) <- c("Dez","Jan","Fev","Mar","Abr","Mai",
                  "Jun","Jul","Ago","Set","Out","Nov")
names(LAMBDA) <- c("Dez","Jan","Fev","Mar","Abr","Mai",
                  "Jun","Jul","Ago","Set","Out","Nov")

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))
A <- 
boxplot(BETA0,outline = F,col="grey80",
        ylab=expression(beta[0]),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,plot = F)

A$stats[1,] <- colQuantiles(as.matrix(BETA0),probs = .025)
A$stats[2,] <- colQuantiles(as.matrix(BETA0),probs = .025)
A$stats[3,] <- colQuantiles(as.matrix(BETA0),probs = .5)
A$stats[4,] <- colQuantiles(as.matrix(BETA0),probs = .975)
A$stats[5,] <- colQuantiles(as.matrix(BETA0),probs = .975)
boxplot(A$stats,outline = F,names=c("Dez","Jan","Fev","Mar","Abr","Mai",
                                    "Jun","Jul","Ago","Set","Out","Nov"),
        ylab=expression(beta[0]),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,ylim=c(min(A$stats),1.2))

for (i in 1:12)
{
  segments(x0 = i,y0 = A$stats[2,i],x1 = i,y1 = A$stats[4,i],lty=2)
  segments(x0 = i+.4,y0 = A$stats[2,i],x1 = i+.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[2,i],x1 = i-.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[3,i],x1 = i+.4,y1 = A$stats[3,i],col="white",lwd=3)
  segments(x0 = i-.4,y0 = A$stats[1,i],x1 = i+.4,y1 = A$stats[1,i],lwd=3)
  segments(x0 = i-.4,y0 = A$stats[4,i],x1 = i+.4,y1 = A$stats[4,i],lwd=3)
  points(i,A$stats[3,i],pch=4,lwd=3,col=4)
}

abline(v=seq(0.5,12.5,3),col="grey80")
text(2,1.2,"Verão",cex=1.5)
text(5,1.2,"Outono",cex=1.5)
text(8,1.2,"Inverno",cex=1.5)
text(11,1.2,"Primavera",cex=1.5)


boxplot(BETA1,outline = F,col="grey80",
        ylab=expression(beta[1]),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5)

A <- 
  boxplot(BETA1,outline = F,col="grey80",
          ylab=expression(beta[0]),xlab="Meses",
          cex.lab=1.5,cex.axis=1.5,plot = F)

A$stats[1,] <- colQuantiles(as.matrix(BETA1),probs = .025)
A$stats[2,] <- colQuantiles(as.matrix(BETA1),probs = .025)
A$stats[3,] <- colQuantiles(as.matrix(BETA1),probs = .5)
A$stats[4,] <- colQuantiles(as.matrix(BETA1),probs = .975)
A$stats[5,] <- colQuantiles(as.matrix(BETA1),probs = .975)
boxplot(A$stats,outline = F,names=c("Dez","Jan","Fev","Mar","Abr","Mai",
                                    "Jun","Jul","Ago","Set","Out","Nov"),
        ylab=expression(beta[1]),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,ylim=c(min(A$stats),0.3))

for (i in 1:12)
{
  segments(x0 = i,y0 = A$stats[2,i],x1 = i,y1 = A$stats[4,i],lty=2)
  segments(x0 = i+.4,y0 = A$stats[2,i],x1 = i+.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[2,i],x1 = i-.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[3,i],x1 = i+.4,y1 = A$stats[3,i],col="white",lwd=3)
  segments(x0 = i-.4,y0 = A$stats[1,i],x1 = i+.4,y1 = A$stats[1,i],lwd=3)
  segments(x0 = i-.4,y0 = A$stats[4,i],x1 = i+.4,y1 = A$stats[4,i],lwd=3)
  points(i,A$stats[3,i],pch=4,lwd=3,col=4)
}

abline(v=seq(0.5,12.5,3),col="grey80")
text(2,0.3,"Verão",cex=1.5)
text(5,0.3,"Outono",cex=1.5)
text(8,0.3,"Inverno",cex=1.5)
text(11,0.3,"Primavera",cex=1.5)

boxplot(PHI,outline = F,col="grey80",
        ylab=expression(phi),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5)

A <- 
  boxplot(PHI,outline = F,col="grey80",
          ylab=expression(phi),xlab="Meses",
          cex.lab=1.5,cex.axis=1.5,plot = F)

A$stats[1,] <- colQuantiles(as.matrix(PHI),probs = .025)
A$stats[2,] <- colQuantiles(as.matrix(PHI),probs = .025)
A$stats[3,] <- colQuantiles(as.matrix(PHI),probs = .5)
A$stats[4,] <- colQuantiles(as.matrix(PHI),probs = .975)
A$stats[5,] <- colQuantiles(as.matrix(PHI),probs = .975)
boxplot(A$stats,outline = F,names=c("Dez","Jan","Fev","Mar","Abr","Mai",
                                    "Jun","Jul","Ago","Set","Out","Nov"),
        ylab=expression(phi),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,ylim=c(min(A$stats),3.5))

for (i in 1:12)
{
  segments(x0 = i,y0 = A$stats[2,i],x1 = i,y1 = A$stats[4,i],lty=2)
  segments(x0 = i+.4,y0 = A$stats[2,i],x1 = i+.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[2,i],x1 = i-.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[3,i],x1 = i+.4,y1 = A$stats[3,i],col="white",lwd=3)
  segments(x0 = i-.4,y0 = A$stats[1,i],x1 = i+.4,y1 = A$stats[1,i],lwd=3)
  segments(x0 = i-.4,y0 = A$stats[4,i],x1 = i+.4,y1 = A$stats[4,i],lwd=3)
  points(i,A$stats[3,i],pch=4,lwd=3,col=4)
}

abline(v=seq(0.5,12.5,3),col="grey80")
text(2,3.5,"Verão",cex=1.5)
text(5,3.5,"Outono",cex=1.5)
text(8,3.5,"Inverno",cex=1.5)
text(11,3.5,"Primavera",cex=1.5)

boxplot(LAMBDA,outline = F,col="grey80",
        ylab=expression(lambda),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5)

A <- 
  boxplot(LAMBDA,outline = F,col="grey80",
          ylab=expression(lambda),xlab="Meses",
          cex.lab=1.5,cex.axis=1.5,plot = F)

A$stats[1,] <- colQuantiles(as.matrix(LAMBDA),probs = .025)
A$stats[2,] <- colQuantiles(as.matrix(LAMBDA),probs = .025)
A$stats[3,] <- colQuantiles(as.matrix(LAMBDA),probs = .5)
A$stats[4,] <- colQuantiles(as.matrix(LAMBDA),probs = .975)
A$stats[5,] <- colQuantiles(as.matrix(LAMBDA),probs = .975)
boxplot(A$stats,outline = F,names=c("Dez","Jan","Fev","Mar","Abr","Mai",
                                    "Jun","Jul","Ago","Set","Out","Nov"),
        ylab=expression(lambda),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,ylim=c(min(A$stats),.65))

for (i in 1:12)
{
  segments(x0 = i,y0 = A$stats[2,i],x1 = i,y1 = A$stats[4,i],lty=2)
  segments(x0 = i+.4,y0 = A$stats[2,i],x1 = i+.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[2,i],x1 = i-.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[3,i],x1 = i+.4,y1 = A$stats[3,i],col="white",lwd=3)
  segments(x0 = i-.4,y0 = A$stats[1,i],x1 = i+.4,y1 = A$stats[1,i],lwd=3)
  segments(x0 = i-.4,y0 = A$stats[4,i],x1 = i+.4,y1 = A$stats[4,i],lwd=3)
  points(i,A$stats[3,i],pch=4,lwd=3,col=4)
}

abline(v=seq(0.5,12.5,3),col="grey80")
text(2,.65,"Verão",cex=1.5)
text(5,.65,"Outono",cex=1.5)
text(8,.65,"Inverno",cex=1.5)
text(11,.65,"Primavera",cex=1.5)


rm(list=ls())



datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
# load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
#             datas[data.c],"STEMOS2.Rda"))
# sampMH2 <- sampMH

load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))

#------ b0 ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:50000,1]),
       sampMH$samples[20000:56000,1]),
     type="l",ylab=expression(beta[0]),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,1],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,1],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ b1 ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:50000,2]),
       sampMH$samples[20000:56000,2]),
     type="l",ylab=expression(beta[1]),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,2],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,2],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ phi ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[30000:60000,3]),
       sampMH$samples[20000:56000,3]),
     type="l",ylab=expression(phi),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,3],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,3],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ lambda ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:50000,4]),
       sampMH$samples[20000:56000,4]),
     type="l",ylab=expression(lambda),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,4],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,4],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))

#------ theta0 ------#
d <- density(x = Thetat.chain.final[1,1,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[0]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,1,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,1,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,1,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,1,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta1 ------#
d <- density(x = Thetat.chain.final[1,2,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[1]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,2,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,2,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,2,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,2,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta2 ------#
d <- density(x = Thetat.chain.final[1,3,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[2]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,3,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,3,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,3,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,3,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta3 ------#
d <- density(x = Thetat.chain.final[1,4,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[3]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,4,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,4,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,4,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,4,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta4 ------#
d <- density(x = Thetat.chain.final[1,5,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[4]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,5,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,5,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,5,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,5,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ theta5 ------#
d <- density(x = Thetat.chain.final[1,6,-(1:burnin)], adjust = 5)
plot(d, xlab=expression(theta[5]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(Thetat.chain.final[1,6,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,6,-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(Thetat.chain.final[1,6,-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(Thetat.chain.final[1,6,-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

sigma.b0.prior <-  10
sigma.b1.prior <-  10
sigma.lambda.prior <-  10


#------ b0 ------#
par(mar=c(4.5,4.5,.5,.5))
acf(PSI.chain$b0,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$b0, adjust = 5)
plot(d, xlab=expression(beta[0]), main="", ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.b0.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ b1 ------#
acf(PSI.chain$b1,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$b1, adjust = 5,from = 0)
plot(d, xlab=expression(beta[1]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.b1.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ phi ------#
acf(PSI.chain$phi,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$phi, adjust = 5)
plot(d, xlab=expression(phi), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ lambda ------#
acf(PSI.chain$lambda,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$lambda, adjust = 5)
plot(d, xlab=expression(lambda), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)



rm(list=ls())

# GRAFICO - PREVISOES (diario) ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"STEMOS.Rda"))

#------ Graficos - Concordancia ------#
postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes4/",
                  "cap4_concord_prevnum_diarioD.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado", cex.lab=2.5,
     cex.axis=2.5)
points(as.numeric(Ft.numforecast[2,1:n,2:(T.new+1)]),
       as.numeric(Yt.k[,,2:(T.new+1)]), pch=19, cex=1.5)

dev.off()

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes4/",
                  "cap4_concord_prevnumcalib_diarioD.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado",cex.lab=2.5,
     cex.axis=2.5)
for (k in 2:(T.new+1))
{
  points(rowMedians(Yt.prev[1:n,,k]),
         as.numeric(Yt.k[,,k]),pch=19,cex=1.5)
}

dev.off()

#------ Graficos ------#
for (l in 1:59)
{
  print(l)
  
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes4/",
                    "cap4_previsoes_",as.character(coords.stations$codigo[l]),
                    "diarioD.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 1))
  plot(0,type="n",ylim=c(0,6),xlim=c(1,(T.new+2)),
       ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
       cex.lab=2.5)
  axis(1,1:6,c("-24","0","+24","+48","+72","+96"),cex.axis=2.5)
  axis(2,seq(0,6,2),seq(0,6,2),cex.axis=2.5)
  abline(h=0)
  polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
          c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
            colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
          col="grey80",border=NA)
  lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
  lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
  points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
  points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
  points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
  lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  
  if (l == 24) {
    legend("topleft", 
           legend = c("Eta", "Prev.", "Obs."), 
           horiz = TRUE, 
           lwd = c(2, 2, NA),
           lty = c(1, 2, NA),
           pch = c(NA, NA, 19),
           cex = 2.2)
  }
  
  dev.off()
}

rm(list=ls())

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"SEMOS.Rda"))

#------ Graficos - Concordancia ------#
postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes5/",
                  "cap4_concord_prevnum_diarioE.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado", cex.lab=2.5,
     cex.axis=2.5)
points(as.numeric(Ft.numforecast[2,1:n,2:(T.new+1)]),
       as.numeric(Yt.k[,,2:(T.new+1)]), pch=19, cex=1.5)

dev.off()

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes5/",
                  "cap4_concord_prevnumcalib_diarioE.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado",cex.lab=2.5,
     cex.axis=2.5)
for (k in 2:(T.new+1))
{
  points(rowMedians(Yt.prev[1:n,,k]),
         as.numeric(Yt.k[,,k]),pch=19,cex=1.5)
}

dev.off()

#------ Graficos ------#
for (l in 1:59)
{
  print(l)
  
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes5/",
                    "cap4_previsoes_",as.character(coords.stations$codigo[l]),
                    "diarioE.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 1))
  plot(0,type="n",ylim=c(0,6),xlim=c(1,(T.new+2)),
       ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
       cex.lab=2.5)
  axis(1,1:6,c("-24","0","+24","+48","+72","+96"),cex.axis=2.5)
  axis(2,seq(0,6,2),seq(0,6,2),cex.axis=2.5)
  abline(h=0)
  polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
          c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
            colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
          col="grey80",border=NA)
  lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
  lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
  points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
  points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
  points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
  lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  
  if (l == 24) {
    legend("topleft", 
           legend = c("Eta", "Prev.", "Obs."), 
           horiz = TRUE, 
           lwd = c(2, 2, NA),
           lty = c(1, 2, NA),
           pch = c(NA, NA, 19),
           cex = 2.2)
  }
  
  dev.off()
}

rm(list=ls())

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"DGOP.Rda"))

#------ Graficos - Concordancia ------#
postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes6/",
                  "cap4_concord_prevnum_diarioF.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado", cex.lab=2.5,
     cex.axis=2.5)
points(as.numeric(Ft.numforecast[2,1:n,2:(T.new+1)]),
       as.numeric(Yt.k[,,2:(T.new+1)]), pch=19, cex=1.5)

dev.off()

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes6/",
                  "cap4_concord_prevnumcalib_diarioF.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado",cex.lab=2.5,
     cex.axis=2.5)
for (k in 2:(T.new+1))
{
  points(rowMedians(Yt.prev[1:n,,k]),
         as.numeric(Yt.k[,,k]),pch=19,cex=1.5)
}

dev.off()

#------ Graficos ------#
for (l in 1:59)
{
  print(l)
  
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes6/",
                    "cap4_previsoes_",as.character(coords.stations$codigo[l]),
                    "diarioF.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 1))
  plot(0,type="n",ylim=c(0,6),xlim=c(1,(T.new+2)),
       ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
       cex.lab=2.5)
  axis(1,1:6,c("-24","0","+24","+48","+72","+96"),cex.axis=2.5)
  axis(2,seq(0,6,2),seq(0,6,2),cex.axis=2.5)
  abline(h=0)
  polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
          c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
            colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
          col="grey80",border=NA)
  lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
  lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
  points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
  points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
  points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
  lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  
  if (l == 24) {
    legend("topleft", 
           legend = c("Eta", "Prev.", "Obs."), 
           horiz = TRUE, 
           lwd = c(2, 2, NA),
           lty = c(1, 2, NA),
           pch = c(NA, NA, 19),
           cex = 2.2)
  }
  
  dev.off()
}

rm(list=ls())

# GRAFICO - CADEIAS (horario) ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
# load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
#             datas[data.c],"STEMOS2.Rda"))
# sampMH2 <- sampMH

load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"DGOP.Rda"))

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))

#------ phi ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,1]),
       sampMH$samples[30000:52500,1]),
     type="l",ylab=expression(phi),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,1],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)

#------ lambda ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[20000:60000,2]),
       sampMH$samples[30000:52500,2]),
     type="l",ylab=expression(lambda),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,2],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)

#------ Thetat ------#
par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))
for (m in (1:r))
{
  plot(1,type="n",
       ylim=c(min(graf.Thetat[3,m,-(1:30)],na.rm=T),max(graf.Thetat[1,m,-(1:30)],na.rm=T)),
       xlim=c(1,T),ylab=bquote(theta[.(m-1)]),xlab="Tempo",cex.lab=1.5,cex.axis=1.5)
  polygon(c(rev(1:T),1:T),c(rev(graf.Thetat[3,m,-1]),graf.Thetat[1,m,-1]),col="grey80",border=NA)
  lines(graf.Thetat[2,m,-1],lwd=2)
  abline(h=0,lty=3,lwd=2)
}

sigma.lambda.prior <-  10


#------ phi ------#
acf(PSI.chain$phi,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$phi, adjust = 5)
plot(d, xlab=expression(phi), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ lambda ------#
acf(PSI.chain$lambda,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$lambda, adjust = 5)
plot(d, xlab=expression(lambda), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)










rm(list=ls())



datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 11
# load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
#             datas[data.c],"STEMOS2.Rda"))
# sampMH2 <- sampMH

load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"STEMOS.Rda"))

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))

#------ b0 ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[50000:120000,1]),
       sampMH$samples[60000:140000,1]),
     type="l",ylab=expression(beta[0]),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,1],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,1],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ b1 ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[50000:120000,2]),
       sampMH$samples[60000:140000,2]),
     type="l",ylab=expression(beta[1]),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,2],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,2],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ phi ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[50000:120000,3]),
       sampMH$samples[60000:140000,3]),
     type="l",ylab=expression(phi),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,3],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,3],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ lambda ------#
plot(c(c(-.5,.2,1,1.5,.8),
       rev(sampMH$samples[50000:120000,4]),
       sampMH$samples[60000:140000,4]),
     type="l",ylab=expression(lambda),xlab="Iteracoes",cex.lab=1.5,cex.axis=1.5)
lines(sampMH$samples[,4],col="grey70")
abline(v=burnin*salto,lty=3,lwd=2)
acf(sampMH$samples[,4],lag.max = 25,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)

#------ Thetat ------#
par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))
for (m in (1:r))
{
  plot(1,type="n",
       ylim=c(min(graf.Thetat[3,m,-(1:30)],na.rm=T),max(graf.Thetat[1,m,-(1:30)],na.rm=T)),
       xlim=c(1,T),ylab=bquote(theta[.(m-1)]),xlab="Tempo",cex.lab=1.5,cex.axis=1.5)
  polygon(c(rev(1:T),1:T),c(rev(graf.Thetat[3,m,-1]),graf.Thetat[1,m,-1]),col="grey80",border=NA)
  lines(graf.Thetat[2,m,-1],lwd=2)
  abline(h=0,lty=3,lwd=2)
}

sigma.b0.prior <-  10
sigma.b1.prior <-  10
sigma.lambda.prior <-  10


#------ b0 ------#
par(mar=c(4.5,4.5,.5,.5))
acf(PSI.chain$b0,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$b0, adjust = 5)
plot(d, xlab=expression(beta[0]), main="", ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.b0.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ b1 ------#
acf(PSI.chain$b1,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$b1, adjust = 5,from = 0)
plot(d, xlab=expression(beta[1]), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.b1.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ phi ------#
acf(PSI.chain$phi,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$phi, adjust = 5)
plot(d, xlab=expression(phi), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

#------ lambda ------#
acf(PSI.chain$lambda,lag.max = 20,main="",ylab="FAC",xlab="Defasagem",cex.lab=1.5,cex.axis=1.5)
d <- density(x = PSI.chain$lambda, adjust = 5)
plot(d, xlab=expression(lambda), main="",ylab="Densidade",cex.lab=1.5,cex.axis=1.5)
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)

BETA0 <- data.frame(NULL)
BETA1 <- data.frame(NULL)
PHI <- data.frame(NULL)
LAMBDA <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
              datas[data.c],"STEMOS.Rda"))
  BETA0[1:2000,data.c] <- PSI.chain$b0
  BETA1[1:2000,data.c] <- PSI.chain$b1
  PHI[1:2000,data.c] <- PSI.chain$phi
  LAMBDA[1:2000,data.c] <- PSI.chain$lambda
}

names(BETA0) <- c("Dez","Jan","Fev","Mar","Abr","Mai",
                  "Jun","Jul","Ago","Set","Out","Nov")
names(BETA1) <- c("Dez","Jan","Fev","Mar","Abr","Mai",
                  "Jun","Jul","Ago","Set","Out","Nov")
names(PHI) <- c("Dez","Jan","Fev","Mar","Abr","Mai",
                "Jun","Jul","Ago","Set","Out","Nov")
names(LAMBDA) <- c("Dez","Jan","Fev","Mar","Abr","Mai",
                   "Jun","Jul","Ago","Set","Out","Nov")

par(mar=c(4.5,4.5,.5,.5))
par(mfrow=c(1,1))
A <- 
  boxplot(BETA0,outline = F,col="grey80",
          ylab=expression(beta[0]),xlab="Meses",
          cex.lab=1.5,cex.axis=1.5,plot = F)

A$stats[1,] <- colQuantiles(as.matrix(BETA0),probs = .025)
A$stats[2,] <- colQuantiles(as.matrix(BETA0),probs = .025)
A$stats[3,] <- colQuantiles(as.matrix(BETA0),probs = .5)
A$stats[4,] <- colQuantiles(as.matrix(BETA0),probs = .975)
A$stats[5,] <- colQuantiles(as.matrix(BETA0),probs = .975)
boxplot(A$stats,outline = F,names=c("Dez","Jan","Fev","Mar","Abr","Mai",
                                    "Jun","Jul","Ago","Set","Out","Nov"),
        ylab=expression(beta[0]),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,ylim=c(min(A$stats),1.05))

for (i in 1:12)
{
  segments(x0 = i,y0 = A$stats[2,i],x1 = i,y1 = A$stats[4,i],lty=2)
  segments(x0 = i+.4,y0 = A$stats[2,i],x1 = i+.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[2,i],x1 = i-.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[3,i],x1 = i+.4,y1 = A$stats[3,i],col="white",lwd=3)
  segments(x0 = i-.4,y0 = A$stats[1,i],x1 = i+.4,y1 = A$stats[1,i],lwd=3)
  segments(x0 = i-.4,y0 = A$stats[4,i],x1 = i+.4,y1 = A$stats[4,i],lwd=3)
  points(i,A$stats[3,i],pch=4,lwd=3,col=4)
}

abline(v=seq(0.5,12.5,3),col="grey80")
text(2,1.05,"Verão",cex=1.5)
text(5,1.05,"Outono",cex=1.5)
text(8,1.05,"Inverno",cex=1.5)
text(11,1.05,"Primavera",cex=1.5)

boxplot(BETA1,outline = F,col="grey80",
        ylab=expression(beta[1]),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5)

A <- 
  boxplot(BETA1,outline = F,col="grey80",
          ylab=expression(beta[0]),xlab="Meses",
          cex.lab=1.5,cex.axis=1.5,plot = F)

A$stats[1,] <- colQuantiles(as.matrix(BETA1),probs = .025)
A$stats[2,] <- colQuantiles(as.matrix(BETA1),probs = .025)
A$stats[3,] <- colQuantiles(as.matrix(BETA1),probs = .5)
A$stats[4,] <- colQuantiles(as.matrix(BETA1),probs = .975)
A$stats[5,] <- colQuantiles(as.matrix(BETA1),probs = .975)
boxplot(A$stats,outline = F,names=c("Dez","Jan","Fev","Mar","Abr","Mai",
                                    "Jun","Jul","Ago","Set","Out","Nov"),
        ylab=expression(beta[1]),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,ylim=c(min(A$stats),0.095))

for (i in 1:12)
{
  segments(x0 = i,y0 = A$stats[2,i],x1 = i,y1 = A$stats[4,i],lty=2)
  segments(x0 = i+.4,y0 = A$stats[2,i],x1 = i+.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[2,i],x1 = i-.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[3,i],x1 = i+.4,y1 = A$stats[3,i],col="white",lwd=3)
  segments(x0 = i-.4,y0 = A$stats[1,i],x1 = i+.4,y1 = A$stats[1,i],lwd=3)
  segments(x0 = i-.4,y0 = A$stats[4,i],x1 = i+.4,y1 = A$stats[4,i],lwd=3)
  points(i,A$stats[3,i],pch=4,lwd=3,col=4)
}

abline(v=seq(0.5,12.5,3),col="grey80")
text(2,0.095,"Verão",cex=1.5)
text(5,0.095,"Outono",cex=1.5)
text(8,0.095,"Inverno",cex=1.5)
text(11,0.095,"Primavera",cex=1.5)

boxplot(PHI,outline = F,col="grey80",
        ylab=expression(phi),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5)

A <- 
  boxplot(PHI,outline = F,col="grey80",
          ylab=expression(phi),xlab="Meses",
          cex.lab=1.5,cex.axis=1.5,plot = F)

A$stats[1,] <- colQuantiles(as.matrix(PHI),probs = .025)
A$stats[2,] <- colQuantiles(as.matrix(PHI),probs = .025)
A$stats[3,] <- colQuantiles(as.matrix(PHI),probs = .5)
A$stats[4,] <- colQuantiles(as.matrix(PHI),probs = .975)
A$stats[5,] <- colQuantiles(as.matrix(PHI),probs = .975)
boxplot(A$stats,outline = F,names=c("Dez","Jan","Fev","Mar","Abr","Mai",
                                    "Jun","Jul","Ago","Set","Out","Nov"),
        ylab=expression(phi),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,ylim=c(min(A$stats),11.5))

for (i in 1:12)
{
  segments(x0 = i,y0 = A$stats[2,i],x1 = i,y1 = A$stats[4,i],lty=2)
  segments(x0 = i+.4,y0 = A$stats[2,i],x1 = i+.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[2,i],x1 = i-.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[3,i],x1 = i+.4,y1 = A$stats[3,i],col="white",lwd=3)
  segments(x0 = i-.4,y0 = A$stats[1,i],x1 = i+.4,y1 = A$stats[1,i],lwd=3)
  segments(x0 = i-.4,y0 = A$stats[4,i],x1 = i+.4,y1 = A$stats[4,i],lwd=3)
  points(i,A$stats[3,i],pch=4,lwd=3,col=4)
}

abline(v=seq(0.5,12.5,3),col="grey80")
text(2,11.5,"Verão",cex=1.5)
text(5,11.5,"Outono",cex=1.5)
text(8,11.5,"Inverno",cex=1.5)
text(11,11.5,"Primavera",cex=1.5)

boxplot(LAMBDA,outline = F,col="grey80",
        ylab=expression(lambda),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5)

A <- 
  boxplot(LAMBDA,outline = F,col="grey80",
          ylab=expression(lambda),xlab="Meses",
          cex.lab=1.5,cex.axis=1.5,plot = F)

A$stats[1,] <- colQuantiles(as.matrix(LAMBDA),probs = .025)
A$stats[2,] <- colQuantiles(as.matrix(LAMBDA),probs = .025)
A$stats[3,] <- colQuantiles(as.matrix(LAMBDA),probs = .5)
A$stats[4,] <- colQuantiles(as.matrix(LAMBDA),probs = .975)
A$stats[5,] <- colQuantiles(as.matrix(LAMBDA),probs = .975)
boxplot(A$stats,outline = F,names=c("Dez","Jan","Fev","Mar","Abr","Mai",
                                    "Jun","Jul","Ago","Set","Out","Nov"),
        ylab=expression(lambda),xlab="Meses",
        cex.lab=1.5,cex.axis=1.5,ylim=c(min(A$stats),0.67))

for (i in 1:12)
{
  segments(x0 = i,y0 = A$stats[2,i],x1 = i,y1 = A$stats[4,i],lty=2)
  segments(x0 = i+.4,y0 = A$stats[2,i],x1 = i+.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[2,i],x1 = i-.4,y1 = A$stats[4,i],col="white")
  segments(x0 = i-.4,y0 = A$stats[3,i],x1 = i+.4,y1 = A$stats[3,i],col="white",lwd=3)
  segments(x0 = i-.4,y0 = A$stats[1,i],x1 = i+.4,y1 = A$stats[1,i],lwd=3)
  segments(x0 = i-.4,y0 = A$stats[4,i],x1 = i+.4,y1 = A$stats[4,i],lwd=3)
  points(i,A$stats[3,i],pch=4,lwd=3,col=4)
}

abline(v=seq(0.5,12.5,3),col="grey80")
text(2,0.67,"Verão",cex=1.5)
text(5,0.67,"Outono",cex=1.5)
text(8,0.67,"Inverno",cex=1.5)
text(11,0.67,"Primavera",cex=1.5)

rm(list=ls())

# GRAFICO - PREVISOES (horario) ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 9
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"STEMOS.Rda"))

#------ Graficos - Concordancia ------#
postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes4/",
                  "cap4_concord_prevnum_horarioD.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado", cex.lab=2.5,
     cex.axis=2.5)
points(as.numeric(Ft.numforecast[2,1:n,2:(T.new+1)]),
       as.numeric(Yt.k[,,2:(T.new+1)]), pch=19, cex=1.5)

dev.off()

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes4/",
                  "cap4_concord_prevnumcalib_horarioD.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado",cex.lab=2.5,
     cex.axis=2.5)
for (k in 2:(T.new+1))
{
  points(rowMedians(Yt.prev[1:n,,k]),
         as.numeric(Yt.k[,,k]),pch=19,cex=1.5)
}

dev.off()

#------ Graficos ------#
for (l in 1:59)
{
  print(l)
  
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes4/",
                    "cap4_previsoes_",as.character(coords.stations$codigo[l]),
                    "horarioD.eps"),
             width = 8, height = 6)
  
  if (l != 42 & l != 59) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,9),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,9,3),seq(0,9,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
    
    if (l == 8) {
      legend("topleft", 
             legend = c("Eta", "Prev.", "Obs."), 
             horiz = TRUE, 
             lwd = c(2, 2, NA),
             lty = c(1, 2, NA),
             pch = c(NA, NA, 19),
             cex = 2.2)
    }
    
  }
  
  if (l == 42) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,12),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,12,3),seq(0,12,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  }
  
  if (l == 59) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,15),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,15,3),seq(0,15,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  }
  
  dev.off()
}

rm(list=ls())

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")
data.c = 9
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ Graficos - Concordancia ------#
postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes5/",
                  "cap4_concord_prevnum_horarioE.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado", cex.lab=2.5,
     cex.axis=2.5)
points(as.numeric(Ft.numforecast[2,1:n,2:(T.new+1)]),
       as.numeric(Yt.k[,,2:(T.new+1)]), pch=19, cex=1.5)

dev.off()

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes5/",
                  "cap4_concord_prevnumcalib_horarioE.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado",cex.lab=2.5,
     cex.axis=2.5)
for (k in 2:(T.new+1))
{
  points(rowMedians(Yt.prev[1:n,,k]),
         as.numeric(Yt.k[,,k]),pch=19,cex=1.5)
}

dev.off()

#------ Graficos ------#
for (l in 1:59)
{
  print(l)
  
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes5/",
                    "cap4_previsoes_",as.character(coords.stations$codigo[l]),
                    "horarioE.eps"),
             width = 8, height = 6)
  
  if (l != 42 & l != 59) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,9),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,9,3),seq(0,9,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
    
    if (l == 8) {
      legend("topleft", 
             legend = c("Eta", "Prev.", "Obs."), 
             horiz = TRUE, 
             lwd = c(2, 2, NA),
             lty = c(1, 2, NA),
             pch = c(NA, NA, 19),
             cex = 2.2)
    }
    
  }
  
  if (l == 42) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,12),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,12,3),seq(0,12,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  }
  
  if (l == 59) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,15),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,15,3),seq(0,15,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  }
  
  dev.off()
}

rm(list=ls())

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")
data.c = 9
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"DGOP.Rda"))

#------ Graficos - Concordancia ------#
postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes6/",
                  "cap4_concord_prevnum_horarioF.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado", cex.lab=2.5,
     cex.axis=2.5)
points(as.numeric(Ft.numforecast[2,1:n,2:(T.new+1)]),
       as.numeric(Yt.k[,,2:(T.new+1)]), pch=19, cex=1.5)

dev.off()

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes6/",
                  "cap4_concord_prevnumcalib_horarioF.eps"),
           width = 6, height = 6,
           horizontal = TRUE)

par(mar=c(5, 5, .5, 1))
plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
     0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
     xlab="Previsto", 
     ylab="Observado",cex.lab=2.5,
     cex.axis=2.5)
for (k in 2:(T.new+1))
{
  points(rowMedians(Yt.prev[1:n,,k]),
         as.numeric(Yt.k[,,k]),pch=19,cex=1.5)
}

dev.off()

#------ Graficos ------#
for (l in 1:59)
{
  print(l)
  
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Previsoes6/",
                    "cap4_previsoes_",as.character(coords.stations$codigo[l]),
                    "horarioF.eps"),
             width = 8, height = 6)
  
  if (l != 42 & l != 59) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,9),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,9,3),seq(0,9,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
    
    if (l == 8) {
      legend("topleft", 
             legend = c("Eta", "Prev.", "Obs."), 
             horiz = TRUE, 
             lwd = c(2, 2, NA),
             lty = c(1, 2, NA),
             pch = c(NA, NA, 19),
             cex = 2.2)
    }
    
  }
  
  if (l == 42) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,12),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,12,3),seq(0,12,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  }
  
  if (l == 59) {
    
    par(mar=c(5, 5, 3, 1))
    plot(0,type="n",ylim=c(0,15),xlim=c(1,(T.new+2)),
         ylab="Vel. do Vento (m/s)",xlab="Horizonte (h)", axes= F,
         cex.lab=2.5)
    axis(1,seq(2,26,6),c("0",paste0("+",seq(6,24,6))),cex.axis=2.5)
    axis(2,seq(0,15,3),seq(0,15,3),cex.axis=2.5)
    axis(3,seq(2,26,6),c(paste0(seq(12,18,6),"h"),paste0(seq(0,12,6),"h")),cex.axis=2.5)
    abline(h=0)
    polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
            c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
              colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
            col="grey80",border=NA)
    lines(1:3,c(Ft[2,l,T],Ft[2,l,T+1],Ft.numforecast[2,l,2]),lwd=2,lty=2) #,t="o",pch=17,cex=2.5)
    lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lwd=2,lty=2) # ,t="o",pch=17,cex=2.5)
    points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,cex=2.5)
    points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,cex=2.5)
    points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,cex=2.5)
    lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lwd=2) #,t="o",pch="x",cex=2.5)
  }
  
  dev.off()
}

rm(list=ls())


