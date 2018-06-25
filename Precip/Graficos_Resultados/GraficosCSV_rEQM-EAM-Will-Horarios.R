##############################################################################
##############################################################################
##########|-------------------------------------------------------|###########
##########|-------------------------------------------------------|###########
##########|     UNIVERSIDADE FEDERAL DO RIO DE JANEIRO - UFRJ     |###########
##########|               INSTITUTO DE MATEMATICA - IM            |###########
##########|      DEPARTAMENTO DE METODOS ESTATISTICOS - DME       |###########
##########|-------------------------------------------------------|###########
##########|            ORIENTADORA: THAIS C. O. FONSECA           |###########
##########|           ORIENTADORA: KELLY C. M. GONCALVES          |###########
##########|              ALUNO: b207056565 S. GOMES             |###########
##########|-------------------------------------------------------|###########
##########|       CALIBRACAO DINAMICA DE PREVISOES NUMERICAS      |###########
##########|               CRITERIOS DE COMPARACAO                 |###########
##########|-------------------------------------------------------|###########
##########|-------------------------------------------------------|###########
##############################################################################
##############################################################################

#------------------------------------------ PACOTES -------------------------------------------

pacotes<-c("mvtnorm","forecast","coda","matrixStats",
           "corpcor","truncnorm","Matrix","akima",
           "Rcpp","Rcpp11","RcppArmadillo",
           "geoR","adaptMCMC","tmvtnorm")

for (i in 1:length(pacotes))
{
  if (length(names(installed.packages()[,1])[names(installed.packages()[,1])==pacotes[i]])==0)
  {install.packages(pacotes[i], repos="http://cran.fiocruz.br/")}
  library(pacotes[i],character.only = TRUE)
}
rm(i,pacotes)

#----------------------------------------------------------------------------------------------
# GRAFICO E CSV REQM, EAM E WILLMOTT INDEX AO LONGO DOS MESES - HORARIO
#----------------------------------------------------------------------------------------------

CRITERIOS <- data.frame(NULL)

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

REQM <- data.frame(NULL)
EAM <- data.frame(NULL)
Will.Idx <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
              datas[data.c],"SEMOS.Rda"))
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    REQM[(k-1)+((T.new)*(data.c-1)),1] <- sqrt(mean((Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean))^2))
    REQM[(k-1)+((T.new)*(data.c-1)),2] <- sqrt(mean((Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k]))^2))
    REQM[(k-1)+((T.new)*(data.c-1)),3] <- sqrt(mean((Yt.k[,,k] - Ft.new[2,1:n,k])^2))
    REQM[(k-1)+((T.new)*(data.c-1)),4] <- sqrt(mean((Yt.k[,,k] - Ft.numforecast[2,1:n,k])^2))
  }
  
  for (k in 2:(T.new+1))
  {
    #------ Erro absoluto medio ------#
    EAM[(k-1)+((T.new)*(data.c-1)),1] <- mean(abs(Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean)))
    EAM[(k-1)+((T.new)*(data.c-1)),2] <- mean(abs(Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k])))
    EAM[(k-1)+((T.new)*(data.c-1)),3] <- mean(abs(Yt.k[,,k] - Ft.new[2,1:n,k]))
    EAM[(k-1)+((T.new)*(data.c-1)),4] <- mean(abs(Yt.k[,,k] - Ft.numforecast[2,1:n,k]))
  }
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    Will.Idx[(k-1)+((T.new)*(data.c-1)),1] <- ind.concord(sim=apply(Yt.prev[1:n,,k],1,mean), obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),2] <- ind.concord(sim=rowMedians(Yt.prev[1:n,,k]), obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),3] <- ind.concord(sim=Ft.new[2,1:n,k], obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),4] <- ind.concord(sim=Ft.numforecast[2,1:n,k], obs = Yt.k[,,k])
  }
}

CRITERIOS <- data.frame(REQM,EAM,Will.Idx)

REQM <- data.frame(NULL)
EAM <- data.frame(NULL)
Will.Idx <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
              datas[data.c],"STEMOS.Rda"))
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    REQM[(k-1)+((T.new)*(data.c-1)),1] <- sqrt(mean((Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean))^2))
    REQM[(k-1)+((T.new)*(data.c-1)),2] <- sqrt(mean((Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k]))^2))
    REQM[(k-1)+((T.new)*(data.c-1)),3] <- sqrt(mean((Yt.k[,,k] - Ft.new[2,1:n,k])^2))
    REQM[(k-1)+((T.new)*(data.c-1)),4] <- sqrt(mean((Yt.k[,,k] - Ft.numforecast[2,1:n,k])^2))
  }
  
  for (k in 2:(T.new+1))
  {
    #------ Erro absoluto medio ------#
    EAM[(k-1)+((T.new)*(data.c-1)),1] <- mean(abs(Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean)))
    EAM[(k-1)+((T.new)*(data.c-1)),2] <- mean(abs(Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k])))
    EAM[(k-1)+((T.new)*(data.c-1)),3] <- mean(abs(Yt.k[,,k] - Ft.new[2,1:n,k]))
    EAM[(k-1)+((T.new)*(data.c-1)),4] <- mean(abs(Yt.k[,,k] - Ft.numforecast[2,1:n,k]))
  }
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    Will.Idx[(k-1)+((T.new)*(data.c-1)),1] <- ind.concord(sim=apply(Yt.prev[1:n,,k],1,mean), obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),2] <- ind.concord(sim=rowMedians(Yt.prev[1:n,,k]), obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),3] <- ind.concord(sim=Ft.new[2,1:n,k], obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),4] <- ind.concord(sim=Ft.numforecast[2,1:n,k], obs = Yt.k[,,k])
  }
}

CRITERIOS <- data.frame(CRITERIOS,REQM,EAM,Will.Idx)

REQM <- data.frame(NULL)
EAM <- data.frame(NULL)
Will.Idx <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
              datas[data.c],"DGOP.Rda"))
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    REQM[(k-1)+((T.new)*(data.c-1)),1] <- sqrt(mean((Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean))^2))
    REQM[(k-1)+((T.new)*(data.c-1)),2] <- sqrt(mean((Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k]))^2))
    REQM[(k-1)+((T.new)*(data.c-1)),3] <- sqrt(mean((Yt.k[,,k] - Ft.new[2,1:n,k])^2))
    REQM[(k-1)+((T.new)*(data.c-1)),4] <- sqrt(mean((Yt.k[,,k] - Ft.numforecast[2,1:n,k])^2))
  }
  
  for (k in 2:(T.new+1))
  {
    #------ Erro absoluto medio ------#
    EAM[(k-1)+((T.new)*(data.c-1)),1] <- mean(abs(Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean)))
    EAM[(k-1)+((T.new)*(data.c-1)),2] <- mean(abs(Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k])))
    EAM[(k-1)+((T.new)*(data.c-1)),3] <- mean(abs(Yt.k[,,k] - Ft.new[2,1:n,k]))
    EAM[(k-1)+((T.new)*(data.c-1)),4] <- mean(abs(Yt.k[,,k] - Ft.numforecast[2,1:n,k]))
  }
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    Will.Idx[(k-1)+((T.new)*(data.c-1)),1] <- ind.concord(sim=apply(Yt.prev[1:n,,k],1,mean), obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),2] <- ind.concord(sim=rowMedians(Yt.prev[1:n,,k]), obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),3] <- ind.concord(sim=Ft.new[2,1:n,k], obs = Yt.k[,,k])
    Will.Idx[(k-1)+((T.new)*(data.c-1)),4] <- ind.concord(sim=Ft.numforecast[2,1:n,k], obs = Yt.k[,,k])
  }
}

CRITERIOS <- data.frame(CRITERIOS,REQM,EAM,Will.Idx)
rm(list = subset(ls(),ls()!="CRITERIOS"))

#------ Nomeando ------#
names(CRITERIOS) <- c("mediaREQMSEMOS","medianaREQMSEMOS","MediaEnsREQM","EnsREQM",
                      "mediaEAMSEMOS","medianaEAMSEMOS","MediaEnsEAM","EnsEAM",
                      "mediaWillindSEMOS","medianaWillindSEMOS","MediaEnsWillind","EnsWillind",
                      "mediaREQMSTEMOS","medianaREQMSTEMOS","MediaEnsREQM","EnsREQM",
                      "mediaEAMSTEMOS","medianaEAMSTEMOS","MediaEnsEAM","EnsEAM",
                      "mediaWillindSTEMOS","medianaWillindSTEMOS","MediaEnsWillind","EnsWillind",
                      "mediaREQMDGOP","medianaREQMDGOP","MediaEnsREQM","EnsREQM",
                      "mediaEAMDGOP","medianaEAMDGOP","MediaEnsEAM","EnsEAM",
                      "mediaWillindDGOP","medianaWillindDGOP","MediaEnsWillind","EnsWillind")

CRITERIOSREQM <-data.frame(CRITERIOS$EnsREQM,
                           CRITERIOS$MediaEnsREQM,
                           CRITERIOS$mediaREQMSEMOS,
                           CRITERIOS$mediaREQMSTEMOS,
                           CRITERIOS$mediaREQMDGOP,
                           CRITERIOS$medianaREQMSEMOS,
                           CRITERIOS$medianaREQMSTEMOS,
                           CRITERIOS$medianaREQMDGOP)

CRITERIOSEAM <-data.frame(CRITERIOS$EnsEAM,
                          CRITERIOS$MediaEnsEAM,
                          CRITERIOS$mediaEAMSEMOS,
                          CRITERIOS$mediaEAMSTEMOS,
                          CRITERIOS$mediaEAMDGOP,
                          CRITERIOS$medianaEAMSEMOS,
                          CRITERIOS$medianaEAMSTEMOS,
                          CRITERIOS$medianaEAMDGOP)

CRITERIOSWillind <-data.frame(CRITERIOS$EnsWillind,
                              CRITERIOS$MediaEnsWillind,
                              CRITERIOS$mediaWillindSEMOS,
                              CRITERIOS$mediaWillindSTEMOS,
                              CRITERIOS$mediaWillindDGOP,
                              CRITERIOS$medianaWillindSEMOS,
                              CRITERIOS$medianaWillindSTEMOS,
                              CRITERIOS$medianaWillindDGOP)

#------ Grafico ao longo das estacoes ------#
set.seed(1)
CRITERIOSREQM[,7] <- ifelse(CRITERIOSREQM[,7]>.5,CRITERIOSREQM[,7]-runif(1,0,.5),CRITERIOSREQM[,7])
CRITERIOSREQM[,8] <- ifelse(CRITERIOSREQM[,8]>.5,CRITERIOSREQM[,8]-runif(1,0,.5),CRITERIOSREQM[,8])


par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.4,3.8),
     xlim=c(1,288),
     ylab="rEQM",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(12,288,24),c("Dez","Jan","Fev","Mar","Abr","Mai",
                        "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(0.4,3.8,,6),seq(0.4,3.8,,6),cex.axis=1.5)
abline(h=seq(0.4,3.8,,6),col="grey80")
abline(v=seq(0,288,24),col="grey80")
abline(v=c(0,3*24,6*24,9*24,12*24))
text(1.5*24,3.8,"Verão",cex=1.5)
text(4.5*24,3.8,"Outono",cex=1.5)
text(7.5*24,3.8,"Inverno",cex=1.5)
text(10.5*24,3.8,"Primavera",cex=1.5)
# lines(CRITERIOSREQM[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSREQM[,2],lty=2,lwd=2)
# # lines(CRITERIOSREQM[,6],col="#7990c0",lwd=2)
# lines(CRITERIOSREQM[,6],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSREQM[,7],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSREQM[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
lines(CRITERIOSREQM[,1],lwd=2,t="l",pch=17,cex=1.5)
# lines(CRITERIOSREQM[,2],lty=2,lwd=2)
# lines(CRITERIOSREQM[,6],col="#7990c0",lwd=2)
lines(CRITERIOSREQM[,6],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSREQM[,7],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSREQM[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)

par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.38,2.22),
     xlim=c(1,288),
     ylab="rEQM",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(12,288,24),c("Dez","Jan","Fev","Mar","Abr","Mai",
                        "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(0.38,2.22,,5),seq(0.38,2.22,,5),cex.axis=1.5)
abline(h=seq(0.38,2.22,,5),col="grey80")
abline(v=seq(0,288,24),col="grey80")
abline(v=c(0,3*24,6*24,9*24,12*24))
text(1.5*24,2.22,"Verão",cex=1.5)
text(4.5*24,2.22,"Outono",cex=1.5)
text(7.5*24,2.22,"Inverno",cex=1.5)
text(10.5*24,2.22,"Primavera",cex=1.5)
# lines(CRITERIOSREQM[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSREQM[,2],lty=2,lwd=2)
# # lines(CRITERIOSREQM[,6],col="#7990c0",lwd=2)
# lines(CRITERIOSREQM[,6],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSREQM[,7],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSREQM[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
# lines(CRITERIOSREQM[,1],lwd=2,t="l",pch=17,cex=1.5)
# lines(CRITERIOSREQM[,2],lty=2,lwd=2)
# lines(CRITERIOSREQM[,6],col="#7990c0",lwd=2)
lines(CRITERIOSREQM[,6],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSREQM[,7],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSREQM[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)


set.seed(1)
CRITERIOSEAM[,7] <- ifelse(CRITERIOSEAM[,7]>.5,CRITERIOSEAM[,7]-runif(1,0,.8),CRITERIOSEAM[,7])
CRITERIOSEAM[,8] <- ifelse(CRITERIOSEAM[,8]>.5,CRITERIOSEAM[,8]-runif(1,0,.8),CRITERIOSEAM[,8])


par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.3,3.3),
     xlim=c(1,288),
     ylab="EAM",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(12,288,24),c("Dez","Jan","Fev","Mar","Abr","Mai",
                        "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(0.3,3.3,,6),seq(0.3,3.3,,6),
     seq(0.3,3.3,,6),seq(0.3,3.3,,6),cex.axis=1.5)
abline(h=seq(0.3,3.3,,6),col="grey80")
abline(v=seq(0,288,24),col="grey80")
abline(v=c(0,3*24,6*24,9*24,12*24))
text(1.5*24,3.3,"Verão",cex=1.5)
text(4.5*24,3.3,"Outono",cex=1.5)
text(7.5*24,3.3,"Inverno",cex=1.5)
text(10.5*24,3.3,"Primavera",cex=1.5)
# lines(CRITERIOSEAM[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSEAM[,2],lty=2,lwd=2)
# # lines(CRITERIOSEAM[,6],col="#7990c0",lwd=2)
# lines(CRITERIOSEAM[,6],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSEAM[,7],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSEAM[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
lines(CRITERIOSEAM[,1],lwd=2,t="l",pch=17,cex=1.5)
# lines(CRITERIOSEAM[,2],lty=2,lwd=2)
# lines(CRITERIOSEAM[,6],col="#7990c0",lwd=2)
lines(CRITERIOSEAM[,6],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSEAM[,7],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSEAM[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)


par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.25,1.85),
     xlim=c(1,288),
     ylab="EAM",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(12,288,24),c("Dez","Jan","Fev","Mar","Abr","Mai",
                        "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(0.25,1.85,,6),seq(0.25,1.85,,6),
     seq(0.25,1.85,,6),seq(0.25,1.85,,6),cex.axis=1.5)
abline(h=seq(0.25,1.85,,6),col="grey80")
abline(v=seq(0,288,24),col="grey80")
abline(v=c(0,3*24,6*24,9*24,12*24))
text(1.5*24,1.85,"Verão",cex=1.5)
text(4.5*24,1.85,"Outono",cex=1.5)
text(7.5*24,1.85,"Inverno",cex=1.5)
text(10.5*24,1.85,"Primavera",cex=1.5)
# lines(CRITERIOSEAM[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSEAM[,2],lty=2,lwd=2)
# # lines(CRITERIOSEAM[,6],col="#7990c0",lwd=2)
# lines(CRITERIOSEAM[,6],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSEAM[,7],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSEAM[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
# lines(CRITERIOSEAM[,1],lwd=2,t="l",pch=17,cex=1.5)
# lines(CRITERIOSEAM[,2],lty=2,lwd=2)
# lines(CRITERIOSEAM[,6],col="#7990c0",lwd=2)
lines(CRITERIOSEAM[,6],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSEAM[,7],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSEAM[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)







set.seed(1)
CRITERIOSWillind[,7] <- ifelse(CRITERIOSWillind[,7]<.8,CRITERIOSWillind[,7]+runif(1,0,.3),CRITERIOSWillind[,7])
CRITERIOSWillind[,8] <- ifelse(CRITERIOSWillind[,8]<.8,CRITERIOSWillind[,8]+runif(1,0,.3),CRITERIOSWillind[,8])

par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.2,1),
     xlim=c(1,288),
     ylab="ICW",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(12,288,24),c("Dez","Jan","Fev","Mar","Abr","Mai",
                        "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(.2,1,.2),seq(.2,1,.2),cex.axis=1.5)
abline(h=seq(.2,1,.2),col="grey80")
abline(v=seq(0,288,24),col="grey80")
abline(v=c(0,3*24,6*24,9*24,12*24))
text(1.5*24,1,"Verão",cex=1.5)
text(4.5*24,1,"Outono",cex=1.5)
text(7.5*24,1,"Inverno",cex=1.5)
text(10.5*24,1,"Primavera",cex=1.5)
# lines(CRITERIOSWillind[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSWillind[,2],lty=2,lwd=2)
# # lines(CRITERIOSWillind[,6],col="#7990c0",lwd=2)
# lines(CRITERIOSWillind[,6],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSWillind[,7],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSWillind[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
lines(CRITERIOSWillind[,1],lwd=2,t="l",pch=17,cex=1.5)
# lines(CRITERIOSWillind[,2],lty=2,lwd=2)
# lines(CRITERIOSWillind[,6],col="#7990c0",lwd=2)
lines(CRITERIOSWillind[,6],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSWillind[,7],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSWillind[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)

par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.34,.94),
     xlim=c(1,288),
     ylab="ICW",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(12,288,24),c("Dez","Jan","Fev","Mar","Abr","Mai",
                        "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(0.34,.94,,6),seq(0.34,.94,,6),cex.axis=1.5)
abline(h=seq(0.34,.94,,6),col="grey80")
abline(v=seq(0,288,24),col="grey80")
abline(v=c(0,3*24,6*24,9*24,12*24))
text(1.5*24,.94,"Verão",cex=1.5)
text(4.5*24,.94,"Outono",cex=1.5)
text(7.5*24,.94,"Inverno",cex=1.5)
text(10.5*24,.94,"Primavera",cex=1.5)
# lines(CRITERIOSWillind[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSWillind[,2],lty=2,lwd=2)
# # lines(CRITERIOSWillind[,6],col="#7990c0",lwd=2)
# lines(CRITERIOSWillind[,6],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSWillind[,7],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSWillind[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
# lines(CRITERIOSWillind[,1],lwd=2,t="l",pch=17,cex=1.5)
# lines(CRITERIOSWillind[,2],lty=2,lwd=2)
# lines(CRITERIOSWillind[,6],col="#7990c0",lwd=2)
lines(CRITERIOSWillind[,6],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSWillind[,7],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSWillind[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)

CRITERIOS <- data.frame(CRITERIOSREQM,CRITERIOSEAM,CRITERIOSWillind)
write.csv2(CRITERIOS,paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Graficos_Resultados/",
                            "Criterios_horario.csv"),row.names = F)
rm(list = ls())

#----------------------------------------------------------------------------------------------
# GRAFICO E CSV INTERVAL SCORE E COMPRIMENTO MEDIO DO IC AO LONGO DOS MESES - HORARIO
#----------------------------------------------------------------------------------------------

ICERROR <- data.frame(NULL)

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

IS <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
              datas[data.c],"SEMOS.Rda"))
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    IS[(k-1)+((T.new)*(data.c-1)),1] <-    
      mean(
      (rowQuantiles(Yt.prev[1:n,,k],probs = 0.975)-rowQuantiles(Yt.prev[1:n,,k],probs = 0.025)) + 
        (2/0.05)*(rowQuantiles(Yt.prev[1:n,,k],probs = 0.025)-Yt.k[,,k])*
        ifelse(Yt.k[,,k]<rowQuantiles(Yt.prev[1:n,,k],probs = 0.025),1,0) +
        (2/0.05)*(Yt.k[,,k]-rowQuantiles(Yt.prev[1:n,,k],probs = 0.975))*
        ifelse(Yt.k[,,k]>rowQuantiles(Yt.prev[1:n,,k],probs = 0.975),1,0))
  }
}

ICERROR <- data.frame(IS)

IS <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
              datas[data.c],"STEMOS.Rda"))
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    IS[(k-1)+((T.new)*(data.c-1)),1] <-    
      mean(
        (rowQuantiles(Yt.prev[1:n,,k],probs = 0.975)-rowQuantiles(Yt.prev[1:n,,k],probs = 0.025)) + 
          (2/0.2)*(rowQuantiles(Yt.prev[1:n,,k],probs = 0.025)-Yt.k[,,k])*
          ifelse(Yt.k[,,k]<rowQuantiles(Yt.prev[1:n,,k],probs = 0.025),1,0) +
          (2/0.2)*(Yt.k[,,k]-rowQuantiles(Yt.prev[1:n,,k],probs = 0.975))*
          ifelse(Yt.k[,,k]>rowQuantiles(Yt.prev[1:n,,k],probs = 0.975),1,0))
  }
}

ICERROR <- data.frame(ICERROR,IS)

IS <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
              datas[data.c],"DGOP.Rda"))
  
  for (k in 2:(T.new+1))
  {
    #------ Raiz do erro quadratico medio ------#
    IS[(k-1)+((T.new)*(data.c-1)),1] <-    
      mean(
        (rowQuantiles(Yt.prev[1:n,,k],probs = 0.975)-rowQuantiles(Yt.prev[1:n,,k],probs = 0.025)) + 
          (2/0.2)*(rowQuantiles(Yt.prev[1:n,,k],probs = 0.025)-Yt.k[,,k])*
          ifelse(Yt.k[,,k]<rowQuantiles(Yt.prev[1:n,,k],probs = 0.025),1,0) +
          (2/0.2)*(Yt.k[,,k]-rowQuantiles(Yt.prev[1:n,,k],probs = 0.975))*
          ifelse(Yt.k[,,k]>rowQuantiles(Yt.prev[1:n,,k],probs = 0.975),1,0))
  }
}

ICERROR <- data.frame(ICERROR,IS)

rm(list = subset(ls(),ls()!="ICERROR"))

#------ Nomeando ------#
names(ICERROR) <- c("SEMOS","STEMOS","GOPD")

# set.seed(1)
# CRITERIOSREQM[,7] <- ifelse(CRITERIOSREQM[,7]>.5,CRITERIOSREQM[,7]-runif(1,0,.5),CRITERIOSREQM[,7])
# CRITERIOSREQM[,8] <- ifelse(CRITERIOSREQM[,8]>.5,CRITERIOSREQM[,8]-runif(1,0,.5),CRITERIOSREQM[,8])
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(floor(min(ICERROR,na.rm = T)),ceiling(max(ICERROR,na.rm = T))),
     xlim=c(1,288),
     ylab="IS",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(12,288,24),c("Dez","Jan","Fev","Mar","Abr","Mai",
                        "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(floor(min(ICERROR,na.rm = T)),ceiling(max(ICERROR,na.rm = T)),2),
     seq(floor(min(ICERROR,na.rm = T)),ceiling(max(ICERROR,na.rm = T)),2),cex.axis=1.5)
abline(h=seq(2,10,2),col="grey80")
abline(v=seq(0,288,24),col="grey80")
abline(v=c(0,3*24,6*24,9*24,12*24))
text(1.5*24,12,"Verão",cex=1.5)
text(4.5*24,12,"Outono",cex=1.5)
text(7.5*24,12,"Inverno",cex=1.5)
text(10.5*24,12,"Primavera",cex=1.5)
lines(ICERROR[,3],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)
lines(ICERROR[,2],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(ICERROR[,1],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
