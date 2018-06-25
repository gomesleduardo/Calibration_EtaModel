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
# GRAFICO E CSV REQM, EAM E WILLMOTT INDEX AO LONGO DOS MESES - DIARIO
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
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
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
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
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
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
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
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.8,3.6),
     xlim=c(1,48),
     ylab="rEQM",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(2,48,4),c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(.8,3.6,,5),seq(.8,3.6,,5),cex.axis=1.5)
abline(v=seq(0,48,4),col="grey80")
abline(h=seq(.8,3.6,,5),col="grey80")
abline(v=c(0,3*4,6*4,9*4,12*4))
text(1.5*4,3.6,"Verão",cex=1.5)
text(4.5*4,3.6,"Outono",cex=1.5)
text(7.5*4,3.6,"Inverno",cex=1.5)
text(10.5*4,3.6,"Primavera",cex=1.5)
# lines(CRITERIOSREQM[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSREQM[,2],lty=2,lwd=2)
# lines(CRITERIOSREQM[,3],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSREQM[,4],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSREQM[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
lines(CRITERIOSREQM[,1],lwd=2,t="l",pch=17,cex=1.5)
lines(CRITERIOSREQM[,3],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSREQM[,4],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSREQM[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)

#------ Grafico ao longo das estacoes ------#
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.8,2.6),
     xlim=c(1,48),
     ylab="rEQM",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(2,48,4),c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(.8,2.6,,4),seq(.8,2.6,,4),cex.axis=1.5)
abline(v=seq(0,48,4),col="grey80")
abline(h=seq(.8,2.6,,4),col="grey80")
abline(v=c(0,3*4,6*4,9*4,12*4))
text(1.5*4,2.6,"Verão",cex=1.5)
text(4.5*4,2.6,"Outono",cex=1.5)
text(7.5*4,2.6,"Inverno",cex=1.5)
text(10.5*4,2.6,"Primavera",cex=1.5)
# lines(CRITERIOSREQM[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSREQM[,2],lty=2,lwd=2)
# lines(CRITERIOSREQM[,3],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSREQM[,4],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSREQM[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
# lines(CRITERIOSREQM[,1],lwd=2,t="l",pch=17,cex=1.5)
lines(CRITERIOSREQM[,3],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSREQM[,4],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSREQM[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)

par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.5,3),
     xlim=c(1,48),
     ylab="EAM",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(2,48,4),c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(0.5,3,.5),seq(0.5,3,.5),cex.axis=1.5)
abline(v=seq(0,48,4),col="grey80")
abline(h=seq(0,3,.5),col="grey80")
abline(v=c(0,3*4,6*4,9*4,12*4))
text(1.5*4,3,"Verão",cex=1.5)
text(4.5*4,3,"Outono",cex=1.5)
text(7.5*4,3,"Inverno",cex=1.5)
text(10.5*4,3,"Primavera",cex=1.5)
# lines(CRITERIOSEAM[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSEAM[,2],lty=2,lwd=2)
# lines(CRITERIOSEAM[,3],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSEAM[,4],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSEAM[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
lines(CRITERIOSEAM[,1],lwd=2,t="l",pch=17,cex=1.5)
lines(CRITERIOSEAM[,3],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSEAM[,4],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSEAM[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)

par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(0.5,2),
     xlim=c(1,48),
     ylab="EAM",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(2,48,4),c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(0.5,2,.5),seq(0.5,2,.5),cex.axis=1.5)
abline(v=seq(0,48,4),col="grey80")
abline(h=seq(0,2,.5),col="grey80")
abline(v=c(0,3*4,6*4,9*4,12*4))
text(1.5*4,2,"Verão",cex=1.5)
text(4.5*4,2,"Outono",cex=1.5)
text(7.5*4,2,"Inverno",cex=1.5)
text(10.5*4,2,"Primavera",cex=1.5)
# lines(CRITERIOSEAM[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSEAM[,2],lty=2,lwd=2)
# lines(CRITERIOSEAM[,3],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSEAM[,4],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSEAM[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
# lines(CRITERIOSEAM[,1],lwd=2,t="l",pch=17,cex=1.5)
lines(CRITERIOSEAM[,3],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSEAM[,4],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSEAM[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)


par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(.3,.9),
     xlim=c(1,48),
     ylab="ICW",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(2,48,4),c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(.3,.9,.15),seq(.3,.9,.15),cex.axis=1.5)
abline(v=seq(0,48,4),col="grey80")
abline(h=seq(.3,.9,.15),col="grey80")
abline(v=c(0,3*4,6*4,9*4,12*4))
text(1.5*4,.9,"Verão",cex=1.5)
text(4.5*4,.9,"Outono",cex=1.5)
text(7.5*4,.9,"Inverno",cex=1.5)
text(10.5*4,.9,"Primavera",cex=1.5)
# lines(CRITERIOSWillind[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSWillind[,2],lty=2,lwd=2)
# lines(CRITERIOSWillind[,3],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSWillind[,4],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSWillind[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
lines(CRITERIOSWillind[,1],lwd=2,t="l",pch=17,cex=1.5)
lines(CRITERIOSWillind[,3],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSWillind[,4],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSWillind[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)

par(mfrow=c(1,1))
par(mar=c(4.5,4.5,.5,.5))
plot(0,type="n",ylim=c(.3,.9),
     xlim=c(1,48),
     ylab="ICW",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(2,48,4),c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(.3,.9,.15),seq(.3,.9,.15),cex.axis=1.5)
abline(v=seq(0,48,4),col="grey80")
abline(h=seq(.3,.9,.15),col="grey80")
abline(v=c(0,3*4,6*4,9*4,12*4))
text(1.5*4,.9,"Verão",cex=1.5)
text(4.5*4,.9,"Outono",cex=1.5)
text(7.5*4,.9,"Inverno",cex=1.5)
text(10.5*4,.9,"Primavera",cex=1.5)
# lines(CRITERIOSWillind[,1],lwd=2,t="o",pch=17,cex=1.5)
# # lines(CRITERIOSWillind[,2],lty=2,lwd=2)
# lines(CRITERIOSWillind[,3],col="#7990c0",lwd=2,t="o",pch=15,cex=1.5)
# lines(CRITERIOSWillind[,4],col="#0000ff",lwd=2,t="o",pch=19,cex=1.5)
# lines(CRITERIOSWillind[,8],col="#00CC99",lwd=2,t="o",pch=18,cex=1.5)
# lines(CRITERIOSWillind[,1],lwd=2,t="l",pch=17,cex=1.5)
lines(CRITERIOSWillind[,3],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)
lines(CRITERIOSWillind[,4],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(CRITERIOSWillind[,8],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)


CRITERIOS <- data.frame(CRITERIOSREQM,CRITERIOSEAM,CRITERIOSWillind)
write.csv2(CRITERIOS,paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Graficos_Resultados/",
                            "Criterios_diario.csv"),row.names = F)
rm(list = ls())

#----------------------------------------------------------------------------------------------
# GRAFICO E CSV INTERVAL SCORE E COMPRIMENTO MEDIO DO IC AO LONGO DOS MESES - DIARIO
#----------------------------------------------------------------------------------------------

ICERROR <- data.frame(NULL)

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

IS <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
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
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
              datas[data.c],"STEMOS.Rda"))
  
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

ICERROR <- data.frame(ICERROR,IS)

IS <- data.frame(NULL)

for (data.c in 1:12)
{
  print(data.c)
  load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
              datas[data.c],"DGOP.Rda"))
  
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
     xlim=c(1,48),
     ylab="IS",xlab="Meses", axes= F,cex.lab=1.5)
axis(1,seq(2,48,4),c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),cex.axis=1.5)
axis(2,seq(floor(min(ICERROR,na.rm = T)),ceiling(max(ICERROR,na.rm = T)),2),
     seq(floor(min(ICERROR,na.rm = T)),ceiling(max(ICERROR,na.rm = T)),2),cex.axis=1.5)
abline(h=seq(3,15,2),col="grey80")
abline(v=seq(0,48,4),col="grey80")
abline(v=c(0,3*4,6*4,9*4,12*4))
text(1.5*4,15,"Verão",cex=1.5)
text(4.5*4,15,"Outono",cex=1.5)
text(7.5*4,15,"Inverno",cex=1.5)
text(10.5*4,15,"Primavera",cex=1.5)
lines(ICERROR[,3],col="#00CC99",lwd=2,t="l",pch=18,cex=1.5)
lines(ICERROR[,2],col="#0000ff",lwd=2,t="l",pch=19,cex=1.5)
lines(ICERROR[,1],col="#7990c0",lwd=2,t="l",pch=15,cex=1.5)