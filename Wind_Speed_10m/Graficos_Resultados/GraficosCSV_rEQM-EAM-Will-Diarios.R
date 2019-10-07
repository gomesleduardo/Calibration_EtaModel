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
           "geoR","adaptMCMC","tmvtnorm")

for (i in 1:length(pacotes))
{
  if (length(names(installed.packages()[,1])[names(installed.packages()[,1])==pacotes[i]])==0)
  {install.packages(pacotes[i], repos="http://cran.fiocruz.br/")}
  library(pacotes[i],character.only = TRUE)
}
rm(i,pacotes)

# GRAFICO E CSV REQM, EAM E WILLMOTT INDEX AO LONGO DOS MESES - DIARIO ----

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

#------ rEQM
REQM <- data.frame(reqm = c(CRITERIOSREQM$CRITERIOS.EnsREQM,
                            CRITERIOSREQM$CRITERIOS.mediaREQMSEMOS,
                            CRITERIOSREQM$CRITERIOS.mediaREQMSTEMOS,
                            CRITERIOSREQM$CRITERIOS.mediaREQMDGOP),
                   model = c(rep("Eta", 48), rep("SEMOS", 48),
                             rep("STEMOS", 48), rep("DGOP", 48)),
                   mes = rep(c(rep("Dez", 4),
                               rep("Jan", 4),
                               rep("Fev", 4),
                               rep("Mar", 4),
                               rep("Abr", 4),
                               rep("Mai", 4),
                               rep("Jun", 4),
                               rep("Jul", 4),
                               rep("Ago", 4),
                               rep("Set", 4),
                               rep("Out", 4),
                               rep("Nov", 4)), 4))

REQM$model <- factor(REQM$model, 
                     c("Eta", "SEMOS", "STEMOS", "DGOP"),
                     c("Eta", "SEMOS", "STEMOS", "DGOP"),
                     ordered = TRUE)

REQM$mes <- factor(REQM$mes, 
                   c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),
                   c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),
                   ordered = TRUE)

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/",
                  "cap4_rEQM_diario.eps"),
           width = 16, height = 3,
           horizontal = TRUE)

par(mar=c(4.5, 4.5, 2.5, 0))
layout(matrix(c(1,1,1,1,1,1,2), byrow = TRUE, 1, 7))
boxplots.triple <- boxplot(reqm ~ model + mes, data = REQM,
                           # subset = (season_c == seasons_list[i]),
                           main = "", outline = FALSE,
                           col = c("white",
                                   "grey90",
                                   "grey60",
                                   "grey30"),
                           xaxt ='n', axes = FALSE,
                           # ylim = c(.8, 3.5),
                           ylab = "rEQM (m/s)",
                           xlab = "", cex.axis = 1.8, 
                           cex.lab = 1.8,
                           pch = 19)
axis(1, seq(-1.5, 51, 4), 
     c("","Dez","Jan","Fev","Mar","Abr","Mai",
       "Jun","Jul","Ago","Set","Out","Nov",""), 
      las = 2, 
     cex.axis = 1.8, cex.lab = 1.8)
axis(2, cex.axis = 1.8)
axis(3, seq(-5.5, 60, 12),
     c("","Verão", "Outono", "Inverno", "Primavera",""),
     cex.axis = 1.8)
abline(v = seq(.5, 50, 4), lty = 2, col = "grey60")
abline(v = seq(.5, 50, 12), lty = 2, col = "grey30")


par(mar=c(0, 0, 0, 0))
plot(0,type='n', axes = FALSE)

legend("center",
       legend = c("Eta", "SEMOS", "STEMOS", "DGOP"), 
       fill = c("white", "grey90", "grey60", "grey30"),
       cex = 1.5)

dev.off()

#------ EAM
EAM <- data.frame(eam = c(CRITERIOSEAM$CRITERIOS.EnsEAM,
                           CRITERIOSEAM$CRITERIOS.mediaEAMSEMOS,
                           CRITERIOSEAM$CRITERIOS.mediaEAMSTEMOS,
                           CRITERIOSEAM$CRITERIOS.mediaEAMDGOP),
                   model = c(rep("Eta", 48), rep("SEMOS", 48),
                             rep("STEMOS", 48), rep("DGOP", 48)),
                   mes = rep(c(rep("Dez", 4),
                               rep("Jan", 4),
                               rep("Fev", 4),
                               rep("Mar", 4),
                               rep("Abr", 4),
                               rep("Mai", 4),
                               rep("Jun", 4),
                               rep("Jul", 4),
                               rep("Ago", 4),
                               rep("Set", 4),
                               rep("Out", 4),
                               rep("Nov", 4)), 4))

EAM$model <- factor(EAM$model, 
                     c("Eta", "SEMOS", "STEMOS", "DGOP"),
                     c("Eta", "SEMOS", "STEMOS", "DGOP"),
                     ordered = TRUE)

EAM$mes <- factor(EAM$mes, 
                   c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),
                   c("Dez","Jan","Fev","Mar","Abr","Mai",
                     "Jun","Jul","Ago","Set","Out","Nov"),
                   ordered = TRUE)

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/",
                  "cap4_EAM_diario.eps"),
           width = 16, height = 3,
           horizontal = TRUE)

par(mar=c(4.5, 4.5, 2.5, 0))
layout(matrix(c(1,1,1,1,1,1,2), byrow = TRUE, 1, 7))
boxplots.triple <- boxplot(eam ~ model + mes, data = EAM,
                           # subset = (season_c == seasons_list[i]),
                           main = "", outline = FALSE,
                           col = c("white",
                                   "grey90",
                                   "grey60",
                                   "grey30"),
                           xaxt ='n', axes = FALSE,
                           # ylim = c(.8, 3.5),
                           ylab = "EAM (m/s)",
                           xlab = "", cex.axis = 1.8, 
                           cex.lab = 1.8,
                           pch = 19)
axis(1, seq(-1.5, 51, 4), 
     c("","Dez","Jan","Fev","Mar","Abr","Mai",
       "Jun","Jul","Ago","Set","Out","Nov",""), 
     las = 2, 
     cex.axis = 1.8, cex.lab = 1.8)
axis(2, cex.axis = 1.8)
axis(3, seq(-5.5, 60, 12),
     c("","Verão", "Outono", "Inverno", "Primavera",""),
     cex.axis = 1.8)
abline(v = seq(.5, 50, 4), lty = 2, col = "grey60")
abline(v = seq(.5, 50, 12), lty = 2, col = "grey30")


par(mar=c(0, 0, 0, 0))
plot(0,type='n', axes = FALSE)

legend("center",
       legend = c("Eta", "SEMOS", "STEMOS", "DGOP"), 
       fill = c("white", "grey90", "grey60", "grey30"),
       cex = 1.5)

dev.off()

#------ ICW
ICW <- data.frame(icw = c(CRITERIOSWillind$CRITERIOS.MediaEnsWillind,
                          CRITERIOSWillind$CRITERIOS.mediaWillindSEMOS,
                          CRITERIOSWillind$CRITERIOS.mediaWillindSTEMOS,
                          CRITERIOSWillind$CRITERIOS.mediaWillindDGOP),
                  model = c(rep("Eta", 48), rep("SEMOS", 48),
                            rep("STEMOS", 48), rep("DGOP", 48)),
                  mes = rep(c(rep("Dez", 4),
                              rep("Jan", 4),
                              rep("Fev", 4),
                              rep("Mar", 4),
                              rep("Abr", 4),
                              rep("Mai", 4),
                              rep("Jun", 4),
                              rep("Jul", 4),
                              rep("Ago", 4),
                              rep("Set", 4),
                              rep("Out", 4),
                              rep("Nov", 4)), 4))

ICW$model <- factor(ICW$model, 
                    c("Eta", "SEMOS", "STEMOS", "DGOP"),
                    c("Eta", "SEMOS", "STEMOS", "DGOP"),
                    ordered = TRUE)

ICW$mes <- factor(ICW$mes, 
                  c("Dez","Jan","Fev","Mar","Abr","Mai",
                    "Jun","Jul","Ago","Set","Out","Nov"),
                  c("Dez","Jan","Fev","Mar","Abr","Mai",
                    "Jun","Jul","Ago","Set","Out","Nov"),
                  ordered = TRUE)

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/",
                  "cap4_ICW_diario.eps"),
           width = 16, height = 3,
           horizontal = TRUE)

par(mar=c(4.5, 4.5, 2.5, 0))
layout(matrix(c(1,1,1,1,1,1,2), byrow = TRUE, 1, 7))
boxplots.triple <- boxplot(icw ~ model + mes, data = ICW,
                           # subset = (season_c == seasons_list[i]),
                           main = "", outline = FALSE,
                           col = c("white",
                                   "grey90",
                                   "grey60",
                                   "grey30"),
                           xaxt ='n', axes = FALSE,
                           # ylim = c(.8, 3.5),
                           ylab = "ICW",
                           xlab = "", cex.axis = 1.8, 
                           cex.lab = 1.8,
                           pch = 19)
axis(1, seq(-1.5, 51, 4), 
     c("","Dez","Jan","Fev","Mar","Abr","Mai",
       "Jun","Jul","Ago","Set","Out","Nov",""), 
     las = 2, 
     cex.axis = 1.8, cex.lab = 1.8)
axis(2, cex.axis = 1.8)
axis(3, seq(-5.5, 60, 12),
     c("","Verão", "Outono", "Inverno", "Primavera",""),
     cex.axis = 1.8)
abline(v = seq(.5, 50, 4), lty = 2, col = "grey60")
abline(v = seq(.5, 50, 12), lty = 2, col = "grey30")


par(mar=c(0, 0, 0, 0))
plot(0,type='n', axes = FALSE)

legend("center",
       legend = c("Eta", "SEMOS", "STEMOS", "DGOP"), 
       fill = c("white", "grey90", "grey60", "grey30"),
       cex = 1.5)

dev.off()

CRITERIOS <- data.frame(CRITERIOSREQM,CRITERIOSEAM,CRITERIOSWillind)
write.csv2(CRITERIOS,paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/",
                            "Criterios_diario.csv"),row.names = F)
rm(list = ls())

# GRAFICO E CSV INTERVAL SCORE E COMPRIMENTO MEDIO DO IC AO LONGO DOS MESES - DIARIO ----

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
names(ICERROR) <- c("SEMOS","STEMOS","DGOP")

#------ IS
IS <- data.frame(is = c(ICERROR$SEMOS,
                        ICERROR$STEMOS,
                        ICERROR$DGOP),
                  model = c(rep("SEMOS", 48),
                            rep("STEMOS", 48), 
                            rep("DGOP", 48)),
                  mes = rep(c(rep("Dez", 4),
                              rep("Jan", 4),
                              rep("Fev", 4),
                              rep("Mar", 4),
                              rep("Abr", 4),
                              rep("Mai", 4),
                              rep("Jun", 4),
                              rep("Jul", 4),
                              rep("Ago", 4),
                              rep("Set", 4),
                              rep("Out", 4),
                              rep("Nov", 4)), 3))

IS$model <- factor(IS$model, 
                    c("SEMOS", "STEMOS", "DGOP"),
                    c("SEMOS", "STEMOS", "DGOP"),
                    ordered = TRUE)

IS$mes <- factor(IS$mes, 
                  c("Dez","Jan","Fev","Mar","Abr","Mai",
                    "Jun","Jul","Ago","Set","Out","Nov"),
                  c("Dez","Jan","Fev","Mar","Abr","Mai",
                    "Jun","Jul","Ago","Set","Out","Nov"),
                  ordered = TRUE)

postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/",
                  "cap4_IS_diario.eps"),
           width = 16, height = 3,
           horizontal = TRUE)

par(mar=c(4.5, 4.5, 2.5, 0))
layout(matrix(c(1,1,1,1,1,1,2), byrow = TRUE, 1, 7))
boxplots.triple <- boxplot(is ~ model + mes, data = IS,
                           # subset = (season_c == seasons_list[i]),
                           main = "", outline = FALSE,
                           col = c("grey90",
                                   "grey60",
                                   "grey30"),
                           xaxt ='n', axes = FALSE,
                           # ylim = c(.8, 3.5),
                           ylab = "IS",
                           xlab = "", cex.axis = 1.8, 
                           cex.lab = 1.8,
                           pch = 19)
axis(1, seq(-1, 39, 3), 
     c("","Dez","Jan","Fev","Mar","Abr","Mai",
       "Jun","Jul","Ago","Set","Out","Nov",""), 
     las = 2, 
     cex.axis = 1.8, cex.lab = 1.8)
axis(2, cex.axis = 1.8)
axis(3, seq(-4, 48, 9),
     c("","Verão", "Outono", "Inverno", "Primavera",""),
     cex.axis = 1.8)
abline(v = seq(.5, 38, 3), lty = 2, col = "grey60")
abline(v = seq(.5, 50, 12), lty = 2, col = "grey30")


par(mar=c(0, 0, 0, 0))
plot(0,type='n', axes = FALSE)

legend("center",
       legend = c("SEMOS", "STEMOS", "DGOP"), 
       fill = c("grey90", "grey60", "grey30"),
       cex = 1.5)

dev.off()

rm(list = ls())

