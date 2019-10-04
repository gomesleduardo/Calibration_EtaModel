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
           "geoR","adaptMCMC","tmvtnorm", "fields", "raster",
           "sp","colorRamps","RColorBrewer")

for (i in 1:length(pacotes))
{
  if (length(names(installed.packages()[,1])[names(installed.packages()[,1])==pacotes[i]])==0)
  {install.packages(pacotes[i], repos="http://cran.fiocruz.br/")}
  library(pacotes[i],character.only = TRUE)
}
rm(i,pacotes)

# GRAFICO - FAC ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 2
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ FAC ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/FAC/",
                    "cap1_FAC_",as.character(coords.stations$codigo[i]),
                    "verao.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  plot(0,type="n",ylim=c(-.5,1),
       xlim=c(0,120),
       ylab="FAC",
       xlab="Defasagem (h)", 
       axes= F, cex.lab = 2.5)
  axis(1,seq(0,120,24),seq(0,120,24),cex.axis=2.5)
  axis(2,cex.axis=2.5)
  points(0:120,as.numeric(acf(Yt[i,,-1],lag.max = 120,plot=F)[[1]]),type="h",lwd=2)
  
  dev.off()
}

data.c = 5
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ FAC ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/FAC/",
                    "cap1_FAC_",as.character(coords.stations$codigo[i]),
                    "outono.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  plot(0,type="n",ylim=c(-.5,1),
       xlim=c(0,120),
       ylab="FAC",
       xlab="Defasagem (h)", 
       axes= F, cex.lab = 2.5)
  axis(1,seq(0,120,24),seq(0,120,24),cex.axis=2.5)
  axis(2,cex.axis=2.5)
  points(0:120,as.numeric(acf(Yt[i,,-1],lag.max = 120,plot=F)[[1]]),type="h",lwd=2)
  
  dev.off()
}

data.c = 8
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ FAC ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/FAC/",
                    "cap1_FAC_",as.character(coords.stations$codigo[i]),
                    "inverno.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  plot(0,type="n",ylim=c(-.5,1),
       xlim=c(0,120),
       ylab="FAC",
       xlab="Defasagem (h)", 
       axes= F, cex.lab = 2.5)
  axis(1,seq(0,120,24),seq(0,120,24),cex.axis=2.5)
  axis(2,cex.axis=2.5)
  points(0:120,as.numeric(acf(Yt[i,,-1],lag.max = 120,plot=F)[[1]]),type="h",lwd=2)
  
  dev.off()
}

data.c = 11
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ FAC ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/FAC/",
                    "cap1_FAC_",as.character(coords.stations$codigo[i]),
                    "primavera.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  plot(0,type="n",ylim=c(-.5,1),
       xlim=c(0,120),
       ylab="FAC",
       xlab="Defasagem (h)", 
       axes= F, cex.lab = 2.5)
  axis(1,seq(0,120,24),seq(0,120,24),cex.axis=2.5)
  axis(2,cex.axis=2.5)
  points(0:120,as.numeric(acf(Yt[i,,-1],lag.max = 120,plot=F)[[1]]),type="h",lwd=2)
  
  dev.off()
}

rm(list = ls())

# GRAFICO - HISTOGRAMA ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 2
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ Histograma ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Histograma/",
                    "cap1_histo_",as.character(coords.stations$codigo[i]),
                    "verao.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  hist(Yt[i,,-1],freq = F,main="",ylab="Densidade",
       xlab="Vel. do Vento (m/s)",col="grey80",ylim=c(0,.8),
       cex.axis=2.5,cex.lab=2.5,lwd=2,xlim=c(0,12))
  
  dev.off()
}

data.c = 5
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ Histograma ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Histograma/",
                    "cap1_histo_",as.character(coords.stations$codigo[i]),
                    "outono.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  hist(Yt[i,,-1],freq = F,main="",ylab="Densidade",
       xlab="Vel. do Vento (m/s)",col="grey80",ylim=c(0,.8),
       cex.axis=2.5,cex.lab=2.5,lwd=2,xlim=c(0,12))
  
  dev.off()
}

data.c = 8
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ Histograma ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Histograma/",
                    "cap1_histo_",as.character(coords.stations$codigo[i]),
                    "inverno.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  hist(Yt[i,,-1],freq = F,main="",ylab="Densidade",
       xlab="Vel. do Vento (m/s)",col="grey80",ylim=c(0,.8),
       cex.axis=2.5,cex.lab=2.5,lwd=2,xlim=c(0,12))
  
  dev.off()
}

data.c = 11
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ Histograma ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Histograma/",
                    "cap1_histo_",as.character(coords.stations$codigo[i]),
                    "primavera.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  hist(Yt[i,,-1],freq = F,main="",ylab="Densidade",
       xlab="Vel. do Vento (m/s)",col="grey80",ylim=c(0,.8),
       cex.axis=2.5,cex.lab=2.5,lwd=2,xlim=c(0,12))
  
  dev.off()
}

rm(list = ls())

# GRAFICO - ST (comparativos) ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 2
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ ST ------#
Yt <- ifelse(Yt<.1,0,Yt)
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/ST/",
                    "cap1_ST_",as.character(coords.stations$codigo[i]),
                    "verao.eps"),
             width = 8, height = 6)

  par(mar=c(5, 5, .5, 2))
  plot(0,type="n",ylim=c(0,9),
       xlim=c(0,120),
       ylab="Vel. do Vento (m/s)",
       xlab="Horizonte (h)", axes = F,cex.lab = 2.5)
  axis(1,seq(0,120,24),c(0,paste0("+",seq(24,120,24))),cex.axis=2.5)
  axis(2,seq(0,9,3),seq(0,9,3),cex.axis=3)
  lines(Yt[i,,2:121],t="l",lwd=2)
  lines(Ft[2,i,2:121],t="l",lwd=2, lty = 2)
  
  legend("topleft",
         legend = c("Obs.", 
                    "Eta"), 
         lwd = c(2, 2),
         lty = c(1, 2),
         horiz = TRUE,
         cex = 2.2)
  
  dev.off()
}

data.c = 5
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ ST ------#
Yt <- ifelse(Yt<.1,0,Yt)
par(mar=c(6.5, 5.5, .5, 1))
for (i in 1:n)
{
  
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/ST/",
                    "cap1_ST_",as.character(coords.stations$codigo[i]),
                    "outono.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  plot(0,type="n",ylim=c(0,9),
       xlim=c(0,120),
       ylab="Vel. do Vento (m/s)",
       xlab="Horizonte (h)", axes = F,cex.lab = 2.5)
  axis(1,seq(0,120,24),c(0,paste0("+",seq(24,120,24))),cex.axis=2.5)
  axis(2,seq(0,9,3),seq(0,9,3),cex.axis=3)
  lines(Yt[i,,2:121],t="l",lwd=2)
  lines(Ft[2,i,2:121],t="l",lwd=2, lty = 2)

  dev.off()
}

data.c = 8
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ ST ------#
Yt <- ifelse(Yt<.1,0,Yt)
par(mar=c(6.5, 5.5, .5, 1))
for (i in 1:n)
{
  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/ST/",
                    "cap1_ST_",as.character(coords.stations$codigo[i]),
                    "inverno.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  plot(0,type="n",ylim=c(0,9),
       xlim=c(0,120),
       ylab="Vel. do Vento (m/s)",
       xlab="Horizonte (h)", axes = F,cex.lab = 2.5)
  axis(1,seq(0,120,24),c(0,paste0("+",seq(24,120,24))),cex.axis=2.5)
  axis(2,seq(0,9,3),seq(0,9,3),cex.axis=3)
  lines(Yt[i,,2:121],t="l",lwd=2)
  lines(Ft[2,i,2:121],t="l",lwd=2, lty = 2)
  
  dev.off()
}

data.c = 11
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

#------ ST ------#
Yt <- ifelse(Yt<.1,0,Yt)
par(mar=c(6.5, 5.5, .5, 1))
for (i in 1:n)
{

  print(i)
  postscript(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/ST/",
                    "cap1_ST_",as.character(coords.stations$codigo[i]),
                    "primavera.eps"),
             width = 8, height = 6)
  
  par(mar=c(5, 5, .5, 2))
  plot(0,type="n",ylim=c(0,9),
       xlim=c(0,120),
       ylab="Vel. do Vento (m/s)",
       xlab="Horizonte (h)", axes = F,cex.lab = 2.5)
  axis(1,seq(0,120,24),c(0,paste0("+",seq(24,120,24))),cex.axis=2.5)
  axis(2,seq(0,9,3),seq(0,9,3),cex.axis=3)
  lines(Yt[i,,2:121],t="l",lwd=2)
  lines(Ft[2,i,2:121],t="l",lwd=2, lty = 2)
  
  dev.off()
}

rm(list = ls())

# CORRELACAO DO ENSEMBLE ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 8
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"SEMOS.Rda"))

AUX <- subset(wind10m.ETA,wind10m.ETA$data == datas[data.c])[,1:(n+3)]
AUX <- AUX[order(AUX$horizonte),]
S2.aux <- data.frame(t(data.frame(c("Varianca amostral","",NA,colVars(as.matrix(AUX[,-(1:3)]))),t(AUX))))

AUX2 <- data.frame(t(AUX[,-(1:3)]))
names(AUX2) <- paste0("-",seq(96,240,24))

plot(AUX2,pch=19,cex.lab=1.5,cex.axis=2.5)
rm(list = ls())

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 7
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"SEMOS.Rda"))

AUX <- subset(wind10m.ETA,wind10m.ETA$data == datas[data.c])[,1:(n+3)]
AUX <- AUX[order(AUX$horizonte),]
S2.aux <- data.frame(t(data.frame(c("Varianca amostral","",NA,colVars(as.matrix(AUX[,-(1:3)]))),t(AUX))))

AUX2 <- data.frame(t(AUX[,-(1:3)]))
names(AUX2) <- paste0("+",seq(96,240,24),"h")


M <- cor(AUX2)

corrplot::corrplot(M, method = "color", type = "lower", diag = FALSE,
                   addCoef.col = "black", 
                   col = c(rep("#E6E6E6rev", 10),rev(grey.colors(10))),
                   tl.col = "black", tl.srt = 0, tl.cex = 1.1, tl.offset = .8,
                   cl.cex = 1.1, cl.lim = c(.5,1))
rect(.5, .5, 6.5, 7.5)

rm(list = ls())

# BOX PLOT HORARIO ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 12
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Horarios/result_horario",
            datas[data.c],"SEMOS.Rda"))

AUX <- data.frame(NULL)
Yt <- ifelse(Yt<.1,0,Yt)
for (w in 0:23)
{
  AUX[1:59,w+1] <- as.numeric(colMeans(subset(obs.wind,
                        as.numeric(substr(obs.wind$data.hora,12,13))==w)[,-1]))
}
names(AUX) <- 0:23

par(mar=c(4.5, 5, .5, 2.1))
boxplot(AUX,outline = F,col="grey80",
        ylab="Velocidade do Vento (m/s)",xlab="Hora",
        cex.lab=1.5,cex.axis=2.5)
rm(list = ls())

# GRAFICO - INTERPOLACAO BILINEAR ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 2
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"SEMOS.Rda"))

AUX <- aux.csv[1,]
wind10m.ETA2 <- data.frame(t(suppressWarnings(interpp(x = coords.ETA$longitude,
                             y = coords.ETA$latitude,
                             z = as.numeric(AUX[,-(1:3)]),
                             xo = gradeMG$longitude,
                             yo = gradeMG$latitude)$z)))

#------ Ajustando DF ------#
df.graf <- data.frame(longitude = coords.ETA$longitude,
                      latitude = coords.ETA$latitude,
                      prevnum = as.numeric(AUX[,-(1:3)]))

#------ Criando SPDF ------#
df.graf2 <- df.graf
coordinates(df.graf2) = c("longitude", "latitude") # promote to SpatialPointsDataFrame
gridded(df.graf2) <- TRUE # promote to SpatialPixelsDataFrame

#------ Ajustando DF ------#
df.graf3 <- data.frame(longitude = gradeMG$longitude,
                       latitude = gradeMG$latitude,
                       previnterpbil = as.numeric(wind10m.ETA2))

#------ Criando SPDF ------#
df.graf4 <- df.graf3
coordinates(df.graf4) = c("longitude", "latitude") # promote to SpatialPointsDataFrame
gridded(df.graf4) <- TRUE # promote to SpatialPixelsDataFrame

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/',
                  'cap1_prevnumA.eps'),
           width = 8, height = 6)

#------ Exemplo de Previsao Numerica ------#
par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c(1:3)])
raster::plot(r, axes = FALSE,  
             zlim = c(0,ceiling(max(df.graf$prevnum))),
             breaks = seq(0, ceiling(max(df.graf$prevnum)), 1),
             col = brewer.pal(ceiling(max(df.graf$prevnum)),"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, ceiling(max(df.graf$prevnum)), 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/',
                  'cap1_prevnumB.eps'),
           width = 8, height = 6)

#------ Exemplo de Previsao Numerica ------#
par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf3[, c(1:3)])
raster::plot(r, axes = FALSE,  
             zlim = c(0,ceiling(max(df.graf3$previnterpbil))),
             breaks = seq(0, ceiling(max(df.graf3$previnterpbil)), 1),
             col = brewer.pal(ceiling(max(df.graf3$previnterpbil)),"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, ceiling(max(df.graf3$previnterpbil)), 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0))
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

#------ Localizacao das estacoes ------#
par(mar=c(4.5,4.5,.5,.5))
plot(geometry(df.graf2), col = grey(.8),ylab = "Latitude",xlab="Longitude",cex.lab=1.2)
lines(coords.MG)
# points(coords.stations[,2:3],pch=17,cex=1.2,col="#7990c0")
points(coords.stations[,2:3],pch=17,cex=1.2,col="grey50")
axis(1,floor(min(coords.ETA$longitude)):ceiling(max(coords.ETA$longitude)),cex.axis=1.2)
axis(2,floor(min(coords.ETA$latitude)):ceiling(max(coords.ETA$latitude)),cex.axis=1.2)
text(x = coords.stations[,2],y = coords.stations[,3],
     labels = coords.stations[,5],
     adj = 0, pos =1,
     cex=.8)

rm(list = ls())

# GRAFICO - LOCALIZACOES, MAPAS DE PREVISAO ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")
data.c = 8
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/mapa2",
            datas[data.c],"DGOP.Rda"))

#------ Ajustando DF ------#
df.graf <- data.frame(longitude = gradeMG$longitude,
                      latitude = gradeMG$latitude,
                      hora18 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                        y = coords.ETA$latitude,
                                                        z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .5),
                                                        xo = gradeMG$longitude,
                                                        yo = gradeMG$latitude)$z),
                      hora18IC1 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .1),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora18IC2 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .9),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora18ME = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                          y = coords.ETA$latitude,
                                                          z = (rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .9)-
                                                                 rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .1))/2,
                                                          xo = gradeMG$longitude,
                                                          yo = gradeMG$latitude)$z),
                      hora00 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                        y = coords.ETA$latitude,
                                                        z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .5),
                                                        xo = gradeMG$longitude,
                                                        yo = gradeMG$latitude)$z),
                      hora00IC1 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .1),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora00IC2 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .9),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora00ME = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                          y = coords.ETA$latitude,
                                                          z = (rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .9)-
                                                                 rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .1))/2,
                                                          xo = gradeMG$longitude,
                                                          yo = gradeMG$latitude)$z),
                      hora06 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                        y = coords.ETA$latitude,
                                                        z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .5),
                                                        xo = gradeMG$longitude,
                                                        yo = gradeMG$latitude)$z),
                      hora06IC1 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .1),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora06IC2 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .9),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora06ME = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                          y = coords.ETA$latitude,
                                                          z = (rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .9)-
                                                                 rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .1))/2,
                                                          xo = gradeMG$longitude,
                                                          yo = gradeMG$latitude)$z),
                      hora12 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                        y = coords.ETA$latitude,
                                                        z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .5),
                                                        xo = gradeMG$longitude,
                                                        yo = gradeMG$latitude)$z),
                      hora12IC1 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .1),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora12IC2 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .9),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora12ME = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                          y = coords.ETA$latitude,
                                                          z = (rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .9)-
                                                                 rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .1))/2,
                                                          xo = gradeMG$longitude,
                                                          yo = gradeMG$latitude)$z),
                      ETA18 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                       y = coords.ETA$latitude,
                                                       z = Ft.numforecast[2,-(1:59),7],
                                                       xo = gradeMG$longitude,
                                                       yo = gradeMG$latitude)$z),
                      ETA00 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                       y = coords.ETA$latitude,
                                                       z = Ft.numforecast[2,-(1:59),13],
                                                       xo = gradeMG$longitude,
                                                       yo = gradeMG$latitude)$z),
                      ETA06 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                       y = coords.ETA$latitude,
                                                       z = Ft.numforecast[2,-(1:59),19],
                                                       xo = gradeMG$longitude,
                                                       yo = gradeMG$latitude)$z),
                      ETA12 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                       y = coords.ETA$latitude,
                                                       z = Ft.numforecast[2,-(1:59),25],
                                                       xo = gradeMG$longitude,
                                                       yo = gradeMG$latitude)$z))

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_ETA18.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "ETA18")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_ETA00.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "ETA00")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_ETA06.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "ETA06")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_ETA12.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "ETA12")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_hora18.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "hora18")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_hora00.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "hora00")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_hora06.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "hora06")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_hora12.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "hora12")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_hora18ME.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "hora18ME")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_hora00ME.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "hora00ME")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_hora06ME.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "hora06ME")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

postscript(paste0('C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Graficos_Resultados/Mapas/',
                  'cap4_mapas_hora12ME.eps'),
           width = 8, height = 6)

par(mar=c(4.5, 4.5, .5, .5))
r <- raster::rasterFromXYZ(df.graf[, c("longitude", "latitude", "hora12ME")])
raster::plot(r, axes = FALSE,  
             zlim = c(0, 9),
             breaks = seq(0, 9, 1),
             col = brewer.pal(9,"Greys"),
             xlab = "Longitude", ylab = "Latitude", cex.lab = 1.8,
             legend.width = 1.5, legend.shrink = 1, 
             axis.args = list(at = seq(0, 9, 2), 
                              cex.axis = 1.8))
box(col = 'white', lwd = 3)
# plot(geometry(df.graf2), add = TRUE, col = grey(0), lwd = .01, lty = 2)
lines(coords.MG)
axis(1, seq(-50,-40,2), cex.axis = 1.8)
axis(2,seq(-22,-16,2), cex.axis = 1.8)

dev.off()

rm(list = ls())

# GRAFICO - LOCALIZACOES GOOGLE MAPS ----

library(ggmap,RgoogleMaps,scales)

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

data.c = 2
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
            datas[data.c],"SEMOS.Rda"))

AUX <- aux.csv[1,]
wind10m.ETA2 <- data.frame(t(suppressWarnings(interpp(x = coords.ETA$longitude,
                                                      y = coords.ETA$latitude,
                                                      z = as.numeric(AUX[,-(1:3)]),
                                                      xo = gradeMG$longitude,
                                                      yo = gradeMG$latitude)$z)))

#------ Ajustando DF ------#
df.graf <- data.frame(longitude = coords.ETA$longitude,
                      latitude = coords.ETA$latitude,
                      prevnum = as.numeric(AUX[,-(1:3)]))

#------ Criando SPDF ------#
df.graf2 <- df.graf
coordinates(df.graf2) = c("longitude", "latitude") # promote to SpatialPointsDataFrame
gridded(df.graf2) <- TRUE # promote to SpatialPixelsDataFrame

map <-
  get_map(location = c(lon = -44.80567, lat = -18.4764), zoom = 6,
          scale = "auto", maptype = c("terrain"), source = c("google"), 
          force = ifelse(source == "google", TRUE, TRUE),
          messaging = FALSE, urlonly = FALSE, filename = "ggmapTemp",
          crop = TRUE, color = c("color"), language = "pt-BR")



mapPoints <- ggmap(map)

# mapPoints
# 
# mapPoints +
#   geom_point(aes(x = longitude,
#                  y = latitude), 
#              col = "orange",data = coords.stations)

mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "black",pch=22,cex=2.3,data = df.graf,alpha=1/10) +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "blue",pch=17,cex=2.5,data = coords.stations) +
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_text(aes(x = longitude,
                y = latitude-.2,
                label = codigo), size = 2.5,data = coords.stations)

##

map <-
  get_map(location = c(lon = -44.80567, lat = -18.4764), zoom = 6,
          scale = "auto", maptype = c("satellite"), source = c("google"), 
          force = ifelse(source == "google", TRUE, TRUE),
          messaging = FALSE, urlonly = FALSE, filename = "ggmapTemp",
          crop = TRUE, color = c("color"), language = "pt-BR")
mapPoints <- ggmap(map)
mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "grey80",pch=22,cex=2.3,data = df.graf,alpha=1/10) +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "blue",pch=17,cex=2.5,data = coords.stations) +
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_text(aes(x = longitude,
                y = latitude-.2,
                label = codigo), size = 2.5,data = coords.stations)

###

map <-
  get_map(location = c(lon = -44.80567, lat = -18.4764), zoom = 6,
          scale = "auto", maptype = c("roadmap"), source = c("google"), 
          force = ifelse(source == "google", TRUE, TRUE),
          messaging = FALSE, urlonly = FALSE, filename = "ggmapTemp",
          crop = TRUE, color = c("color"), language = "pt-BR")
mapPoints <- ggmap(map)
mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "black",pch=22,cex=2.3,data = df.graf,alpha=1/10) +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "blue",pch=17,cex=2.5,data = coords.stations) +
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_text(aes(x = longitude,
                y = latitude-.2,
                label = codigo), size = 2.5,data = coords.stations)

###

map <-
  get_map(location = c(lon = -44.80567, lat = -18.4764), zoom = 6,
          scale = "auto", maptype = c("hybrid"), source = c("google"), 
          force = ifelse(source == "google", TRUE, TRUE),
          messaging = FALSE, urlonly = FALSE, filename = "ggmapTemp",
          crop = TRUE, color = c("color"), language = "pt-BR")
mapPoints <- ggmap(map)
mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "grey80",pch=22,cex=2.3,data = df.graf,alpha=1/10) +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "blue",pch=17,cex=2.5,data = coords.stations) +
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_text(aes(x = longitude,
                y = latitude-.2,
                label = codigo), size = 2.5,data = coords.stations)

###

# GRAFICO -  MAPAS DE PREVISAO RGOOGLEMAPS ----

datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")
data.c = 8
load(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/mapa2",
            datas[data.c],"DGOP.Rda"))

#------ Ajustando DF ------#
df.graf <- data.frame(longitude = gradeMG$longitude,
                      latitude = gradeMG$latitude,
                      hora18 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                        y = coords.ETA$latitude,
                                                        z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .5),
                                                        xo = gradeMG$longitude,
                                                        yo = gradeMG$latitude)$z),
                      hora18IC1 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .1),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora18IC2 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .9),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora18ME = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                          y = coords.ETA$latitude,
                                                          z = (rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .9)-
                                                                 rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,7]), probs = .1))/2,
                                                          xo = gradeMG$longitude,
                                                          yo = gradeMG$latitude)$z),
                      hora00 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                        y = coords.ETA$latitude,
                                                        z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .5),
                                                        xo = gradeMG$longitude,
                                                        yo = gradeMG$latitude)$z),
                      hora00IC1 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .1),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora00IC2 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .9),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora00ME = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                          y = coords.ETA$latitude,
                                                          z = (rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .9)-
                                                                 rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,13]), probs = .1))/2,
                                                          xo = gradeMG$longitude,
                                                          yo = gradeMG$latitude)$z),
                      hora06 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                        y = coords.ETA$latitude,
                                                        z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .5),
                                                        xo = gradeMG$longitude,
                                                        yo = gradeMG$latitude)$z),
                      hora06IC1 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .1),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora06IC2 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .9),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora06ME = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                          y = coords.ETA$latitude,
                                                          z = (rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .9)-
                                                                 rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,19]), probs = .1))/2,
                                                          xo = gradeMG$longitude,
                                                          yo = gradeMG$latitude)$z),
                      hora12 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                        y = coords.ETA$latitude,
                                                        z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .5),
                                                        xo = gradeMG$longitude,
                                                        yo = gradeMG$latitude)$z),
                      hora12IC1 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .1),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora12IC2 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                           y = coords.ETA$latitude,
                                                           z = rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .9),
                                                           xo = gradeMG$longitude,
                                                           yo = gradeMG$latitude)$z),
                      hora12ME = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                          y = coords.ETA$latitude,
                                                          z = (rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .9)-
                                                                 rowQuantiles(x = as.matrix(Yt.prev[-(1:59),,25]), probs = .1))/2,
                                                          xo = gradeMG$longitude,
                                                          yo = gradeMG$latitude)$z),
                      ETA18 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                       y = coords.ETA$latitude,
                                                       z = Ft.numforecast[2,-(1:59),7],
                                                       xo = gradeMG$longitude,
                                                       yo = gradeMG$latitude)$z),
                      ETA00 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                       y = coords.ETA$latitude,
                                                       z = Ft.numforecast[2,-(1:59),13],
                                                       xo = gradeMG$longitude,
                                                       yo = gradeMG$latitude)$z),
                      ETA06 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                       y = coords.ETA$latitude,
                                                       z = Ft.numforecast[2,-(1:59),19],
                                                       xo = gradeMG$longitude,
                                                       yo = gradeMG$latitude)$z),
                      ETA12 = suppressWarnings(interpp(x = coords.ETA$longitude,
                                                       y = coords.ETA$latitude,
                                                       z = Ft.numforecast[2,-(1:59),25],
                                                       xo = gradeMG$longitude,
                                                       yo = gradeMG$latitude)$z))

#------ Criando SPDF ------#
df.graf2 <- df.graf
coordinates(df.graf2) = c("longitude", "latitude") # promote to SpatialPointsDataFrame
gridded(df.graf2) <- TRUE # promote to SpatialPixelsDataFrame

map <-
  get_map(location = c(lon = -44.80567, lat = -18.4764), zoom = 6,
          scale = "auto", maptype = c("terrain"), source = c("google"), 
          force = ifelse(source == "google", TRUE, TRUE),
          messaging = FALSE, urlonly = FALSE, filename = "ggmapTemp",
          crop = TRUE, color = c("color"), language = "pt-BR")

mapPoints <- ggmap(map)

mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "black",pch=22,cex=2.3,data = df.graf,alpha=1/10) +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  geom_point(aes(x = longitude,
                 y = latitude), 
             col = "blue",pch=17,cex=2.5,data = coords.stations) +
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_text(aes(x = longitude,
                y = latitude-.2,
                label = codigo), size = 2.5,data = coords.stations)


map <-
  get_map(location = c(lon = -44.80567, lat = -18.4764), zoom = 6,
          scale = "auto", maptype = c("terrain"), source = c("google"), 
          force = ifelse(source == "google", TRUE, TRUE),
          messaging = FALSE, urlonly = FALSE, filename = "ggmapTemp",
          crop = TRUE, color = c("color"), language = "pt-BR")

mapPoints <- ggmap(map)

 # ggplot(df.graf, aes(longitude, latitude)) +
   # mapPoints +
   # geom_raster(aes(fill = hora18),alpha=1,data=df.graf) +
   # scale_fill_gradientn(breaks = seq(0,9,.5),limits=c(0,9),
   #                        colours = colorRampPalette(brewer.pal(9,"YlGnBu"))(18),
   #                        name = "m/s", guide="legend") +
 
mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude,
                 colour = hora18), data = df.graf,
             alpha=1/5) +
  scale_colour_gradientn(breaks = seq(0,9,.5),limits=c(0,9),
                         colours = colorRampPalette(brewer.pal(9,"YlGnBu"))(18),
                         name = "m/s") +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.key.height = unit(2,"cm"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))

mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude,
                 colour = hora18), data = df.graf,
             alpha=1/5) +
  scale_colour_gradientn(breaks = seq(0,9,.5),limits=c(0,9),
                         colours = colorRampPalette(brewer.pal(9,"YlGnBu"))(18),
                         name = "m/s") +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.key.height = unit(2,"cm"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))

mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude,
                 colour = ETA18), data = df.graf,
             alpha=1/5) +
  scale_colour_gradientn(breaks = seq(0,9,.5),limits=c(0,9),
                         colours = colorRampPalette(brewer.pal(9,"YlGnBu"))(18),
                         name = "m/s") +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.key.height = unit(2,"cm"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))

mapPoints +
  geom_point(aes(x = longitude,
                 y = latitude,
                 colour = hora18ME), data = df.graf,
             alpha=1/15) +
  scale_colour_gradientn(breaks = seq(0,5,1),limits=c(0,5),
                         colours = colorRampPalette(brewer.pal(5,"BuPu"))(10),
                         name = "m/s") +
  geom_path(aes(x = V1,
                y = V2),
            col = "black", data = coords.MG) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.key.height = unit(2,"cm"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12))
     