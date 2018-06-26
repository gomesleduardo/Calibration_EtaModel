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
##########|              ALUNO: LUIZ EDUARDO S. GOMES             |###########
##########|-------------------------------------------------------|###########
##########|       CALIBRACAO DINAMICA DE PREVISOES NUMERICAS      |###########
##########|                      DYNAMIC GOP                      |###########
##########|-------------------------------------------------------|###########
##########|-------------------------------------------------------|###########
##############################################################################
##############################################################################

#----------------------------------------------------------------------------------------------
# CARREGANDO PACOTES
#----------------------------------------------------------------------------------------------

pacotes<-c("mvtnorm","forecast","coda","matrixStats",
           "corpcor","Matrix","akima","Rcpp","Rcpp11",
           "RcppArmadillo","geoR","adaptMCMC",
           "tmvtnorm","mvnfast")

for (i in 1:length(pacotes))
{
  if (length(names(installed.packages()[,1])[names(installed.packages()[,1])==pacotes[i]])==0)
  {install.packages(pacotes[i], repos="http://cran.fiocruz.br/")}
  library(pacotes[i],character.only = TRUE) 
}
rm(i,pacotes)

#----------------------------------------------------------------------------------------------
# FUNCOES
#----------------------------------------------------------------------------------------------

#################################
# Na pratica, data e hora de previsao maxima nao entrarao pois isto sera calculado
# operacionalmente a partir da ultima data disponivel + horizonte requerido
# Em geral, a qtd de novas localizacoes tambem sera fixa. Aqui serve so pra exemplificar.
# SAIDA: matrizes necessarias para os calculos (Yt, Ft, Gt, etc... )
#################################

definicoes.usuario <- function(forecast.date, # data de previsao maxima
                               forecast.hour, # hora da previsao maxima
                               horizon, # horizonte
                               lag.hours) # tamanho do periodo de treinamento (janela passada)
{
  #------ Contador de tempo ------#
  a <- Sys.time()
  
  #---------------------------------------- LENDO DADOS ------------------------------------------
  obs.windCE <- read.csv("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_100mCE/dados_parajuru/velvento85m_aerogeradores_parajuru_20170101_20171231.csv")
  obs.windCE <- data.frame(data.hora = as.character(strptime(paste0(substr(obs.windCE$data,7,10),"-",
                                                                    substr(obs.windCE$data,4,5),"-",
                                                                    substr(obs.windCE$data,1,2)," ",
                                                                    obs.windCE$hora,":01"),
                                                             format = "%Y-%m-%d %H:%M:%S",
                                                             tz = "GMT")),
                           obs.windCE[,-(1:2)])
  obs.windCE <- subset(obs.windCE,
                       substr(obs.windCE$data.hora,15,16)=="00" |
                         substr(obs.windCE$data.hora,15,16)=="30")
  obs.windtorreCE <- read.csv("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_100mCE/dados_parajuru/dados_torre_energimp_parajuru_30min_20170101_20171231.csv")
  obs.windtorreCE <- data.frame(data.hora = as.character(strptime(paste0(as.character(obs.windtorreCE$data)," ",
                                                                         ifelse(nchar(obs.windtorreCE$hora)==2,
                                                                                obs.windtorreCE$hora,
                                                                                paste0("0",obs.windtorreCE$hora)),":",
                                                                         ifelse(nchar(obs.windtorreCE$minutos)==2,
                                                                                obs.windtorreCE$minutos,
                                                                                paste0("0",obs.windtorreCE$minutos)),
                                                                         ":01"),
                                                                  format = "%Y-%m-%d %H:%M:%S",
                                                                  tz = "GMT")),
                                torre = obs.windtorreCE$vel85)
  obs.windCE <- merge(x = obs.windCE,y = obs.windtorreCE,by = "data.hora",all = TRUE)
  rm(obs.windtorreCE)
  obs.windCE$data.hora <- as.character(obs.windCE$data.hora)
  coords.windfarm <- read.csv2("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_100mCE/dados_parajuru/coords_parque_praias_de_parajuru_energimp.csv")
  coords.windfarm[,1] <- as.character(coords.windfarm[,1])
  coords.windfarm[,2] <- as.numeric(as.character(coords.windfarm[,2]))
  coords.windfarm[,3] <- as.numeric(as.character(coords.windfarm[,3]))
  names(coords.windfarm)[2:3] <- tolower(names(coords.windfarm)[2:3])
  locat.aux <- NULL
  for (l in 2:length(names(obs.windCE)))
  {
    locat.aux <- c(locat.aux,
                   which(coords.windfarm$AG==names(obs.windCE)[l]))
  }
  coords.windfarm <- coords.windfarm[locat.aux,]
  rm(locat.aux,l)
  coords.ETACE <- read.csv("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_100mCE/coords_Eta5km_operacional_CE.csv")
  prevnumext <- read.csv("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_100mCE/dados_parajuru/previsoes_numericas_parajuru_20170102_20171229.csv")
  prevnumext$runtime1 <- strptime(paste0(substr(as.character(prevnumext$runtime1),1,4),"-",
                                         substr(as.character(prevnumext$runtime1),5,6),"-",
                                         substr(as.character(prevnumext$runtime1),7,8)," ",
                                         substr(as.character(prevnumext$runtime1),9,10),":",
                                         substr(as.character(prevnumext$runtime1),11,12),":01"),
                                  format = "%Y-%m-%d %H:%M:%S",
                                  tz = "GMT")
  prevnumext$prevtime1 <- strptime(paste0(substr(as.character(prevnumext$prevtime1),1,4),"-",
                                          substr(as.character(prevnumext$prevtime1),5,6),"-",
                                          substr(as.character(prevnumext$prevtime1),7,8)," ",
                                          substr(as.character(prevnumext$prevtime1),9,10),":",
                                          substr(as.character(prevnumext$prevtime1),11,12),":01"),
                                   format = "%Y-%m-%d %H:%M:%S",
                                   tz = "GMT")
  
  prevnumext <- data.frame(runtime1 = prevnumext$runtime1,
                           prevtime1 = prevnumext$prevtime1,
                           horizonte = as.numeric((prevnumext$prevtime1-prevnumext$runtime1)/60/60),
                           prevnumext[,-(1:2)])
  
  prevnumext$runtime1 <- as.character(prevnumext$runtime1)
  prevnumext$prevtime1 <- as.character(prevnumext$prevtime1)
  
  #------ Periodo de treinamento ------#
  train.period <- as.character(seq(strptime(paste0(forecast.date," ",
                                                   forecast.hour,":00:01"),
                                            format = "%Y-%m-%d %H:%M:%S",
                                            tz = "GMT")-(horizon-24)*60*60-30*60-(lag.hours-1)*30*60,
                                   strptime(paste0(forecast.date," ",
                                                   forecast.hour,":00:01"),
                                            format = "%Y-%m-%d %H:%M:%S",
                                            tz = "GMT")-(horizon-24)*60*60-30*60,30*60))
  
  # #------ Recuperando as informacoes nos arquivos .csv ------#
  # #----------------------------------------------------------------------------------------------
  # # VELOCIDADE DO VENTO - 10 m
  # #----------------------------------------------------------------------------------------------
  # wind10m.ETA <- data.frame(NULL)
  # aux.date <- seq(as.Date(head(unique(substr(train.period,1,10)),n=1))-9,
  #                 as.Date(tail(unique(substr(train.period,1,10)),n=1)), by = "1 day")
  # 
  # #------ Removendo dias que nao houve funcionamento do Eta ------#
  # aux.date.rm <- NULL
  # if (sum(aux.date < as.Date("2015-10-11"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date < as.Date("2015-10-11")))
  # if (sum(aux.date > as.Date("2017-09-28"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date > as.Date("2017-09-28")))
  # if (sum(aux.date == as.Date("2015-11-09"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-11-09")))
  # if (sum(aux.date == as.Date("2015-12-08"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-12-08")))
  # if (sum(aux.date == as.Date("2015-12-25"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-12-25")))
  # if (sum(aux.date == as.Date("2016-01-10"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-01-10")))
  # if (sum(aux.date == as.Date("2016-03-21"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-03-21")))
  # if (sum(aux.date == as.Date("2016-08-07"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-08-07")))
  # if (sum(aux.date == as.Date("2017-01-02"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-01-02")))
  # if (sum(aux.date == as.Date("2017-02-05"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-02-05")))
  # if (sum(aux.date == as.Date("2017-03-28"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-03-28")))
  # if (sum(aux.date == as.Date("2017-09-13"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-09-13")))
  # if (length(aux.date.rm)>0) aux.date <- aux.date[-aux.date.rm]
  # 
  # for (i in 1:length(aux.date))
  # {
  #   print(paste0(i,"/",length(aux.date)))
  #   #------ lendo .csv auxiliar ------#
  #   aux.csv <- read.csv(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/velvento10m/",
  #                              "velvento10m_Eta15km_MG_",
  #                              substr(aux.date[i],1,4),
  #                              substr(aux.date[i],6,7),
  #                              substr(aux.date[i],9,10),
  #                              "12.csv"))
  #   
  #   #------ organizando dia/hora origem ------#
  #   aux.csv[,1] <- as.character(strptime(paste0(substr(aux.csv[,1],1,4),"-",
  #                                               substr(aux.csv[,1],5,6),"-",
  #                                               substr(aux.csv[,1],7,8)," ",
  #                                               substr(aux.csv[,1],9,10),":00:01"),
  #                                        format = "%Y-%m-%d %H:%M:%S",
  #                                        tz = "GMT"))
  #   
  #   #------ organizando dia/hora data de previsao ------#
  #   aux.csv[,2] <- as.character(strptime(paste0(substr(aux.csv[,2],1,4),"-",
  #                                               substr(aux.csv[,2],5,6),"-",
  #                                               substr(aux.csv[,2],7,8)," ",
  #                                               substr(aux.csv[,2],9,10),":00:01"),
  #                                        format = "%Y-%m-%d %H:%M:%S",
  #                                        tz = "GMT"))
  #   
  #   #------ criando a variavel "horizonte" ------#
  #   aux.csv <- data.frame(aux.csv[,1:2],
  #                         horizonte = as.numeric(difftime(time1 = strptime(aux.csv[,2], format = "%Y-%m-%d %H:%M:%S",tz = "GMT"),
  #                                                         time2 = strptime(aux.csv[,1], format = "%Y-%m-%d %H:%M:%S",tz = "GMT"),
  #                                                         tz = "GMT", units = "hours")),
  #                         aux.csv[,-(1:2)])
  #   
  #   #------ mantendo apenas o horario desejado e o horizonte de interesse ------#
  #   aux.csv <- subset(aux.csv,aux.csv$horizonte >= horizon)
  #   
  #   #------ interpolando para as estacoes (interpolacao bilinear) ------#
  #   for (l in 1:dim(aux.csv)[1])
  #   {
  #     wind10m.ETA <- rbind(wind10m.ETA,data.frame(t(c(as.character(aux.csv[l,1:3]),
  #                                                     suppressWarnings(interpp(x = coords.ETA$longitude,
  #                                                                              y = coords.ETA$latitude,
  #                                                                              z = as.numeric(aux.csv[l,-(1:3)]),
  #                                                                              xo = coords.stations$longitude,
  #                                                                              yo = coords.stations$latitude)$z)))))
  #   }
  # }
  # 
  # #------ formatando as variaveis ------#
  # wind10m.ETA[,1] <- as.character(wind10m.ETA[,1])
  # wind10m.ETA[,2] <- as.character(wind10m.ETA[,2])
  # for (i in 3:dim(wind10m.ETA)[2]){wind10m.ETA[,i] <- as.numeric(as.character(wind10m.ETA[,i]))}
  # 
  # #------ nomeando as variaveis ------#
  # names(wind10m.ETA) <- c(names(aux.csv)[1:3],as.character(coords.stations$codigo))
  
  #------ Matriz de distancias D ------#
  D <- as.matrix(dist(data.frame(coords.windfarm$longitude,coords.windfarm$latitude)))
  D <- matrix(as.numeric(D),dim(D)[1],dim(D)[2])
  
  #------ Matriz de evolucao Gt ------#
  Gt <- as.matrix(bdiag(diag(1,1),
                        matrix(c(cos(2*pi*1/48),-sin(2*pi*1/48),
                                 sin(2*pi*1/48),cos(2*pi*1/48)),2)))
  Gt <- matrix(as.numeric(Gt),dim(Gt)[1],dim(Gt)[2])
  # Gt = matriz bloco diagonal: diag(qtd de cov + 1 - 2) + harmonicos (fc de sen e cos)
  
  #------ Comprimento do vetor Theta (Tendencia espacial + Sazonalidade) ------#
  r <- dim(Gt)[1]
  
  #------ Total de tempos ------#
  T <- length(train.period)
  
  #------ Total de localizacoes ------#
  n <-  dim(coords.windfarm)[1]
  
  #------ Matriz de covariaveis Ft ------#
  #------ Valor observado ------#
  Yt <- array(NA,dim=c(n,1,T+1))
  Ft <- array(NA,dim=c(r,n,T+1))
  for (t in 2:(T+1))
  {
    Yt[,,t] <- as.numeric(subset(obs.windCE, obs.windCE$data.hora == train.period[t-1])[,-1])
    
    Ft[,,t] <- t(data.frame(# A restricao do horizonte de previsao nas previsoes passadas nao pode ser cumprida aqui
      # Nao ha previsao sab, dom e feriados. Assim, o jeito e usar sem restricao
      ensemble.mean = rep(as.numeric(apply(subset(prevnumext, 
                                                  prevnumext$prevtime1 == train.period[t-1],
                                                  "wind_speed"),2,mean)),n),
      saz1 = 1,
      saz2 = 0))
  }
  
  #------ Periodo de previsao ------#
  forecast.period <- as.character(seq(strptime(paste0(forecast.date," ",
                                                      forecast.hour,":00:01"),
                                               format = "%Y-%m-%d %H:%M:%S",
                                               tz = "GMT")-(horizon-24)*60*60,
                                      strptime(paste0(forecast.date," ",
                                                      forecast.hour,":00:01"),
                                               format = "%Y-%m-%d %H:%M:%S",
                                               tz = "GMT"),30*60))
  
  # #------ Recuperando as informacoes nos arquivos .csv ------#
  # #----------------------------------------------------------------------------------------------
  # # VELOCIDADE DO VENTO - 10 m
  # #----------------------------------------------------------------------------------------------
  # 
  # wind10m.ETA2 <- data.frame(NULL)
  # aux.date <- seq(as.Date(head(unique(substr(forecast.period,1,10)),n=1))-9,
  #                 as.Date(tail(unique(substr(forecast.period,1,10)),n=1))-1, by = "1 day")
  # 
  # #------ Removendo dias que nao houve funcionamento do Eta ------#
  # aux.date.rm <- NULL
  # if (sum(aux.date < as.Date("2015-10-11"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date < as.Date("2015-10-11")))
  # if (sum(aux.date > as.Date("2017-09-28"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date > as.Date("2017-09-28")))
  # if (sum(aux.date == as.Date("2015-11-09"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-11-09")))
  # if (sum(aux.date == as.Date("2015-12-08"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-12-08")))
  # if (sum(aux.date == as.Date("2015-12-25"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-12-25")))
  # if (sum(aux.date == as.Date("2016-01-10"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-01-10")))
  # if (sum(aux.date == as.Date("2016-03-21"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-03-21")))
  # if (sum(aux.date == as.Date("2017-01-02"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-01-02")))
  # if (sum(aux.date == as.Date("2017-02-05"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-02-05")))
  # if (sum(aux.date == as.Date("2017-03-28"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-03-28")))
  # if (sum(aux.date == as.Date("2017-09-13"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-09-13")))
  # if (length(aux.date.rm)>0) aux.date <- aux.date[-aux.date.rm]
  # 
  # for (i in 1:length(aux.date))
  # {
  #   print(paste0(i,"/",length(aux.date)))
  #   #------ lendo .csv auxiliar
  #   aux.csv <- read.csv(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/velvento10m/",
  #                              "velvento10m_Eta15km_MG_",
  #                              substr(aux.date[i],1,4),
  #                              substr(aux.date[i],6,7),
  #                              substr(aux.date[i],9,10),
  #                              "12.csv"))
  #   
  #   #------ organizando dia/hora origem ------#
  #   aux.csv[,1] <- as.character(strptime(paste0(substr(aux.csv[,1],1,4),"-",
  #                                               substr(aux.csv[,1],5,6),"-",
  #                                               substr(aux.csv[,1],7,8)," ",
  #                                               substr(aux.csv[,1],9,10),":00:01"),
  #                                        format = "%Y-%m-%d %H:%M:%S",
  #                                        tz = "GMT"))
  #   
  #   #------ organizando dia/hora data de previsao ------#
  #   aux.csv[,2] <- as.character(strptime(paste0(substr(aux.csv[,2],1,4),"-",
  #                                               substr(aux.csv[,2],5,6),"-",
  #                                               substr(aux.csv[,2],7,8)," ",
  #                                               substr(aux.csv[,2],9,10),":00:01"),
  #                                        format = "%Y-%m-%d %H:%M:%S",
  #                                        tz = "GMT"))
  #   
  #   #------ criando a variavel "horizonte" ------#
  #   aux.csv <- data.frame(aux.csv[,1:2],
  #                         horizonte = as.numeric(difftime(time1 = strptime(aux.csv[,2], format = "%Y-%m-%d %H:%M:%S",tz = "GMT"),
  #                                                         time2 = strptime(aux.csv[,1], format = "%Y-%m-%d %H:%M:%S",tz = "GMT"),
  #                                                         tz = "GMT", units = "hours")),
  #                         aux.csv[,-(1:2)])
  #   
  #   #------ mantendo apenas o horario desejado e o horizonte de interesse ------#
  #   aux.csv <- subset(aux.csv,aux.csv$horizonte >= horizon)
  #   
  #   #------ interpolando para as estacoes (interpolacao bilinear) ------#
  #   for (l in 1:dim(aux.csv)[1])
  #   {
  #     wind10m.ETA2 <- rbind(wind10m.ETA2,data.frame(t(c(as.character(aux.csv[l,1:3]),
  #                                                       suppressWarnings(interpp(x = coords.ETA$longitude,
  #                                                                                y = coords.ETA$latitude,
  #                                                                                z = as.numeric(aux.csv[l,-(1:3)]),
  #                                                                                xo = coords.stations$longitude,
  #                                                                                yo = coords.stations$latitude)$z),
  #                                                       as.character(aux.csv[l,4:(qgrid.new+3)])))))
  #   }
  # }
  # 
  # #------ formatando as variaveis ------#
  # wind10m.ETA2[,1] <- as.character(wind10m.ETA2[,1])
  # wind10m.ETA2[,2] <- as.character(wind10m.ETA2[,2])
  # for (i in 3:dim(wind10m.ETA2)[2]){wind10m.ETA2[,i] <- as.numeric(as.character(wind10m.ETA2[,i]))}
  # 
  # #------ nomeando as variaveis ------#
  # names(wind10m.ETA2) <- c(names(aux.csv)[1:3],as.character(coords.stations$codigo),as.character(coords.ETA$ID)[1:qgrid.new])
  
  #------ Matriz de distancias D ------#
  # D.new <- as.matrix(dist(data.frame(c(coords.stations$longitude,coords.ETA$longitude[1:qgrid.new]),
  #                                    c(coords.stations$latitude,coords.ETA$latitude[1:qgrid.new]))))
  # D.new <- matrix(as.numeric(D.new),dim(D.new)[1],dim(D.new)[2])
  D.new <- D
  
  #------ Total de localizacoes - Completo ------#
  n.new <-  dim(D.new)[1]
  
  #------ Previsoes a frente ------#
  T.new <- length(forecast.period)
  
  #------ Matriz de covariaveis Ft ------#
  #------ Valor observado ------#
  Yt.k <- array(NA,dim=c(n,1,T.new+1))
  Ft.new <- array(NA,dim=c(r,n.new,T.new+1))
  for (t in 2:(T.new+1))
  {
    Yt.k[,,t] <- as.numeric(subset(obs.windCE, obs.windCE$data.hora == forecast.period[t-1])[,-1])
    
    Ft.new[,,t] <- t(data.frame(# A restricao do horizonte de previsao nas previsoes passadas nao pode ser cumprida aqui
      # Nao ha previsao sab, dom e feriados. Assim, o jeito e usar sem restricao
      ensemble.mean = rep(as.numeric(apply(subset(prevnumext, 
                                                  prevnumext$prevtime1 == forecast.period[t-1],
                                                  "wind_speed"),2,mean)),n),
      saz1 = 1,
      saz2 = 0))
  }
  
  resultado <- list(Yt = Yt,
                    Ft = Ft,
                    Gt = as.matrix(Gt),
                    D = as.matrix(D),
                    Ft.k = Ft.new,
                    D.full = as.matrix(D.new),
                    relatorio = data.frame(codigo = coords.windfarm$AG,
                                           longitude = coords.windfarm$longitude,
                                           latitude = coords.windfarm$latitude),
                    forecast.period = forecast.period,
                    Yt.k = Yt.k)
  
  #------ Tempo total de processamento ------#
  b <- Sys.time()-a
  print(b)
  
  return(resultado)
}

#################################
# Aqui e onde os parametros sao estimados por metodos computacionais intensivos
# O total de repeticoes e calculado da seguinte forma: T = M*salto
# Se desejo uma amostra posteriori com 2000 e um burnin de 500, entao:
# M = 2500, burnin = 500
# O salto depende da autocorrelacao da cadeia MCMC
# Nesta configuracao, se o salto for 10, entao
# Sera repetido o MCMC M*salto = 2500*10 = 25000 vezes
#################################

MCMC <- function(M, # Qtd de amostras a posteriori (+ burnin)
                 burnin, # Qtd de amostras a posteriori descartadas inicialmente
                 salto, # Espacamento das cadeias MCMC (thinning)
                 Desc, # Vetor de descontos (comprimento = qtd de covariaveis + 1,
                 # sendo os 2 valores finais tendo de ser iguais sob esta configuracao)
                 delta.sig2, # Fator de desconto para o sigma2
                 Yt, # Vetor observado
                 Ft, # Matriz de design
                 Gt, # Matriz de evolucao
                 Dist) # Matriz de distancia
{
  #------ Contador de tempo ------#
  a <- Sys.time()
  
  #------ Valores missing ------#
  YtNA <- which(is.na(Yt),arr.ind = TRUE)
  YtNA <- subset(YtNA,YtNA[,3]!=1)
  
  #------ Inputando valores ------#
  Yt[YtNA] <- 1
  
  #---------------------#
  #------ Funcoes ------#
  #---------------------#
  
  #------ Fc M-H adapt ------#
  MHadapt <- function (p, n, init, scale = rep(1, length(init)), adapt = !is.null(acc.rate), 
                       acc.rate = NULL, gamma = 2/3, list = TRUE, showProgressBar = interactive(), 
                       n.start = 0, ...) 
  {
    # if (adapt & !is.numeric(acc.rate)) 
    #   stop("Argument acc.rate is missing!")
    # if (gamma <= 0.5 | gamma > 1) 
    #   stop("Argument gamma must be in (0.5, 1]!")
    
    # parar a adaptacao depois de algumas iterações
    if (is.numeric(adapt)) 
      n.adapt <- adapt
    if (adapt == TRUE) 
      n.adapt <- Inf
    if (adapt == FALSE) 
      n.adapt <- 0
    
    # dimensao
    d <- length(init)
    
    # matriz vazia com dimencao d x n
    X <- matrix(NA, ncol = d, nrow = n)
    
    # nome das variaveis
    colnames(X) <- names(init)
    
    # vetor inicial
    X[1, ] <- init
    
    # vetor com a log aplicada salva
    p.val <- rep(NA, n)
    
    # aplica a funcao no vetor inicial
    # val <- p(X[1, ], ...)
    val <- p(X[1, ])
    
    # caso tenha mais argumentos
    if (is.list(val)) {
      returns.list <- TRUE
      extras <- list()
      if (!"log.density" %in% names(val)) {
        stop("The list returned by 'p' must contain an element named 'log.density!'")
      }
      if (length(val$log.density) > 1) 
        stop("The list element 'log.density' must be a scalar value!")
      p.val[1] <- val$log.density
      extras[[1]] <- val["log.density" != names(val)]
    }
    else {
      returns.list <- FALSE
      if (length(val) > 1) 
        stop("The function 'p' must return a scalar value or a named list! See ?MCMC.!")
      p.val[1] <- val
    }
    
    # matriz de escala
    # se a dimensao e maior que 1 (as contas mudam caso d>1)
    if (d > 1) {
      if (length(scale) == d) {
        M <- diag(scale)
      }
      else {
        M <- scale
      }
    }
    else {
      M <- matrix(scale)
    }
    
    # # erro da matriz de escala
    # if (ncol(M) != length(init)) 
    #   stop("Length or dimension of 'init' and 'scale' do not match!")
    
    # contas !!!!!!!!!!!!!
    
    # S e a matriz que ajusta a taxa de aceitacao
    # cholesky da matriz de escala
    S <- t(chol(M))
    
    # cat("  generate", n, "samples \\n")
    # # barra de progresso
    # if (showProgressBar) {
    #   pb <- txtProgressBar(min = 0, max = n, style = 3)
    # }
    
    # print na tela
    update.step <- max(5, floor(n/100))
    
    # contador de aceitacao
    k <- 0
    
    # loop
    for (i in 2:n) {
      
      # # mostrar barra de progresso
      # if (showProgressBar && i%%update.step == 0) {
      #   setTxtProgressBar(pb, i)
      # }
      
      # amostra de uma t com d graus de liberdade
      # U <- rt(d, df = d)
      U <- rnorm(d)
      
      # valor proposto
      X.prop <- c(X[i - 1, ] + S %*% U)
      
      # nomeando
      names(X.prop) <- names(init)
      
      # aplicando na fc de log vero
      # val <- p(X.prop, ...)
      val <- p(X.prop)
      
      # caso tenha mais argumentos
      if (returns.list) {
        p.val.prop <- val$log.density
        extras.prop <- val["log.density" != names(val)]
      }
      else {
        # densidade proposta
        p.val.prop <- val
      }
      
      # probabilidade de aceitacao
      alpha <- min(1, exp(p.val.prop - p.val[i - 1]))
      
      # caso seja infinito
      if (!is.finite(alpha)) 
        alpha <- 0
      
      # aceitacao - rejeicao
      if (runif(1) < alpha) {
        # X recebe valor proposto
        X[i, ] <- X.prop
        # p.val recebe densidade no valor proposto
        p.val[i] <- p.val.prop
        
        # caso tenha mais argumentos
        if (returns.list) {
          extras[[i]] <- extras.prop
        }
        # contador de aceitacao
        k <- k + 1
      }
      else {
        X[i, ] <- X[i - 1, ]
        p.val[i] <- p.val[i - 1]
        if (returns.list) {
          extras[[i]] <- extras[[i - 1]]
        }
      }
      
      # adaptacao
      # ii e onde a adapatcao comeca
      ii <- i + n.start
      
      # n.adapt e onde a adaptacao termina
      if (ii < n.adapt) {
        # taxa de adaptacao
        adapt.rate <- min(1, d * ii^(-gamma))
        
        # matriz de escala M recebe
        M <- S %*% (diag(d) + adapt.rate * (alpha - acc.rate) * 
                      U %*% t(U)/sum(U^2)) %*% t(S)
        
        # Com o passar das iteracoes, adapt.rate -> 0
        # Se a probabilidade de aceitacao alpha for 0.234, a matriz de escala nao muda
        
        # autovalores da matriz escala
        eig <- eigen(M, only.values = TRUE)$values
        
        # caso a matriz nao seja simetrica
        tol <- ncol(M) * max(abs(eig)) * .Machine$double.eps
        if (!isSymmetric(M) | is.complex(eig) | !all(Re(eig) > 
                                                     tol)) {
          M <- as.matrix(Matrix::nearPD(M)$mat)
        }
        
        # S recebe cholesky de M
        S <- t(chol(M))
      }
    }
    
    # fechar conexao da barra  
    # if (showProgressBar) {
    #   close(pb)
    # }
    
    # taxa de aceitacao 
    acceptance.rate <- round(k/(n - 1), 3)
    
    # juntando resultado
    if (list) {
      res <- list(samples = X, log.p = p.val, cov.jump = M, 
                  n.sample = n, acceptance.rate = acceptance.rate, 
                  adaption = adapt, sampling.parameters = list(sample.density = p, 
                                                               acc.rate = acc.rate, gamma = gamma))
      if (returns.list) {
        res$extra.values = extras
      }
      return(res)
    }
    else {
      cat("Acceptance rate:", acceptance.rate, "\\n")
      return(X)
    }
  }
  
  MHadapt.add <- function (MCMC.object, n.update, ...) 
  {
    if (!is.null(names(MCMC.object))) {
      if (is.matrix(MCMC.object)) 
        stop("Only MCMC objects generated with option 'list=TRUE' can be updated!")
      p <- MCMC.object$sampling.parameters$sample.density
      init <- MCMC.object$samples[nrow(MCMC.object$samples), 
                                  ]
      scale <- MCMC.object$cov.jump
      acc.rate <- MCMC.object$sampling.parameters$acc.rate
      gamma <- MCMC.object$sampling.parameters$gamma
      n.before <- MCMC.object$n.sample
      adapt <- MCMC.object$adaption
      samp.update <- MHadapt(p = p, n = n.update, init = init, 
                             scale = scale, adapt = adapt, acc.rate = acc.rate, 
                             gamma = gamma, list = TRUE, n.start = n.before, ...)
      MCMC.object$cov.jump <- samp.update$cov.jump
      m <- c(MCMC.object$n.sample, samp.update$n.sample)
      MCMC.object$acceptance.rate <- 1/sum(m) * (m[1] * MCMC.object$acceptance.rate + 
                                                   m[2] * samp.update$acceptance.rate)
      MCMC.object$n.sample <- MCMC.object$n.sample + n.update
      MCMC.object$samples <- rbind(MCMC.object$samples, samp.update$samples)
      MCMC.object$log.p <- c(MCMC.object$log.p, samp.update$log.p)
      if ("extra.values" %in% names(MCMC.object)) {
        MCMC.object$extra.values <- c(MCMC.object$extra.values, 
                                      samp.update$extra.values)
      }
      return(MCMC.object)
    }
    if (is.null(names(MCMC.object))) {
      MCMC.object <- lapply(MCMC.object, function(x) MCMC.add.samples(x, 
                                                                      n.update = n.update, ...))
      return(MCMC.object)
    }
  }
  
  #------ Total de tempos ------#
  T <- dim(Ft)[3]-1
  
  #------ Total de covariaveis + 1 ------#
  r <- dim(Ft)[1]
  
  #------ Total de localizacoes ------#
  n <- dim(Ft)[2]
  
  #------ valor da censura ------#
  csrd <- 0
  
  #------ Matriz de desconto ------#
  Desc <- diag(1/sqrt(Desc))
  
  #------------------------------------#
  #------ Distribuicoes a priori ------#
  #------------------------------------#
  
  # #------ phi ------#
  # a.phi.prior <- 2
  # b.phi.prior <- max(Dist)/6
  # log.phi.prior <- function(x){dgamma(x,shape = a.phi.prior,rate = b.phi.prior,log = T)}
  
  # #------ lambda ------#
  # mu.lambda.prior <- 1
  # sigma.lambda.prior <- 1
  # log.lambda.prior <- function(x){dnorm(x,mean = mu.lambda.prior,sd = sqrt(sigma.lambda.prior),log = T)}
  
  #------------------------------------------#
  #------ Valores iniciais das cadeias ------#
  #------------------------------------------#
  
  #------ Thetat final ------#
  Thetat.chain.final <- array(NA, dim = c(T+1,r,M))
  
  #------ mt e Ct final ------#
  ht.final <- matrix(0,M,r)
  Ht.final <- array(diag(r), dim = c(r,r,M))
  
  #------ sigma2t final ------#
  sigma2t.chain.final <- matrix(NA,T+1,M)
  
  #------ phi e lambda = PSI ------#
  # PSI.chain <- data.frame(t(c(1,1)))
  # names(PSI.chain) <- c("phi","lambda")
  # PSI.chain <- data.frame(t(c(1)))
  # names(PSI.chain) <- c("lambda")
  
  #------ Matriz de correlacao da proposta M-H ------#
  # EPS.PSI <- diag(2.4/sqrt(dim(PSI.chain)[2]),dim(PSI.chain)[2])
  
  #------ Vetor de PGs latentes correntes ------#
  # Xt <- BC(Yt = Yt, lambda = PSI.chain$lambda, csrd = csrd)
  Xt <- BC(Yt = Yt, lambda = 1, csrd = csrd)
  
  #------ Matriz de correlacao corrente ------#
  # mat.cor.c <- MatCor(phi = PSI.chain$phi, DistMatrix = Dist)
  mat.cor.c <- MatCor(phi = 100000, DistMatrix = Dist)
  
  #------ Fc de verossimilhanca corrente ------#
  # p.log <- function(x)
  # {
  #   if (all(x>0))
  #   {
  #     loglikep(Yt = Yt,
  #              Ft = Ft,
  #              Thetat = auxFFBS$Thetatchain,
  #              sigma2t = auxFFBS$sigma2tchain,
  #              phi = x[1],
  #              lambda = x[2],
  #              DistMatrix = Dist,
  #              csrd = csrd) +
  #       log.phi.prior(x[1]) + log.lambda.prior(x[2])
  #   } else {-9.99e+300}
  # }
  # p.log <- function(x)
  # {
  #   if (all(x>0))
  #   {
  #     loglikep(Yt = Yt,
  #              Ft = Ft,
  #              Thetat = auxFFBS$Thetatchain,
  #              sigma2t = auxFFBS$sigma2tchain,
  #              phi = 100000,
  #              lambda = x,
  #              DistMatrix = Dist,
  #              csrd = csrd) +
  #       log.lambda.prior(x)
  #   } else {-9.99e+300}
  # }
  
  #------ Atualizacao da barra de progresso ------#
  update.step <- max(5, floor(M/100))
  
  #------ Barra de progresso ------#
  pb <- txtProgressBar(max = M, style = 3)
  
  for(i in 2:M)
  {
    #--------------------------------------------------------------#
    #------ Forward Filtering Backward Smoothing para Thetat ------#
    #--------------------------------------------------------------#
    auxFFBS <- FFBSvar(Xt = Xt,
                       Ft = Ft,
                       Gt = Gt, 
                       CorMatrix = mat.cor.c,
                       DescMatrix = Desc, 
                       delta = delta.sig2)
    
    ht.final[i,] <- auxFFBS$ht[,,T+1]
    Ht.final[,,i] <- auxFFBS$Ht[,,T+1]
    Thetat.chain.final[,,i] <- auxFFBS$Thetatchain
    Xt.mu <- auxFFBS$Xtmu
    sigma2t.chain.final[,i] <- auxFFBS$sigma2tchain
    
    #--------------------------------------------------------------------#
    #------ Metropolis-Hastings com propostas adaptativas para PSI ------#
    #--------------------------------------------------------------------#
    
    # if (i==2)
    # {
    #   sampMH <- MHadapt(p.log, n=salto, init=as.numeric(PSI.chain[i-1,]),
    #                     scale = EPS.PSI,
    #                     adapt=TRUE, acc.rate=0.234, showProgressBar=FALSE)
    # } else {
    #   sampMH <- MHadapt.add(sampMH, n.update=salto,showProgressBar=FALSE)
    # }
    # 
    # #------ Thinning ------# 
    # PSI.chain[i,] <- sampMH$samples[(i-1)*salto,]
    
    #------ Atualizacao da matriz de correlacao corrente Vt(k) ------#
    # if (all.equal(as.numeric(PSI.chain[i,]),as.numeric(PSI.chain[i-1,]))!=T)
    # {
    #   #------ Matriz de correlacao corrente ------#
    #   mat.cor.c <- MatCor(phi = PSI.chain$phi[i], DistMatrix = Dist)
    # }
    
    #--------------------------------------------------#
    #------ Amostrador de Gibbs p/ Z var latente ------#
    #--------------------------------------------------#
    for (t in 2:(T+1))
    {
      if (length(which(Yt[,,t] == csrd))>0)
      {
        cens.c <- which(Yt[,,t] == csrd)
        restr <- rep(Inf,n)
        restr[cens.c] <- rep(((csrd^PSI.chain$lambda[i])-1)/PSI.chain$lambda[i],length(cens.c))
        Zt.c <- tmvtnorm::rtmvnorm(n=1,
                                   mean = Xt.mu[,,t],
                                   sigma= auxFFBS$sigma2tchain[t]*mat.cor.c,
                                   upper=restr,
                                   algorithm = "gibbs",
                                   burn.in.samples = 0)
        # start.value = Yt[,,t])
        Yt[cens.c,,t] <- ifelse(is.nan(Zt.c[,cens.c])==F,
                                Zt.c[,cens.c],
                                Yt[cens.c,,t])
      }
    }
    
    #--------------------------------------------------#
    #------ Amostrador de Gibbs p/ U var latente ------#
    #--------------------------------------------------#
    for (t in unique(YtNA[,3]))
    {
      miss.c <- subset(YtNA,YtNA[,3]==t)[,1]
      Ut.c <- mvtnorm::rmvnorm(n = 1,
                               mean = Xt.mu[,,t],
                               sigma = auxFFBS$sigma2tchain[t]*mat.cor.c,
                               method = "chol")
      Yt[miss.c,,t] <- ifelse(is.nan(Ut.c[,miss.c])==F,
                              Ut.c[,miss.c],
                              Yt[miss.c,,t])
    }
    
    #------ Ut ------#
    Ut <- Yt
    # Evitar transformacao dos valores
    Yt[YtNA] <- -1
    
    #------ Vetor de PGs latentes correntes ------#
    # Xt <- BC(Yt = Yt, lambda = PSI.chain$lambda[i], csrd = csrd)
    Xt <- BC(Yt = Yt, lambda = 1, csrd = csrd)
    Xt[YtNA] <- Ut[YtNA]
    Yt[YtNA] <- Ut[YtNA]
    
    #------ Barra de progresso na tela ------#
    if((i%%update.step)==0){setTxtProgressBar(pb, i)}
  }
  
  #------ Fechando conexao da barra de progresso ------#
  close(pb)
  
  #------ Tempo total de processamento ------#
  b <- Sys.time()-a
  print(b)
  
  #-----------------------------------#
  #------ Tracos - Cadeias MCMC ------#
  #-----------------------------------#
  
  par(mar=c(4,4,.5,.5))
  par(mfrow=c(2,2))
  
  # #------ phi ------#
  # plot(sampMH$samples[,1],type="l",ylab=expression(phi),xlab="Iteracoes")
  # abline(v=burnin*salto,lty=2)
  # acf(sampMH$samples[,1],lag.max = 25,main="",ylab="",xlab="Defasagem")
  
  # #------ lambda ------#
  # plot(sampMH$samples[,2],type="l",ylab=expression(lambda),xlab="Iteracoes")
  # abline(v=burnin*salto,lty=2)
  # acf(sampMH$samples[,2],lag.max = 25,main="",ylab="",xlab="Defasagem")
  #------ lambda ------#
  # plot(sampMH$samples[,1],type="l",ylab=expression(lambda),xlab="Iteracoes")
  # abline(v=burnin*salto,lty=2)
  # acf(sampMH$samples[,1],lag.max = 25,main="",ylab="",xlab="Defasagem")
  
  #-------------------------------------------------#
  #------ Amostras da distribuicao posteriori ------#
  #-------------------------------------------------#
  
  #------ Thetat ------#
  par(mar=c(4,4,.5,.5))
  par(mfrow=c(2,2))
  graf.Thetat <- array(NA,dim=c(3,r,T+1))
  for(l in 1:(T+1))
  {
    for (c in (1:r))
    {
      graf.Thetat[,c,l] <- c(as.numeric(quantile(Thetat.chain.final[l,c,-(1:burnin)],.975,na.rm=T)),
                             as.numeric(quantile(Thetat.chain.final[l,c,-(1:burnin)],.5,na.rm=T)),
                             as.numeric(quantile(Thetat.chain.final[l,c,-(1:burnin)],.025,na.rm=T)))
    }
  }
  
  for (m in (1:(r-1)))
  {
    plot(1,type="n",
         ylim=c(min(graf.Thetat[3,m,-1],na.rm=T),max(graf.Thetat[1,m,-1],na.rm=T)),
         xlim=c(1,T),ylab=bquote(theta[.(m-1)]),xlab="Tempo")
    polygon(c(rev(1:T),1:T),c(rev(graf.Thetat[3,m,-1]),graf.Thetat[1,m,-1]),col="grey80",border=NA)
    lines(graf.Thetat[2,m,-1],lty=2,lwd=2)
    abline(h=0,lty=3)
  }
  
  #------ sigma2t ------#
  graf.sigma2t <- array(NA,dim=c(3,T+1))
  for(l in 1:(T+1))
  {
    graf.sigma2t[,l] <- c(as.numeric(quantile(sigma2t.chain.final[l,-(1:burnin)],.975,na.rm=T)),
                          as.numeric(quantile(sigma2t.chain.final[l,-(1:burnin)],.5,na.rm=T)),
                          as.numeric(quantile(sigma2t.chain.final[l,-(1:burnin)],.025,na.rm=T)))
  }
  plot(1,type="n",
       ylim=c(min(graf.sigma2t[3,-1],na.rm=T),max(graf.sigma2t[1,-1],na.rm=T)),
       xlim=c(1,T),ylab=expression(sigma^2),xlab="Tempo")
  polygon(c(rev(1:T),1:T),c(rev(graf.sigma2t[3,-1]),graf.sigma2t[1,-1]),col="grey80",border=NA)
  lines(graf.sigma2t[2,-1],lty=2,lwd=2)
  abline(h=0,lty=3)
  
  #------ phi ------#
  par(mar=c(4,4,.5,.5))
  par(mfrow=c(2,3))
  # plot(PSI.chain$phi[-(1:burnin)],type="l",ylab=expression(phi),xlab="Iteracoes")
  # acf(PSI.chain$phi[-(1:burnin)],lag.max = 20,main="",ylab="",xlab="Defasagem")
  # d <- density(x = PSI.chain$phi[-(1:burnin)], adjust = 2)
  # plot(d, xlab=expression(phi), main="",ylab="Densidade")
  # h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
  # h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
  # polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
  # lines(d,lwd=2)
  # lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
  # abline(h=0)
  
  #------ lambda ------#
  # plot(PSI.chain$lambda[-(1:burnin)],type="l",ylab=expression(lambda),xlab="Iteracoes")
  # acf(PSI.chain$lambda[-(1:burnin)],lag.max = 20,main="",ylab="",xlab="Defasagem")
  # d <- density(x = PSI.chain$lambda[-(1:burnin)], adjust = 2)
  # plot(d, xlab=expression(lambda), main="",ylab="Densidade")
  # h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
  # h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
  # polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
  # lines(d,lwd=2)
  # lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
  # abline(h=0)
  
  #------ Salvando resultado em forma de lista ------#
  resultado <- list(
    # taxa.aceit = round(sampMH$acceptance.rate,4),
    ht.post = ht.final,
    Ht.post = Ht.final,
    Thetat.post = Thetat.chain.final,
    sigma2t.post = sigma2t.chain.final,
    # phi.post = PSI.chain$phi,
    # lambda.post = PSI.chain$lambda,
    M = M,
    burnin = burnin,
    salto = salto,
    Desc = Desc,
    delta.sig2 = delta.sig2,
    T = T,
    r = r,
    Gt = Gt)
  
  #------ Retornando lista ------#
  return(resultado)
}

#################################
# Com os valores obtidos na fc anterior e feita a previsao
#################################

previsao.posteriori <- function(Ft.k,
                                D.full,
                                estimacao.resultados,
                                forecast.period,
                                Yt.k)
{
  #------ Contador de tempo ------#
  a <- Sys.time()
  
  #------ Transformacao Box-Cox inversa ------#
  BCinv <- function(X_BCinv,lambda_BCinv){(1+(lambda_BCinv*X_BCinv))^(1/lambda_BCinv)}
  
  #------ Censura ------#
  csrd <- 0
  
  #------ Total de localizacoes - Completo ------#
  n.new <-  dim(D.full)[1]
  
  #------ Previsoes a frente ------#
  T.new <- dim(Ft.k)[3]-1
  
  #------ Qtde. de amostras posterioris ------#
  M.post <- estimacao.resultados$M-estimacao.resultados$burnin
  
  #------ Array com prev a posteriori ------#
  Yt.prev <- array(NA,dim=c(n.new,M.post,T.new+1))
  
  #------ Matriz Gt ------#
  Gt <- estimacao.resultados$Gt
  
  #------ Matriz de descontos ------#
  Desc <- estimacao.resultados$Desc
  
  #------ Removendo o burn-in ------#
  # PSI.chain <- data.frame(phi = estimacao.resultados$phi.post,
  #                         lambda = estimacao.resultados$lambda.post)[-(1:estimacao.resultados$burnin),]
  # PSI.chain <- data.frame(lambda = data.frame(lambda = estimacao.resultados$lambda.post)[-(1:estimacao.resultados$burnin),])
  ht.final <- estimacao.resultados$ht.post[-(1:estimacao.resultados$burnin),]
  Ht.final <- estimacao.resultados$Ht.post[,,-(1:estimacao.resultados$burnin)]
  sigma2t.chain.final <- estimacao.resultados$sigma2t.post[,-(1:estimacao.resultados$burnin)]
  
  #------ Atualizando desconto da variancia ------#
  nt.k <- 1
  for (k in 2:(estimacao.resultados$T+1+T.new)){nt.k[k] <- estimacao.resultados$delta.sig2*nt.k[k-1] + 1}
  nt.k <- nt.k[(estimacao.resultados$T+1):(estimacao.resultados$T+1+T.new)]
  
  #------ Atualizacao da barra de progresso ------#
  update.step <- max(5, floor(M.post/100))
  
  #------ Barra de progresso ------#
  pb <- txtProgressBar(max = M.post, style = 3)
  
  for (m in 1:M.post)
  {
    ft.k <- array(NA,dim=c(n.new,1,T.new+1))
    at.k <- array(NA,dim=c(estimacao.resultados$r,1,T.new+1))
    at.k[,,1] <- ht.final[m,]
    Qt.k <- array(NA,dim=c(n.new,n.new,T.new+1))
    Rt.k <- array(NA,dim=c(estimacao.resultados$r,estimacao.resultados$r,T.new+1))
    Rt.k[,,1] <- Ht.final[,,m]
    # mat.cor.new <- MatCor(phi = PSI.chain$phi[m], DistMatrix = D.full)
    mat.cor.new <- MatCor(phi = 100000, DistMatrix = D.full)
    sigma2t.k <- sigma2t.chain.final[(T+1),m]
    for (k in 2:(T.new+1))
    {
      at.k[,,k] <- Gt%*%at.k[,,k-1]
      Rt.k[,,k] <- symmetricMatrix(Desc%*%Gt%*%Rt.k[,,k-1]%*%t(Gt)%*%Desc)
      ft.k <- t(Ft.k[,,k])%*%at.k[,,k]
      Qt.k <- symmetricMatrix(t(Ft.k[,,k])%*%Rt.k[,,k]%*%Ft.k[,,k]+
                                sigma2t.k*mat.cor.new)
      Xt.k <- as.numeric(mvnfast::rmvn(n = 1, mu = ft.k, sigma = Qt.k,isChol = TRUE))
      # Yt.prev[,m,k] <- ifelse(Xt.k>=((csrd^PSI.chain$lambda[m])-1)/PSI.chain$lambda[m],
      #                         BCinv(Xt.k,lambda = PSI.chain$lambda[m]),0)
      Yt.prev[,m,k] <- ifelse(Xt.k>=((csrd^1)-1)/1,
                              BCinv(Xt.k,lambda = 1),0)
    }
    
    #------ Barra de progresso na tela ------#
    if((m%%update.step)==0){setTxtProgressBar(pb, m)}
  }
  
  #------ Fechando conexao da barra de progresso ------#
  close(pb)
  
  #------ Qtd de localizacoes ------#
  n <- dim(Yt.k)[1]
  
  #------ Criterios de comparacao ------#
  #------ Raiz do erro quadratico medio ------#
  rEQM <- data.frame(forecast.period,
                     ensemble = NA,
                     post.mean = NA,
                     ind.mean = NA,
                     post.median = NA,
                     ind.median = NA)
  
  for (k in 2:(T.new+1))
  {
    # Media do ensemble
    rEQM[k-1,2] <- sqrt(mean((Yt.k[,,k] - Ft.k[1,1:n,k])^2,na.rm = T))
    # Media a posteriori
    rEQM[k-1,3] <- sqrt(mean((Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean))^2,na.rm = T))
    # Mediana a posteriori
    rEQM[k-1,5] <- sqrt(mean((Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k]))^2,na.rm = T))
  }
  
  rEQM[,4] <- ifelse(rEQM[,3]>rEQM[,2],1,0)
  rEQM[,6] <- ifelse(rEQM[,5]>rEQM[,2],1,0)
  
  #------ Erro absoluto medio ------#
  EAM <- data.frame(forecast.period,
                    ensemble = NA,
                    post.mean = NA,
                    ind.mean = NA,
                    post.median = NA,
                    ind.median = NA)
  
  for (k in 2:(T.new+1))
  {
    # Media do ensemble
    EAM[k-1,2] <- mean(abs(Yt.k[,,k] - Ft.k[1,1:n,k]),na.rm = T)
    # Media a posteriori
    EAM[k-1,3] <- mean(abs(Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean)),na.rm = T)
    # Mediana a posteriori
    EAM[k-1,5] <- mean(abs(Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k])),na.rm = T)
  }
  
  EAM[,4] <- ifelse(EAM[,3]>EAM[,2],1,0)
  EAM[,6] <- ifelse(EAM[,5]>EAM[,2],1,0)
  
  #------ Indice de concordancia de Willmott ------#
  ind.concord <- function(sim,obs)
  {1 - ((sum((obs-sim)^2,na.rm=T))/sum((abs(sim-mean(obs,na.rm=T)) + abs(obs-mean(obs,na.rm=T)))^2,na.rm=T))}
  
  ICW <- data.frame(forecast.period,
                    ensemble = NA,
                    post.mean = NA,
                    ind.mean = NA,
                    post.median = NA,
                    ind.median = NA)
  
  for (k in 2:(T.new+1))
  {
    # Media do ensemble
    ICW[k-1,2] <- ind.concord(sim=Ft.k[1,1:n,k], obs = Yt.k[,,k])
    # Media a posteriori
    ICW[k-1,3] <- ind.concord(sim=apply(Yt.prev[1:n,,k],1,mean), obs = Yt.k[,,k])
    # Mediana a posteriori
    ICW[k-1,5] <- ind.concord(sim=rowMedians(Yt.prev[1:n,,k]), obs = Yt.k[,,k])
  }
  
  ICW[,4] <- ifelse(ICW[,3]<ICW[,2],1,0)
  ICW[,6] <- ifelse(ICW[,5]<ICW[,2],1,0)
  
  #------ Graficos - Concordancia ------#
  par(mar=c(4,4,.5,.5))
  par(mfrow=c(1,2))
  plot(ceiling(min(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)):ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)),
       ceiling(min(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)):ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)),t="l",
       xlab="Previsão numérica", ylab="Observado")
  points(as.numeric(Ft.k[1,1:n,2:(T.new+1)]),
         as.numeric(Yt.k[,,2:(T.new+1)]),pch=19,col=2)
  abline(v=0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)),
         h=0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)),
         lty=2)
  
  plot(ceiling(min(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)):ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)),
       ceiling(min(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)):ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)),t="l",
       xlab="Previsão numérica calibrada", ylab="Observado")
  for (k in 2:(T.new+1))
  {
    points(rowMedians(Yt.prev[1:n,,k]),
           as.numeric(Yt.k[,,k]),pch=19,col=2)
  }
  abline(v=0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)),
         h=0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]),na.rm=T)),
         lty=2)
  
  #------ Resultado ------#
  resultado <- list(Yt.prev = Yt.prev,
                    rEQM = rEQM,
                    EAM = EAM,
                    ICW = ICW)
  
  #------ Tempo total de processamento ------#
  b <- Sys.time()-a
  print(b)
  
  return(resultado)
}

#-------------------------------------- FUNCOES - Rcpp ---------------------------------------

Rcpp::sourceCpp("C:/Users/b207056565/Desktop/Calibration_EtaModel/Funcs_cpp/Funcs_DGOP.cpp")

#----------------------------------------------------------------------------------------------
# RESULTADOS
#----------------------------------------------------------------------------------------------

DEFINICOES <- definicoes.usuario(forecast.date = as.Date("2017-08-21"),
                                 forecast.hour = 06,
                                 horizon = 96,
                                 lag.hours = 480)    

RESULTADO <- MCMC(M = 600,
                  burnin = 100,
                  salto = 20,
                  Desc = c(.99,
                           rep(.99,2)),
                  delta.sig2 = .95,
                  Yt = DEFINICOES$Yt,
                  Ft = DEFINICOES$Ft,
                  Gt = DEFINICOES$Gt,
                  Dist = DEFINICOES$D)

PREVISAO <- previsao.posteriori(Ft.k = DEFINICOES$Ft.k,
                                D.full = DEFINICOES$D.full,
                                estimacao.resultados = RESULTADO,
                                forecast.period = DEFINICOES$forecast.period,
                                Yt.k = DEFINICOES$Yt.k)

#----------------------------------------------------------------------------------------------
# RELATORIOS
#----------------------------------------------------------------------------------------------

# save.image(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m02/DGOP_teste.Rda"))

View(PREVISAO$rEQM)
View(PREVISAO$EAM)
View(PREVISAO$ICW)

#------ Graficos ------#
Yt <- DEFINICOES$Yt
T.new <- length(DEFINICOES$forecast.period)
Yt.prev <- PREVISAO$Yt.prev
Yt.k <- DEFINICOES$Yt.k
Ft.new <- DEFINICOES$Ft.k
n <- dim(DEFINICOES$Yt)[1]

par(mfrow=c(2,1))
par(mar=c(4,4,3,.5))
for (l in 1:n)
{
  plot(0,type="n",ylim=c(5,15),xlim=c(1,(T.new+1)),
       ylab=DEFINICOES$relatorio$codigo[l],xlab="horizonte", axes= F)
  axis(1,1:(length(DEFINICOES$forecast.period)+1),c("0",paste0("+",1:length(DEFINICOES$forecast.period))))
  axis(3,1:(length(DEFINICOES$forecast.period)+1),
       c("5.5h",paste0(rep(c(seq(6,23.5,.5),seq(0,5.5,.5)),
                           length(DEFINICOES$forecast.period)/48),"h"),
         "6h"))
  abline(h=0)
  polygon(c(rev(2:(T.new+1)),2:(T.new+1)),
          c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.025)),
            colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.975)),
          col="grey80",border=NA)
  polygon(c(rev(2:(T.new+1)),2:(T.new+1)),
          c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
            colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
          col="grey90",border=NA)
  abline(v=1:(length(DEFINICOES$forecast.period)+1),lty=2,col="grey80")
  axis(2,5:15,5:15)
  points(1,Yt[l,,T+1],pch=19,col=2)
  points(2:(T.new+1),Yt.k[l,,2:(T.new+1)],pch=19,col=2)
  lines(2:(T.new+1),colMeans(Yt.prev[l,,2:(T.new+1)]),lty=2,lwd=2)
  lines(2:(T.new+1),Ft.new[1,l,2:(T.new+1)],lty=2,lwd=2,col=4)
}

rm(Yt,T.new,Yt.prev,Yt.k,Ft.new,l,n)

par(mfrow=c(1,1))
par(mar=c(4,4,.5,.5))
plot(PREVISAO$rEQM[,c(2)],t="l",lwd=2,xlab="tempo",ylab="rEQM",
     ylim = c(min(PREVISAO$rEQM[,c(2,5)]),max(PREVISAO$rEQM[,c(2,5)])))
lines(PREVISAO$rEQM[,c(5)],lwd=2,col=2)

par(mfrow=c(1,1))
par(mar=c(4,4,.5,.5))
plot(PREVISAO$EAM[,c(2)],t="l",lwd=2,xlab="tempo",ylab="EAM",
     ylim = c(min(PREVISAO$EAM[,c(2,5)]),max(PREVISAO$EAM[,c(2,5)])))
lines(PREVISAO$EAM[,c(5)],lwd=2,col=2)

par(mfrow=c(1,1))
par(mar=c(4,4,.5,.5))
plot(PREVISAO$ICW[,c(2)],t="l",lwd=2,xlab="tempo",ylab="ICW",
     ylim = c(min(PREVISAO$ICW[,c(2,5)]),max(PREVISAO$ICW[,c(2,5)])))
lines(PREVISAO$ICW[,c(5)],lwd=2,col=2)

data.frame(rEQM = colMeans2(as.matrix(PREVISAO$rEQM[,c(2,5)])),
           EAM = colMeans2(as.matrix(PREVISAO$EAM[,c(2,5)])),
           ICW = colMeans2(as.matrix(PREVISAO$ICW[,c(2,5)])))

data.frame(rEQM = round(sum(as.matrix(PREVISAO$rEQM[,c(6)]))/nrow(PREVISAO$rEQM),2),
           EAM = round(sum(as.matrix(PREVISAO$EAM[,c(6)]))/nrow(PREVISAO$EAM),2),
           ICW = round(sum(as.matrix(PREVISAO$ICW[,c(6)]))/nrow(PREVISAO$ICW),2))
