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
##########|             SPATIOTEMPORAL ENSEMBLE MOS               |###########
##########|-------------------------------------------------------|###########
##########|-------------------------------------------------------|###########
##############################################################################
##############################################################################

#------------------------------------ CARREGANDO PACOTES --------------------------------------
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

#------------------------------------ DATAS PARA RODAR ----------------------------------------
datas <- seq(as.Date("2015-12-21"),
             as.Date("2016-11-21"),
             by="1 month")

for (data.c in 1:12)
{
  print(paste0("RODANDO PARA O MES ",data.c))
  
  #---------------------------------------- LENDO DADOS ------------------------------------------
  coords.ETA <- read.csv("C:/Users/b207056565/Desktop/Calibration_EtaModel/coords_Eta15km_MG.csv")
  coords.MG <- read.csv2("C:/Users/b207056565/Desktop/Calibration_EtaModel/coordMG.csv")
  coords.stations <- read.csv2("C:/Users/b207056565/Desktop/Calibration_EtaModel/Dados_Observados/info104_inmet_MG_e_limites_ALTERADO59.csv")
  gradeMG <- read.csv2("C:/Users/b207056565/Desktop/Calibration_EtaModel/gradeMG.csv")
  obs.wind <- read.csv2("C:/Users/b207056565/Desktop/Calibration_EtaModel/Dados_Observados/velvento104_imput_20151001_20170930 - ALTERADO59.csv")
  
  #----------------------------------------------------------------------------------------------
  # OPCIONAIS DO USUARIO - Definicoes
  #----------------------------------------------------------------------------------------------
  
  #------ Data para previsao ------#
  forecast.date <- datas[data.c]
  forecast.hour <- 12
  
  #------ Horizonte ------#
  horizon <- 96
  
  #------ Defasagem utilizada ------#
  lag.days <- 90
  
  #------ Periodo de treinamento ------#
  train.period <- seq(as.Date(forecast.date-2*(horizon/24)-lag.days+1),
                      as.Date(forecast.date-2*(horizon/24)),by = '1 day')
  
  #------ Recuperando as informacoes nos arquivos .csv ------#
  #----------------------------------------------------------------------------------------------
  # VELOCIDADE DO VENTO - 10 m
  #----------------------------------------------------------------------------------------------
  wind10m.ETA <- data.frame(NULL)
  aux.date <- unique(c(train.period-15,train.period))
  
  #------ Removendo dias que nao houve funcionamento do Eta ------#
  aux.date.rm <- NULL
  if (sum(aux.date < as.Date("2015-10-11"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date < as.Date("2015-10-11")))
  if (sum(aux.date > as.Date("2017-09-28"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date > as.Date("2017-09-28")))
  if (sum(aux.date == as.Date("2015-11-09"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-11-09")))
  if (sum(aux.date == as.Date("2015-12-08"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-12-08")))
  if (sum(aux.date == as.Date("2015-12-25"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-12-25")))
  if (sum(aux.date == as.Date("2016-01-10"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-01-10")))
  if (sum(aux.date == as.Date("2016-03-21"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-03-21")))
  if (sum(aux.date == as.Date("2016-08-07"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-08-07")))
  if (sum(aux.date == as.Date("2017-01-02"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-01-02")))
  if (sum(aux.date == as.Date("2017-02-05"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-02-05")))
  if (sum(aux.date == as.Date("2017-03-28"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-03-28")))
  if (sum(aux.date == as.Date("2017-09-13"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-09-13")))
  if (length(aux.date.rm)>0) aux.date <- aux.date[-aux.date.rm]

  for (i in 1:length(aux.date))
  {
    #------ lendo .csv auxiliar ------#
    aux.csv <- read.csv(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/velvento10m/",
                               "velvento10m_Eta15km_MG_",
                               substr(aux.date[i],1,4),
                               substr(aux.date[i],6,7),
                               substr(aux.date[i],9,10),
                               "12.csv"))
    
    #------ organizando dia/hora origem ------#
    aux.csv[,1] <- as.character(strptime(paste0(substr(aux.csv[,1],1,4),"-",
                                                substr(aux.csv[,1],5,6),"-",
                                                substr(aux.csv[,1],7,8)," ",
                                                substr(aux.csv[,1],9,10),":00:01"),
                                         format = "%Y-%m-%d %H:%M:%S",
                                         tz = "GMT"))
    
    #------ organizando dia/hora data de previsao ------#
    aux.csv[,2] <- as.character(strptime(paste0(substr(aux.csv[,2],1,4),"-",
                                                substr(aux.csv[,2],5,6),"-",
                                                substr(aux.csv[,2],7,8)," ",
                                                substr(aux.csv[,2],9,10),":00:01"),
                                         format = "%Y-%m-%d %H:%M:%S",
                                         tz = "GMT"))
    
    #------ criando a variavel "horizonte" ------#
    aux.csv <- data.frame(aux.csv[,1:2],
                          horizonte = as.numeric(difftime(time1 = strptime(aux.csv[,2], format = "%Y-%m-%d %H:%M:%S",tz = "GMT"),
                                                          time2 = strptime(aux.csv[,1], format = "%Y-%m-%d %H:%M:%S",tz = "GMT"),
                                                          tz = "GMT", units = "hours")),
                          aux.csv[,-(1:2)])
    
    #------ mantendo apenas o horario desejado e o horizonte de interesse ------#
    aux.csv <- subset(aux.csv,substr(aux.csv$data,12,13)==forecast.hour & aux.csv$horizonte >= horizon)
    
    #------ interpolando para as estacoes (interpolacao bilinear) ------#
    for (l in 1:dim(aux.csv)[1])
    {
      wind10m.ETA <- rbind(wind10m.ETA,data.frame(t(c(as.character(aux.csv[l,1:3]),
                                                      suppressWarnings(interpp(x = coords.ETA$longitude,
                                                                               y = coords.ETA$latitude,
                                                                               z = as.numeric(aux.csv[l,-(1:3)]),
                                                                               xo = coords.stations$longitude,
                                                                               yo = coords.stations$latitude)$z)))))
    }
  }
  
  #------ formatando as variaveis ------#
  wind10m.ETA[,1] <- as.character(wind10m.ETA[,1])
  wind10m.ETA[,2] <- as.character(wind10m.ETA[,2])
  for (i in 3:dim(wind10m.ETA)[2]){wind10m.ETA[,i] <- as.numeric(as.character(wind10m.ETA[,i]))}
  
  #------ nomeando as variaveis ------#
  names(wind10m.ETA) <- c(names(aux.csv)[1:3],as.character(coords.stations$codigo))
  
  #-------------------------------------- FUNCOES - Rcpp ---------------------------------------
  
  Rcpp::sourceCpp("C:/Users/b207056565/Desktop/Calibration_EtaModel/Funcs_STEMOS.cpp")
  
  #------ Transformacao Box-Cox inversa ------#
  BCinv <- function(X_BCinv,lambda_BCinv){(1+(lambda_BCinv*X_BCinv))^(1/lambda_BCinv)}
  
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
  
  p.log <- function(x)
  {
    if (all(x>0))
    {
      loglikep(Yt = Yt,
               Ft = Ft,
               Thetat = Thetat.chain.final[,,i],
               Vt = genVt(b0 = x[1],
                          b1 = x[2],
                          s2t = s2t,
                          phi = x[3],
                          DistMatrix = D),
               lambda = x[4],
               csrd = csrd) +
        log.b0.prior(x[1]) + log.b1.prior(x[2]) +
        log.phi.prior(x[3]) + log.lambda.prior(x[4])
    } else {-9.99e+300}
  }
  
  ind.concord <- function(sim,obs)
  {1 - ((sum((obs-sim)^2))/sum((abs(sim-mean(obs)) + abs(obs-mean(obs)))^2))}
  
  #----------------------------------------------------------------------------------------------
  # DADOS
  #----------------------------------------------------------------------------------------------
  
  #------ Matriz de distancias D ------#
  D <- as.matrix(dist(data.frame(coords.stations$longitude,coords.stations$latitude)))
  
  #------ Matriz de evolucao Gt ------#
  Gt <- as.matrix(bdiag(diag(1,6)))
  
  #------ Comprimento do vetor Theta (Tendencia espacial + Sazonalidade) ------#
  r <- dim(Gt)[1]
  
  train.period <- train.period + (horizon/24)
  train.period.rm <- NULL
  
  if (sum(train.period < as.Date("2015-10-16"))>0) train.period.rm <- c( train.period.rm,which( train.period < as.Date("2015-10-16")))
  if (sum(train.period > as.Date("2017-09-28"))>0) train.period.rm <- c(train.period.rm,which(train.period > as.Date("2017-09-28")))
  if (length(train.period.rm)>0) train.period <- train.period[-train.period.rm]
  
  #----------------------------------------------------------------------------------------------
  # ALOCANDO OS VETORES
  #----------------------------------------------------------------------------------------------
  
  #------ Total de tempos ------#
  T <- length(train.period)
  
  #------ Total de localizacoes ------#
  n <-  dim(coords.stations)[1]
  
  #------ valor da censura ------#
  csrd <- 0.1
  
  #------ Matriz de covariaveis Ft ------#
  #------ Valor observado ------#
  Yt <- array(NA,dim=c(n,1,T+1))
  Ft <- array(NA,dim=c(r,n,T+1))
  s2t <- array(NA,dim=c(1,n,T+1))
  for (t in 2:(T+1))
  {
    Yt[,,t] <- as.numeric(subset(obs.wind, substr(obs.wind$data.hora,12,13) == forecast.hour &
                                   substr(obs.wind$data.hora,1,10) == as.Date(train.period[t-1]))[,-1])
    
    Ft[,,t] <- t(data.frame(level = 1,
                            ensemble.mean = apply(subset(wind10m.ETA, 
                                                         wind10m.ETA$data == train.period[t-1])[,-(1:3)],2,mean),
                            AR = as.numeric(subset(obs.wind,obs.wind$data.hora ==
                                                     as.character(strptime(paste0(train.period[t-1]," ",forecast.hour,":00:01"),
                                                                           format = "%Y-%m-%d %H:%M:%S",
                                                                           tz = "GMT")-(horizon/24)*24*60*60))[,-1]),
                            alt = coords.stations$altitude,
                            lat = coords.stations$latitude,
                            long = coords.stations$longitude))
    
    s2t[,,t] <- apply(subset(wind10m.ETA, 
                             wind10m.ETA$data == train.period[t-1])[,-(1:3)],2,var)
  }
  
  #----------------------------------------------------------------------------------------------
  # MCMC
  #----------------------------------------------------------------------------------------------
  
  #------ Num de iteracoes M ------#
  M <- 2100
  
  #------ Burn-in ------#
  burnin <- 100
  
  #------ Salto/Afinamento (thinning) ------#
  salto <- 30
  
  #------ Proporcao de obs censuradas (ideal < 15%) ------#
  print(paste0(round(length(which(Yt[,,-1] < csrd))/length(as.numeric(Yt[,,-1]))*100,2),"%"))
  
  #------------------------------------#
  #------ Distribuicoes a priori ------#
  #------------------------------------#
  
  #------ b0 ------#
  mu.b0.prior <- 0
  sigma.b0.prior <- 1
  log.b0.prior <- function(x){
    if(length(x)>1)
    {
      result <- NULL
      for (i in 1:length(x))
      {
        result <- c(result,tmvtnorm::dtmvnorm(x[i],
                           mean = mu.b0.prior,
                           sigma = sigma.b0.prior,
                           lower = 0,log = TRUE))
      }
      return(result)
    } else {tmvtnorm::dtmvnorm(x,
                               mean = mu.b0.prior,
                               sigma = sigma.b0.prior,
                               lower = 0,log = TRUE)}
      
    }
  
  #------ b1 ------#
  mu.b1.prior <- 0
  sigma.b1.prior <- 1
  log.b1.prior <- function(x){
    if(length(x)>1)
    {
      result <- NULL
      for (i in 1:length(x))
      {
        result <- c(result,tmvtnorm::dtmvnorm(x[i],
                                              mean = mu.b1.prior,
                                              sigma = sigma.b1.prior,
                                              lower = 0,log = TRUE))
      }
      return(result)
    } else {tmvtnorm::dtmvnorm(x,
                               mean = mu.b1.prior,
                               sigma = sigma.b1.prior,
                               lower = 0,log = TRUE)}
    
  }
  
  #------ phi ------#
  a.phi.prior <- 2
  b.phi.prior <- max(D)/6
  log.phi.prior <- function(x){dgamma(x,shape = a.phi.prior,rate = b.phi.prior,log = T)}
  
  #------ lambda ------#
  mu.lambda.prior <- 1
  sigma.lambda.prior <- 1
  log.lambda.prior <- function(x){dnorm(x,mean = mu.lambda.prior,sd = sqrt(sigma.lambda.prior),log = T)}
  
  #------------------------------------------#
  #------ Valores iniciais das cadeias ------#
  #------------------------------------------#
  
  #------ Thetat final ------#
  Thetat.chain.final <- array(NA, dim = c(T+1,r,M))
  
  #------ b0, b1, phi e lambda = PSI ------#
  PSI.chain <- data.frame(t(c(3,2,4,2)))
  names(PSI.chain) <- c("b0","b1","phi","lambda")
  
  #------ Matriz de correlacao da proposta M-H ------#
  EPS.PSI <- diag(2.4/sqrt(dim(PSI.chain)[2]),dim(PSI.chain)[2])
  
  #------ Vetor de PGs latentes correntes ------#
  Xt <- BC(Yt = Yt, lambda = PSI.chain$lambda, csrd = csrd)
  
  #------------------#
  #------ MCMC ------#
  #------------------#
  
  #------ Matriz de correlacao corrente ------#
  Vt.c <- genVt(b0 = PSI.chain$b0,
                b1 = PSI.chain$b1,
                s2t = s2t,
                phi = PSI.chain$phi,
                DistMatrix = D)
  
  #------ Contador de tempo ------#
  a <- Sys.time()
  
  #------ Atualizacao da barra de progresso ------#
  update.step <- max(5, floor(M/100))
  
  #------ Barra de progresso ------#
  pb <- txtProgressBar(max = M, style = 3)
  
  for(i in 2:M)
  {
    #-----------------------------------------------#
    #------ Amostrador de Gibbs p/ para Theta ------#
    #-----------------------------------------------#
    
    mu.Theta.posterior <- matrix(0,nrow = r)
    Sigma.Theta.posterior <- matrix(0,nrow = r,ncol = r)
    for (t in 2:(T+1))
    {
      Sigma.Theta.posterior <- Sigma.Theta.posterior + Ft[,,t]%*%solve(Vt.c[,,t])%*%t(Ft[,,t])
      mu.Theta.posterior <- mu.Theta.posterior + Ft[,,t]%*%solve(Vt.c[,,t])%*%Xt[,,t]
    }
    Sigma.Theta.posterior <- symmetricMatrix(solve(Sigma.Theta.posterior + diag(r)))
    mu.Theta.posterior <- Sigma.Theta.posterior%*%mu.Theta.posterior
    Thetat.chain.final[,,i] <- matrix(rmvn(n = 1,mu = mu.Theta.posterior,sigma = Sigma.Theta.posterior,isChol = TRUE),
                                      nrow = T+1,ncol=r, byrow = TRUE)
    
    Xt.mu <- array(NA,dim=c(n,1,T+1))
    for (t in 2:(T+1))
    {
      Xt.mu[,,t] <- t(Ft[,,t])%*%Thetat.chain.final[1,,i]
    }
    
    #--------------------------------------------------------------------#
    #------ Metropolis-Hastings com propostas adaptativas para PSI ------#
    #--------------------------------------------------------------------#
    
    if (i==2)
    {
      sampMH <- MHadapt(p.log, n=salto, init=as.numeric(PSI.chain[i-1,]),
                        scale = EPS.PSI,
                        adapt=TRUE, acc.rate=0.234, showProgressBar=FALSE)
    } else {
      sampMH <- MHadapt.add(sampMH, n.update=salto,showProgressBar=FALSE)
    }
    
    #------ Thinning ------# 
    PSI.chain[i,] <- sampMH$samples[(i-1)*salto,]
    
    #------ Atualizacao da matriz de correlacao corrente ------#
    if (all.equal(as.numeric(PSI.chain[i,]),as.numeric(PSI.chain[i-1,]))!=T)
    {
      Vt.c <- genVt(b0 = PSI.chain$b0[i],
                    b1 = PSI.chain$b1[i],
                    s2t = s2t,
                    phi = PSI.chain$phi[i],
                    DistMatrix = D)
    }
    
    #--------------------------------------------------#
    #------ Amostrador de Gibbs p/ Z var latente ------#
    #--------------------------------------------------#
    for (t in 2:(T+1))
    {
      if (length(which(Yt[,,t] < csrd))>0)
      {
        cens.c <- which(Yt[,,t] < csrd)
        restr <- rep(Inf,n)
        restr[cens.c] <- rep(((csrd^PSI.chain$lambda[i])-1)/PSI.chain$lambda[i],length(cens.c))
        Zt.c <- tmvtnorm::rtmvnorm(n=1,
                                   mean = Xt.mu[,,t],
                                   sigma= Vt.c[,,t],
                                   upper=restr,
                                   algorithm = "gibbs",
                                   burn.in.samples = 0)
        # start.value = Yt[,,t])
        Yt[cens.c,,t] <- ifelse(is.nan(Zt.c[,cens.c])==F,
                                Zt.c[,cens.c],
                                Yt[cens.c,,t])
      }
    }
    
    #------ Vetor de PGs latentes correntes ------#
    Xt <- BC(Yt = Yt, lambda = PSI.chain$lambda[i], csrd = csrd)
    
    #------ Barra de progresso na tela ------#
    if((i%%update.step)==0){setTxtProgressBar(pb, i)}
  }
  
  #------ Fechando conexao da barra de progresso ------#
  close(pb)
  
  #------ Tempo total de processamento ------#
  b <- Sys.time()-a
  print(b)
  
  #------ Taxa de aceitacao p/ M-H ------#
  print(paste0("Acceptance rate: ",round(sampMH$acceptance.rate,4)*100,"%"))
  
  #-----------------------------------#
  #------ Tracos - Cadeias MCMC ------#
  #-----------------------------------#
  
  par(mar=c(4,4,.5,.5))
  par(mfrow=c(4,2))
  
  #------ b0 ------#
  plot(sampMH$samples[,1],type="l",ylab=expression(beta[0]),xlab="Iteracoes")
  abline(v=burnin*salto,lty=2)
  acf(sampMH$samples[,1],lag.max = 25,main="",ylab="",xlab="Defasagem")
  
  #------ b1 ------#
  plot(sampMH$samples[,2],type="l",ylab=expression(beta[1]),xlab="Iteracoes")
  abline(v=burnin*salto,lty=2)
  acf(sampMH$samples[,2],lag.max = 25,main="",ylab="",xlab="Defasagem")
  
  #------ phi ------#
  plot(sampMH$samples[,3],type="l",ylab=expression(phi),xlab="Iteracoes")
  abline(v=burnin*salto,lty=2)
  acf(sampMH$samples[,3],lag.max = 25,main="",ylab="",xlab="Defasagem")
  
  #------ lambda ------#
  plot(sampMH$samples[,4],type="l",ylab=expression(lambda),xlab="Iteracoes")
  abline(v=burnin*salto,lty=2)
  acf(sampMH$samples[,4],lag.max = 25,main="",ylab="",xlab="Defasagem")
  
  #-------------------------------------------------#
  #------ Amostras da distribuicao posteriori ------#
  #-------------------------------------------------#
  
  #------ Thetat ------#
  par(mar=c(4,4,.5,.5))
  par(mfrow=c(4,2))
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
  
  for (m in (1:r))
  {
    plot(1,type="n",
         ylim=c(min(graf.Thetat[3,m,-1],na.rm=T),max(graf.Thetat[1,m,-1],na.rm=T)),
         xlim=c(1,T),ylab=bquote(theta[.(m-1)]),xlab="Tempo")
    polygon(c(rev(1:T),1:T),c(rev(graf.Thetat[3,m,-1]),graf.Thetat[1,m,-1]),col="grey80",border=NA)
    lines(graf.Thetat[2,m,-1],lty=2,lwd=2)
    abline(h=0,lty=3)
  }
  
  #------ b0 ------#
  par(mar=c(4,4,.5,.5))
  par(mfrow=c(4,3))
  plot(PSI.chain$b0[-(1:burnin)],type="l",ylab=expression(beta[0]),xlab="Iteracoes")
  acf(PSI.chain$b0[-(1:burnin)],lag.max = 20,main="",ylab="",xlab="Defasagem")
  d <- density(x = PSI.chain$b0[-(1:burnin)], adjust = 2)
  plot(d, xlab=expression(beta[0]), main="", ylab="Densidade")
  h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
  h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b0[-(1:burnin)],c(0.975))))]
  polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
  lines(d,lwd=2)
  lines(seq(min(d$x),max(d$x),.01),exp(log.b0.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
  abline(h=0)
  
  #------ b1 ------#
  plot(PSI.chain$b1[-(1:burnin)],type="l",ylab=expression(beta[1]),xlab="Iteracoes")
  acf(PSI.chain$b1[-(1:burnin)],lag.max = 20,main="",ylab="",xlab="Defasagem")
  d <- density(x = PSI.chain$b1[-(1:burnin)], adjust = 2,from = 0)
  plot(d, xlab=expression(beta[1]), main="",ylab="Densidade")
  h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
  h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$b1[-(1:burnin)],c(0.975))))]
  polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
  lines(d,lwd=2)
  lines(seq(min(d$x),max(d$x),.01),exp(log.b1.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
  abline(h=0)
  
  #------ phi ------#
  plot(PSI.chain$phi[-(1:burnin)],type="l",ylab=expression(phi),xlab="Iteracoes")
  acf(PSI.chain$phi[-(1:burnin)],lag.max = 20,main="",ylab="",xlab="Defasagem")
  d <- density(x = PSI.chain$phi[-(1:burnin)], adjust = 2)
  plot(d, xlab=expression(phi), main="",ylab="Densidade")
  h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
  h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
  polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
  lines(d,lwd=2)
  lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
  abline(h=0)
  
  #------ lambda ------#
  plot(PSI.chain$lambda[-(1:burnin)],type="l",ylab=expression(lambda),xlab="Iteracoes")
  acf(PSI.chain$lambda[-(1:burnin)],lag.max = 20,main="",ylab="",xlab="Defasagem")
  d <- density(x = PSI.chain$lambda[-(1:burnin)], adjust = 2)
  plot(d, xlab=expression(lambda), main="",ylab="Densidade")
  h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
  h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
  polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
  lines(d,lwd=2)
  lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
  abline(h=0)
  
  #----------------------------------------------------------------------------------------------
  # ENTRANDO COM DADOS PARA TODA REGIAO
  #----------------------------------------------------------------------------------------------
  
  #------ Periodo de previsao ------#
  forecast.period <- seq(tail(train.period, n=1)+1-(horizon/24),
                         forecast.date-(horizon/24),
                         by="1 day")
  
  #------ Recuperando as informacoes nos arquivos .csv ------#
  
  #----------------------------------------------------------------------------------------------
  # VELOCIDADE DO VENTO - 10 m
  #----------------------------------------------------------------------------------------------
  wind10m.ETA <- data.frame(NULL)
  aux.date <- seq(head(forecast.period, n=1)-15,tail(forecast.period, n=1), by = "1 day")
  
  #------ Removendo dias que nao houve funcionamento do Eta ------#
  aux.date.rm <- NULL
  if (sum(aux.date < as.Date("2015-10-11"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date < as.Date("2015-10-11")))
  if (sum(aux.date > as.Date("2017-09-28"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date > as.Date("2017-09-28")))
  if (sum(aux.date == as.Date("2015-11-09"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-11-09")))
  if (sum(aux.date == as.Date("2015-12-08"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-12-08")))
  if (sum(aux.date == as.Date("2015-12-25"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2015-12-25")))
  if (sum(aux.date == as.Date("2016-01-10"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-01-10")))
  if (sum(aux.date == as.Date("2016-03-21"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2016-03-21")))
  if (sum(aux.date == as.Date("2017-01-02"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-01-02")))
  if (sum(aux.date == as.Date("2017-02-05"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-02-05")))
  if (sum(aux.date == as.Date("2017-03-28"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-03-28")))
  if (sum(aux.date == as.Date("2017-09-13"))>0) aux.date.rm <- c(aux.date.rm,which(aux.date == as.Date("2017-09-13")))
  if (length(aux.date.rm)>0) aux.date <- aux.date[-aux.date.rm]

  for (i in 1:length(aux.date))
  {
    #------ lendo .csv auxiliar
    aux.csv <- read.csv(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/velvento10m/",
                               "velvento10m_Eta15km_MG_",
                               substr(aux.date[i],1,4),
                               substr(aux.date[i],6,7),
                               substr(aux.date[i],9,10),
                               "12.csv"))
    
    #------ organizando dia/hora origem ------#
    aux.csv[,1] <- as.character(strptime(paste0(substr(aux.csv[,1],1,4),"-",
                                                substr(aux.csv[,1],5,6),"-",
                                                substr(aux.csv[,1],7,8)," ",
                                                substr(aux.csv[,1],9,10),":00:01"),
                                         format = "%Y-%m-%d %H:%M:%S",
                                         tz = "GMT"))
    
    #------ organizando dia/hora data de previsao ------#
    aux.csv[,2] <- as.character(strptime(paste0(substr(aux.csv[,2],1,4),"-",
                                                substr(aux.csv[,2],5,6),"-",
                                                substr(aux.csv[,2],7,8)," ",
                                                substr(aux.csv[,2],9,10),":00:01"),
                                         format = "%Y-%m-%d %H:%M:%S",
                                         tz = "GMT"))
    
    #------ criando a variavel "horizonte" ------#
    aux.csv <- data.frame(aux.csv[,1:2],
                          horizonte = as.numeric(difftime(time1 = strptime(aux.csv[,2], format = "%Y-%m-%d %H:%M:%S",tz = "GMT"),
                                                          time2 = strptime(aux.csv[,1], format = "%Y-%m-%d %H:%M:%S",tz = "GMT"),
                                                          tz = "GMT", units = "hours")),
                          aux.csv[,-(1:2)])
    
    #------ mantendo apenas o horario desejado e o horizonte de interesse ------#
    aux.csv <- subset(aux.csv,substr(aux.csv$data,12,13)==forecast.hour & aux.csv$horizonte >= horizon)
    
    #------ interpolando para as estacoes (interpolacao bilinear) ------#
    for (l in 1:dim(aux.csv)[1])
    {
      wind10m.ETA <- rbind(wind10m.ETA,data.frame(t(c(as.character(aux.csv[l,1:3]),
                                                      suppressWarnings(interpp(x = coords.ETA$longitude,
                                                                               y = coords.ETA$latitude,
                                                                               z = as.numeric(aux.csv[l,-(1:3)]),
                                                                               xo = c(coords.stations$longitude,coords.ETA$longitude),
                                                                               yo = c(coords.stations$latitude,coords.ETA$latitude))$z)))))
    }
  }
  
  #------ formatando as variaveis ------#
  wind10m.ETA[,1] <- as.character(wind10m.ETA[,1])
  wind10m.ETA[,2] <- as.character(wind10m.ETA[,2])
  for (i in 3:dim(wind10m.ETA)[2]){wind10m.ETA[,i] <- as.numeric(as.character(wind10m.ETA[,i]))}
  
  #------ nomeando as variaveis ------#
  names(wind10m.ETA) <- c(names(aux.csv)[1:3],as.character(coords.stations$codigo),as.character(coords.ETA$ID))
  
  #------ Matriz de distancias D ------#
  D.new <- as.matrix(dist(data.frame(c(coords.stations$longitude,coords.ETA$longitude),
                                     c(coords.stations$latitude,coords.ETA$latitude))))
  D.new <- matrix(as.numeric(D.new),dim(D.new)[1],dim(D.new)[2])
  
  #------ Total de localizacoes - Completo ------#
  n.new <-  dim(D.new)[1]
  
  #------ Previsoes a frente ------#
  T.new <- horizon/24
  
  #------ Qtde. de amostras posterioris ------#
  M.post <- dim(PSI.chain)[1]-burnin
  
  #------ Removendo a defasagem na data ------#
  forecast.period <- forecast.period + (horizon/24)
  
  #------ Matriz de covariaveis Ft ------#
  #------ Valor observado ------#
  Yt.k <- array(NA,dim=c(n,1,T.new+1))
  Yt.prev <- array(NA,dim=c(n.new,M.post,T.new+1))
  Ft.new <- array(NA,dim=c(r,n.new,T.new+1))
  Ft.numforecast <- array(NA,dim=c(r,n.new,T.new+1))
  s2t.new <- array(NA,dim=c(1,n.new,T.new+1))
  for (t in 2:(T.new+1))
  {
    Yt.k[,,t] <- as.numeric(subset(obs.wind, substr(obs.wind$data.hora,12,13) == forecast.hour &
                                     substr(obs.wind$data.hora,1,10) == as.Date(forecast.period[t-1]))[,-1])
    
    Ft.new[,,t] <- t(data.frame(level = 1,
                                ensemble.mean = apply(subset(wind10m.ETA, 
                                                             wind10m.ETA$data == forecast.period[t-1])[,-(1:3)],
                                                      2,mean),
                                AR = c(as.numeric(subset(obs.wind,obs.wind$data.hora ==
                                                           as.character(strptime(paste0(forecast.period[t-1]," ",forecast.hour,":00:01"),
                                                                                 format = "%Y-%m-%d %H:%M:%S",
                                                                                 tz = "GMT")-(horizon/24)*24*60*60))[,-1]),rep(1,n.new-n)),
                                alt = c(coords.stations$altitude,rep(1,n.new-n)),
                                lat = c(coords.stations$latitude,rep(1,n.new-n)),
                                long = c(coords.stations$longitude,rep(1,n.new-n))))
    
    
    Ft.numforecast[,,t] <- t(data.frame(level = 1,
                                        ensemble.mean = as.numeric(subset(wind10m.ETA, 
                                                                          wind10m.ETA$data == forecast.period[t-1] &
                                                                            wind10m.ETA$horizonte ==  horizon)[,-(1:3)]),
                                        AR = c(as.numeric(subset(obs.wind,obs.wind$data.hora ==
                                                                   as.character(strptime(paste0(forecast.period[t-1]," ",forecast.hour,":00:01"),
                                                                                         format = "%Y-%m-%d %H:%M:%S",
                                                                                         tz = "GMT")-(horizon/24)*24*60*60))[,-1]),rep(1,n.new-n)),
                                        alt = c(coords.stations$altitude,rep(1,n.new-n)),
                                        lat = c(coords.stations$latitude,rep(1,n.new-n)),
                                        long = c(coords.stations$longitude,rep(1,n.new-n))))
    
    s2t.new[,,t] <- apply(subset(wind10m.ETA, 
                                 wind10m.ETA$data == forecast.period[t-1])[,-(1:3)],2,var)
  }
  
  #----------------------------------------------------------------------------------------------
  # RESULTADOS
  #----------------------------------------------------------------------------------------------
 
  #------ Removendo o burn-in ------#
  PSI.chain <- PSI.chain[-(1:burnin),]
  Thetat.chain.final <- Thetat.chain.final[,,-(1:burnin)]

  ##############################################
  ###### PREVENDO APENAS PARA AS ESTAÇÕES ######
  ##############################################
 
  parc.prev <- 1:n
  Yt.prev <- array(NA,dim=c(n,M.post,T.new+1))
  s2t.new2 <- s2t.new
  s2t.new <- array(NA,dim=c(1,n,T.new+1))
  for (l in parc.prev)
  {
    s2t.new[,l,] <- s2t.new2[1,l,]
  }
  rm(s2t.new2)
  
  #------ Contador de tempo ------#
  a1 <- Sys.time()
  
  #------ Atualizacao da barra de progresso ------#
  update.step <- max(5, floor(M.post/100))
  
  #------ Barra de progresso ------#
  pb <- txtProgressBar(max = M.post, style = 3)
  
  for (m in 1:M.post)
  {
    ft.k <- array(NA,dim=c(n.new,1,T.new+1))
    Qt.k <- array(NA,dim=c(n.new,n.new,T.new+1))
    Vt.c.k <- genVt(b0 = PSI.chain$b0[m],
                    b1 = PSI.chain$b1[m],
                    s2t = s2t.new,
                    phi = PSI.chain$phi[m],
                    DistMatrix = D)
    for (k in 2:(T.new+1))
    {
      ft.k <- t(Ft.new[,parc.prev,k])%*%Thetat.chain.final[1,,k]
      Qt.k <- symmetricMatrix(Vt.c.k[parc.prev,parc.prev,k])
      Xt.k <- as.numeric(mvnfast::rmvn(n = 1, mu = ft.k, sigma = Qt.k,isChol = TRUE))
      Yt.prev[,m,k] <- ifelse(Xt.k>=((csrd^PSI.chain$lambda[m])-1)/PSI.chain$lambda[m],
                              BCinv(Xt.k,lambda = PSI.chain$lambda[m]),0)
    }
    
    #------ Barra de progresso na tela ------#
    if((m%%update.step)==0){setTxtProgressBar(pb, m)}
  }
  
  ############################################
  ###### PREVISOES COMPLETAS PARA MAPAS ######
  ############################################
  
  # for (m in 1:M.post)
  # {
  #   print(m)
  #   ft.k <- array(NA,dim=c(n.new,1,T.new+1))
  #   at.k <- array(NA,dim=c(r,1,T.new+1))
  #   at.k[,,1] <- ht.final[m,]
  #   Qt.k <- array(NA,dim=c(n.new,n.new,T.new+1))
  #   Rt.k <- array(NA,dim=c(r,r,T.new+1))
  #   Rt.k[,,1] <- Ht.final[,,m]
  #   Vt.c.k <- genVt(b0 = PSI.chain$b0[m],
  #                   b1 = PSI.chain$b1[m],
  #                   s2t = s2t.new,
  #                   phi = PSI.chain$phi[m],
  #                   DistMatrix = D.new)
  #   for (k in 2:(T.new+1))
  #   {
  #     at.k[,,k] <- Gt%*%at.k[,,k-1]
  #     Rt.k[,,k] <- symmetricMatrix(Desc%*%Gt%*%Rt.k[,,k-1]%*%t(Gt)%*%Desc)
  #     ft.k <- t(Ft.new[,,k])%*%at.k[,,k]
  #     Qt.k <- symmetricMatrix(t(Ft.new[,,k])%*%Rt.k[,,k]%*%Ft.new[,,k]+Vt.c.k[,,k])
  #     Xt.k <- as.numeric(mvnfast::rmvn(n = 1, mu = ft.k, sigma = Qt.k,isChol = TRUE))
  #     Yt.prev[,m,k] <- ifelse(Xt.k>=((csrd^PSI.chain$lambda[m])-1)/PSI.chain$lambda[m],
  #                             BCinv(Xt.k,lambda = PSI.chain$lambda[m]),0)
  #   }
  # # ------ Barra de progresso na tela ------#
  # if((m%%update.step)==0){setTxtProgressBar(pb, m)}
  # }
  
  #------ Tempo total de processamento ------#
  b1 <- Sys.time()-a1
  print(b1)
  
  #----------------------------------------------------------------------------------------------
  # MEDIDAS DE COMPARACAO
  #----------------------------------------------------------------------------------------------
  
  # #------ Raiz do erro quadratico medio ------#
  # for (k in 2:(T.new+1))
  # {
  #   print(k - 1)
  #   # Media a posteriori
  #   print(sqrt(mean((Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean))^2)))
  #   # Mediana a posteriori
  #   print(sqrt(mean((Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k]))^2)))
  #   # Media do ensemble
  #   print(sqrt(mean((Yt.k[,,k] - Ft.new[2,1:n,k])^2)))
  #   # Ensemble
  #   print(sqrt(mean((Yt.k[,,k] - Ft.numforecast[2,1:n,k])^2)))
  # }
  # 
  # #------ Erro absoluto medio ------#
  # for (k in 2:(T.new+1))
  # {
  #   print(k - 1)
  #   # Media a posteriori
  #   print(mean(abs(Yt.k[,,k] - apply(Yt.prev[1:n,,k],1,mean))))
  #   # Mediana a posteriori
  #   print(mean(abs(Yt.k[,,k] - rowMedians(Yt.prev[1:n,,k]))))
  #   # Media do ensemble
  #   print(mean(abs(Yt.k[,,k] - Ft.new[2,1:n,k])))
  #   # Ensemble
  #   print(mean(abs(Yt.k[,,k] - Ft.numforecast[2,1:n,k])))
  # }
  # 
  # #------ Indice de concordancia de Willmott ------#
  # for (k in 2:(T.new+1))
  # {
  #   print(k - 1)
  #   # Media a posteriori
  #   print(ind.concord(sim=apply(Yt.prev[1:n,,k],1,mean), obs = Yt.k[,,k]))
  #   # Mediana a posteriori
  #   print(ind.concord(sim=rowMedians(Yt.prev[1:n,,k]), obs = Yt.k[,,k]))
  #   # Media do ensemble
  #   print(ind.concord(sim=Ft.new[2,1:n,k], obs = Yt.k[,,k]))
  #   # Ensemble
  #   print(ind.concord(sim=Ft.numforecast[2,1:n,k], obs = Yt.k[,,k]))
  # }
  
  # #------ Graficos - Concordancia ------#
  # par(mar=c(4,4,.5,.5))
  # par(mfrow=c(1,2))
  # plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
  #      0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
  #      xlab="Previsão numérica", ylab="Observado")
  # points(as.numeric(Ft.numforecast[2,1:n,2:(T.new+1)]),
  #        as.numeric(Yt.k[,,2:(T.new+1)]),pch=19,col=2)
  # abline(v=1:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
  #        h=1:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
  #        lty=2)
  # 
  # plot(0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
  #      0:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),t="l",
  #      xlab="Previsão numérica calibrada", ylab="Observado")
  # for (k in 2:(T.new+1))
  # {
  #   points(rowMedians(Yt.prev[1:n,,k]),
  #          as.numeric(Yt.k[,,k]),pch=19,col=2)
  # }
  # abline(v=1:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
  #        h=1:ceiling(max(as.numeric(Yt.k[,,2:(T.new+1)]))),
  #        lty=2)
  
  # #------ Graficos ------#
  # par(mfrow=c(2,2))
  # par(mar=c(4,4,.5,.5))
  # for (l in 1:59)
  # {
  #   plot(0,type="n",ylim=c(0,8),xlim=c(1,(T.new+2)),
  #        ylab=coords.stations$codigo[l],xlab="horizonte", axes= F)
  #   axis(1,1:6,c("-24","0","+24","+48","+72","+96"))
  #   axis(2,0:8,0:8)
  #   abline(h=0)
  #   polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
  #           c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.025)),
  #             colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.975)),
  #           col="grey80",border=NA)
  #   polygon(c(rev(3:(T.new+2)),3:(T.new+2)),
  #           c(rev(colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.1)),
  #             colQuantiles(Yt.prev[l,,2:(T.new+1)],probs = 0.9)),
  #           col="grey90",border=NA)
  #   points(1,ifelse(Yt[l,,T]>0,Yt[l,,T],0),pch=19,col=2)
  #   points(2,ifelse(Yt[l,,T+1]>0,Yt[l,,T+1],0),pch=19,col=2)
  #   points(3:(T.new+2),Yt.k[l,,2:(T.new+1)],pch=19,col=2)
  #   lines(3:(T.new+2),colMeans(Yt.prev[l,,2:(T.new+1)]),lty=2,lwd=2)
  #   lines(3:(T.new+2),Ft.numforecast[2,l,2:(T.new+1)],lty=2,lwd=2,col=4)
  # }
  # 
  # plot.new()
  # legend("bottom",legend = c("Obs","Eta","Calibrado","IC95%","IC90%"),
  #        lty = c(NA,2,2,NA,NA), lwd = c(2,2,2,15,15), pch= c(19,NA,NA,15,15),
  #        col = c(2,4,1,"grey80","grey90"))
  
  #----------------------------------------------------------------------------------------------
  # MAPAS - demanda tempo!
  #----------------------------------------------------------------------------------------------
  
  ############################################
  ############################################
  rm(D.new)
  rm(Vt.c.k)
  save.image(paste0("C:/Users/b207056565/Desktop/Calibration_EtaModel/Wind_Speed_10m/Aplicacao/Resultados/Diarios/result_diario",
                    datas[data.c],"SEMOS.Rda"))
  rm(list = subset(ls(),ls()!="data.c" & ls()!="datas"))
  ############################################
  ############################################
}

