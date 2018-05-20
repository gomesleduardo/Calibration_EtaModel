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
##########|                     *[SIMULACAO]*                     |###########
##########|-------------------------------------------------------|###########
##########|-------------------------------------------------------|###########
##############################################################################
##############################################################################

#----------------------------------------------------------------------------------------------
# CARREGANDO PACOTES
#----------------------------------------------------------------------------------------------
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

#------------------------------------------ FUNCOES -------------------------------------------

Rcpp::sourceCpp("C:/Users/b207056565/Desktop/Dissertacao_Dados/Wind_Speed_10m3/Funcs_DGOP.cpp")

#---------------------------------------- DEFINICOES ------------------------------------------

#------ Fc M-H adapt ------#
MHadapt <- function (p, n, init, scale = rep(1, length(init)), adapt = !is.null(acc.rate), 
          acc.rate = NULL, gamma = 2/3, list = TRUE, showProgressBar = interactive(), 
          n.start = 0, ...) 
{
  if (adapt & !is.numeric(acc.rate)) 
    stop("Argument acc.rate is missing!")
  if (gamma <= 0.5 | gamma > 1) 
    stop("Argument gamma must be in (0.5, 1]!")
  if (is.numeric(adapt)) 
    n.adapt <- adapt
  if (adapt == TRUE) 
    n.adapt <- Inf
  if (adapt == FALSE) 
    n.adapt <- 0
  d <- length(init)
  X <- matrix(NA, ncol = d, nrow = n)
  colnames(X) <- names(init)
  X[1, ] <- init
  p.val <- rep(NA, n)
  val <- p(X[1, ], ...)
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
  if (ncol(M) != length(init)) 
    stop("Length or dimension of 'init' and 'scale' do not match!")
  S <- t(chol(M))
  # cat("  generate", n, "samples \\n")
  if (showProgressBar) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
  }
  update.step <- max(5, floor(n/100))
  k <- 0
  for (i in 2:n) {
    if (showProgressBar && i%%update.step == 0) {
      setTxtProgressBar(pb, i)
    }
    U <- rt(d, df = d)
    X.prop <- c(X[i - 1, ] + S %*% U)
    names(X.prop) <- names(init)
    val <- p(X.prop, ...)
    if (returns.list) {
      p.val.prop <- val$log.density
      extras.prop <- val["log.density" != names(val)]
    }
    else {
      p.val.prop <- val
    }
    alpha <- min(1, exp(p.val.prop - p.val[i - 1]))
    if (!is.finite(alpha)) 
      alpha <- 0
    if (runif(1) < alpha) {
      X[i, ] <- X.prop
      p.val[i] <- p.val.prop
      if (returns.list) {
        extras[[i]] <- extras.prop
      }
      k <- k + 1
    }
    else {
      X[i, ] <- X[i - 1, ]
      p.val[i] <- p.val[i - 1]
      if (returns.list) {
        extras[[i]] <- extras[[i - 1]]
      }
    }
    ii <- i + n.start
    if (ii < n.adapt) {
      adapt.rate <- min(1, d * ii^(-gamma))
      M <- S %*% (diag(d) + adapt.rate * (alpha - acc.rate) * 
                    U %*% t(U)/sum(U^2)) %*% t(S)
      eig <- eigen(M, only.values = TRUE)$values
      tol <- ncol(M) * max(abs(eig)) * .Machine$double.eps
      if (!isSymmetric(M) | is.complex(eig) | !all(Re(eig) > 
                                                   tol)) {
        M <- as.matrix(Matrix::nearPD(M)$mat)
      }
      S <- t(chol(M))
    }
  }
  if (showProgressBar) {
    close(pb)
  }
  acceptance.rate <- round(k/(n - 1), 3)
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
             Thetat = auxFFBS$Thetatchain,
             sigma2t = auxFFBS$sigma2tchain,
             phi = x[1],
             lambda = x[2],
             DistMatrix = D,
             csrd = csrd) +
      log.phi.prior(x[1]) + log.lambda.prior(x[2])
  } else {-9.99e+300}
}

#------ Total de tempos ------#
T <-  120

#------ Total de localizacoes ------#
n <-  60

#------ Total de membros no ensemble ------#
q <-  1

#------------------------------------#
#------ Parametros verdadeiros ------#
#------------------------------------#

#------ Fator de desconto para sigma2t ------#
delta.sig2 <- 1

#------ variancia dinamica: sigma2t via progressão beta ------#
set.seed(1)
sigma2t <- 2
for (t in 2:(T+1))
{
  sigma2t[t] <- delta.sig2*sigma2t[t-1]/rbeta(n = 1, shape1 = delta.sig2*((t-1)/2), 
                                              shape2 = (1-delta.sig2)*((t-1)/2))
}

#------ Grafico - sigma2t ------#
par(mar=c(4,5,.5,.5))
plot(sigma2t,t="l",ylab=expression(sigma[t]^2),xlab="Tempo")

#------ taxa de decaimento espacial ------#
phi <- 1.2

#------ transformacao Box-Cox ------#
lambda <- .8

#------ valor da censura ------#
csrd <- .1

#----------------------------------------------------------------------------------------------
# GERANDO DADOS
#----------------------------------------------------------------------------------------------

#------ Matriz de distancias D ------#
set.seed(999)
D <- as.matrix(dist(data.frame(x=runif(n,0,10),y=runif(n,0,10))))

#------ Matriz de evolucao Gt ------#
Gt <- as.matrix(bdiag(diag(1,q+1)))

#------ Comprimento do vetor Theta (Tendencia espacial + Sazonalidade) ------#
r <- dim(Gt)[1]

#------ Matriz de covariaveis Ft ------#
set.seed(1)
Ft <- array(NA,dim=c(r,n,T+1))
for (t in 2:(T+1)){Ft[,,t] <- t(as.matrix(
  data.frame(1,Ft=as.numeric(rmvnorm(1,
                                     mean = rep(4,n),
                                     sigma = sigma2t[t]*MatCor(phi = phi,DistMatrix = D))),
             t(rep(c(1,0),(r-q)/2))),n))}

#------ Matriz de covariancia da evolucao Wt ------#
set.seed(1)
Wt <- symmetricMatrix(matrix(rWishart(1,r,diag(c(0.0005,rep(.0001,r-1)))),r))
# Wt <- diag(c(rep(.0001,1),rep(.001,3),rep(.001,4)))

#------ Thetat ------#
set.seed(1)
Thetat <- data.frame(t(c(-.3,.7)))
for (t in 2:(T+1)){Thetat[t,] <- Gt%*%as.numeric(Thetat[t-1,]) + t(rmvnorm(1,mean=rep(0,r),sigma=Wt))}

#------ Grafico - Thetat ------#
par(mar=c(4,4,.5,.5))
par(mfrow=c(2,1))
for (m in 1:r){
  plot(Thetat[,m],t="l",xlab="Tempo",ylab=bquote(theta[.(m-1)]))}

#------ PGs verdadeiros ------#
set.seed(999)
Xt <- array(NA,dim=c(n,1,T+1))
for (t in 2:(T+1)){Xt[,,t] <- 
  as.numeric(rmvnorm(n = 1, 
                     mean = t(Ft[,,t])%*%as.numeric(Thetat[t,]),
                     sigma = sigma2t[t]*MatCor(phi = phi,DistMatrix = D)))}

#------ Transformacao Box-Cox inversa ------#
BCinv <- function(X_BCinv,lambda_BCinv){(1+(lambda_BCinv*X_BCinv))^(1/lambda_BCinv)}

#------ Valor observado transformado Y ------#
Yt <- array(NA,dim=c(n,1,T+1))
for (t in 2:(T+1)){Yt[,,t] <- ifelse(Xt[,,t]>=((csrd^lambda)-1)/lambda,BCinv(Xt[,,t],lambda = lambda),0)}

# #---------------------------- RECUPERANDO OS VALORES (VERIFICACAO) ----------------------------
# Yt2 <- array(NA,dim=c(n,1,T+1))
# for (t in 2:(T+1)){Yt2[,,t] <- ifelse(Xt[,,t]>=((csrd^lambda)-1)/lambda,BCinv(Xt[,,t],lambda = lambda),Xt[,,t])}
# 
# Xt2 <- BC(Yt = Yt2,lambda = lambda,csrd = csrd)
# 
# par(mar=c(4,4,.5,.5))
# par(mfrow=c(4,1))
# hist(Xt[,,-1],xlab="Latent Gaussian Process",main="",ylab=bquote(X[t]),freq = FALSE)
# hist(Yt[,,-1],xlab="Wind Speed (m/s) -- Censored",main="",ylab=bquote(Y[t]),freq = FALSE)
# hist(Yt2[,,-1],xlab="Wind Speed (m/s) -- Uncensored",main="",ylab=bquote(Y2[t]),freq = FALSE)
# hist(Xt2[,,-1],xlab="Latent Gaussian Process (recovered)",main="",ylab=bquote(X2[t]),freq = FALSE)
# par(mfrow=c(1,1))
# 
# #------ ANALISE DE CORRELACAO ENTRE PARAMETROS VIA FC DE VEROSSIMILHANCA CONJUNTA 2 A 2 ---------
# 
# #--------------------------#
# #------ phi e lambda ------#
# #--------------------------#
# 
# par(mfrow=c(1,1))
# #------ Dominio de variacao dos parametros ------#
# phi.range <-  seq(0.5,2.3, l=20)
# lambda.range <-  seq(.7,1.2, l=20)
# 
# #------ grade e fc. verossimilhanca bivariada ------#
# grid <-  expand.grid(phi.range,lambda.range)
# conj <-  0
# 
# #------ Barra de progresso ------#
# pb <- txtProgressBar(max = nrow(grid), style = 3)
# 
# #------ EMV iterativo sobre a grade ------#
# for (i in 1:nrow(grid))
# {
# 
#   #------ artificio matematico ------#
#   #------ exp(-12345,67) = exp(-12,34567)*exp(100) ------#
#   conj[i] <- exp(loglikep(Yt = Yt2,
#                           Ft = Ft,
#                           Thetat = as.matrix(Thetat),
#                           sigma2t = sigma2t,
#                           phi = grid[i,1],
#                           DistMatrix = D,
#                           lambda = grid[i,2],csrd = csrd)/100)*exp(100)
# 
#   if((i%%50)==0){setTxtProgressBar(pb, i)}
# }
# 
# #------ Fechando conexao da barra de progresso ------#
# close(pb)
# 
# #------ Linhas de contorno - fc. de verossimilhanca ------#
# grid <- data.frame(grid,conj)
# z <-  matrix(conj, ncol=length(phi.range), nrow =length(lambda.range))
# par(mar=c(4,4,.5,.5))
# contour(phi.range,lambda.range,z,nlevels = 10,drawlabels = FALSE,
#         main="",xlab=expression(phi),ylab=expression(lambda))
# abline(v=phi,h=lambda,lty=2,col=2)
# 
# #------ removendo arquivos desnecessarios ------#
# rm(Xt,Xt2,Yt2,grid,z,phi.range,lambda.range,conj,i,pb)
rm(t)

#----------------------------------------------------------------------------------------------
# MCMC
#----------------------------------------------------------------------------------------------

#------ Num de iteracoes M ------#
M <- 1100

#------ Burn-in ------#
burnin <- 100

#------ Salto/Afinamento (thinning) ------#
salto <- 200

#------ Proporcao de obs censuradas (ideal < 15%) ------#
print(paste0(round(length(which(Yt[,,-1] <= csrd))/length(as.numeric(Yt[,,-1]))*100,2),"%"))

#------------------#
#------ FFBS ------#
#------------------#

Desc <- diag(c(1/sqrt(.9),1/sqrt(.85)))

#------------------------------------#
#------ Distribuicoes a priori ------#
#------------------------------------#

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

# #------ mt e Ct final ------#
# ht.final <- matrix(0,M,r)
# Ht.final <- array(diag(r), dim = c(r,r,M))

#------ sigma2t final ------#
sigma2t.chain.final <- matrix(NA,T+1,M)

#------ phi e lambda = PSI ------#
PSI.chain <- data.frame(t(c(3,.5)))
names(PSI.chain) <- c("phi","lambda")

#------ Matriz de correlacao da proposta M-H ------#
EPS.PSI <- diag(((2.4)^2)/dim(PSI.chain)[2],dim(PSI.chain)[2])

#------ Vetor de PGs latentes correntes ------#
Xt <- BC(Yt = Yt, lambda = PSI.chain$lambda, csrd = csrd)

#------------------------------------------------#
#------ MCMC - Inicializacao de argumentos ------#
#------------------------------------------------#

#------ Matriz de correlacao corrente ------#
mat.cor.c <- MatCor(phi = PSI.chain$phi, DistMatrix = D)

#------ Contador de tempo ------#
a <- Sys.time()

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
  
  # ht.final[i,] <- auxFFBS$ht[,,T+1]
  # Ht.final[,,i] <- auxFFBS$Ht[,,T+1]
  Thetat.chain.final[,,i] <- auxFFBS$Thetatchain
  Xt.mu <- auxFFBS$Xtmu
  sigma2t.chain.final[,i] <- auxFFBS$sigma2tchain

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
  
  #------ Atualizacao da matriz de correlacao corrente Vt(k) ------#
  if (all.equal(as.numeric(PSI.chain[i,]),as.numeric(PSI.chain[i-1,]))!=TRUE)
  {
    #------ Matriz de correlacao corrente ------#
    mat.cor.c <- MatCor(phi = PSI.chain$phi[i], DistMatrix = D)
  }
  
  #--------------------------------------------------#
  #------ Amostrador de Gibbs p/ Z var latente ------#
  #--------------------------------------------------#
  for (t in 2:(T+1))
  {
    if (length(which(Yt[,,t] <= csrd))>0)
    {
      cens.c <- which(Yt[,,t] <= csrd)
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

#------ phi ------#
par(mar=c(4,4,.5,.5))
par(mfrow=c(2,2))
plot(sampMH$samples[,1],type="l",ylab=expression(phi),xlab="Iteracoes")
abline(h=phi,col=2,lwd=2)
acf(sampMH$samples[,1],lag.max = 25,main="",ylab="",xlab="Defasagem")

#------ lambda ------#
plot(sampMH$samples[,2],type="l",ylab=expression(lambda),xlab="Iteracoes")
abline(h=lambda,col=2,lwd=2)
acf(sampMH$samples[,2],lag.max = 25,main="",ylab="",xlab="Defasagem")

#-------------------------------------------------#
#------ Amostras da distribuicao posteriori ------#
#-------------------------------------------------#

#------ Thetat ------#
par(mar=c(4,4,.5,.5))
par(mfrow=c(2,1))
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
  lines(Thetat[,m],col=2,lwd=2)
  abline(h=0,lty=3)
}

#------ phi ------#
par(mar=c(4,4,.5,.5))
par(mfrow=c(2,3))
plot(PSI.chain$phi[-(1:burnin)],type="l",ylab=expression(phi),xlab="Iteracoes")
abline(h=phi,col=2,lwd=2)
acf(PSI.chain$phi[-(1:burnin)],lag.max = 20,main="",ylab="",xlab="Defasagem")
d <- density(x = PSI.chain$phi[-(1:burnin)], adjust = 2)
plot(d, xlab=expression(phi), main="",ylab="Densidade")
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$phi[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.phi.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)
abline(v=phi,col=2,lwd=2)

#------ lambda ------#
plot(PSI.chain$lambda[-(1:burnin)],type="l",ylab=expression(lambda),xlab="Iteracoes")
abline(h=lambda,col=2,lwd=2)
acf(PSI.chain$lambda[-(1:burnin)],lag.max = 20,main="",ylab="",xlab="Defasagem")
d <- density(x = PSI.chain$lambda[-(1:burnin)], adjust = 2)
plot(d, xlab=expression(lambda), main="",ylab="Densidade")
h1 <- d$x[d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975)))] # ?rea a ser pintada
h2 <- d$y[which(d$x >= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.025))) & d$x <= as.numeric(quantile(PSI.chain$lambda[-(1:burnin)],c(0.975))))]
polygon(c(min(h1),h1,max(h1)),c(0,h2,0),col="gray80")
lines(d,lwd=2)
lines(seq(min(d$x),max(d$x),.01),exp(log.lambda.prior(seq(min(d$x),max(d$x),.01))),lwd=2,lty=2)
abline(h=0)
abline(v=lambda,col=2,lwd=2)