#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

// FUNCAO: MATRIZ DE CORRELACAO
arma::mat MatCor(double phi,
                 arma::mat DistMatrix)
  
{
  //RESULTADOS
  return arma::exp(-phi*DistMatrix);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

// FUNCAO: TRANSFORMACAO BOX COX
arma::cube BC(arma::cube Yt,
              double lambda,
              double csrd)
  
{
  // DECLARACAO DE ARGUMENTOS
  int n = Yt.n_rows;
  int T = Yt.n_slices - 1;
  arma::cube YtBC((const double*)Yt.begin(),n,1,T+1);
  arma::cube Xt(n,1,T+1);
  
  // CALCULOS
  for (int t = 1; t < (T+1); t++) 
  {
    for (int l = 0; l < n; l++)
    {
      if (arma::as_scalar(YtBC(l,0,t)) > csrd)
        {Xt(l,0,t) = (std::pow(arma::as_scalar(YtBC(l,0,t)),lambda)-1)/lambda;} 
      else {Xt(l,0,t) = arma::as_scalar(Yt(l,0,t));}
    }
  }
  
  // RESULTADOS
  return Xt;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

// FUNCAO: MATRIZ SIMETRICA
arma::mat symmetricMatrix(arma::mat matrix)
  
{
  // DECLARACAO DE ARGUMENTOS
  int n = matrix.n_rows;
  arma::mat  B((const double*)matrix.begin(),n,n);
  
  //RESULTADOS
  return (trimatl(B) + trimatl(B).t() - diagmat(B));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix makeposdef(NumericMatrix matrix){
  
  // NOME DO PACOTE QUE CONTEM A FUNCAO
  Rcpp::Environment package_env("package:corpcor"); 
  
  // FAZENDO A FUNCAO SER CHAMADA
  Rcpp::Function rfunction = package_env["make.positive.definite"];    
  
  // RESULTADO
  return rfunction(matrix);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector RMVN(NumericVector mean, NumericMatrix covmatrix){
  
  // NOME DO PACOTE QUE CONTEM A FUNCAO
  Rcpp::Environment package_env("package:mvtnorm"); 
  
  // FAZENDO A FUNCAO SER CHAMADA
  Rcpp::Function rfunction2 = package_env["rmvnorm"];    
  
  // RESULTADO
  return rfunction2(1,mean,covmatrix,"svd");
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector RGama(NumericVector shape, NumericVector rate){
  
  // NOME DO PACOTE QUE CONTEM A FUNCAO
  Rcpp::Environment package_env("package:stats"); 
  
  // FAZENDO A FUNCAO SER CHAMADA
  Rcpp::Function rfunction4 = package_env["rgamma"];    
  
  // RESULTADO
  return rfunction4(1,shape,rate);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

// FUNCAO: FORWARD FILTERING BACKWARD SAMPLING - Variancia descontada e Matriz de descontos
List FFBSvar(arma::cube Xt, 
             arma::cube Ft,
             arma::mat Gt,
             arma::mat CorMatrix,
             arma::mat DescMatrix,
             double delta)
  
{
  // DECLARACAO DE ARGUMENTOS
  int n = Xt.n_rows;
  int T = Xt.n_slices - 1;
  int r = Ft.n_rows;
  arma::cube XtFFBS((const double*)Xt.begin(),n,1,T+1);
  arma::cube FtFFBS((const double*)Ft.begin(),r,n,T+1);
  arma::mat  GtFFBS((const double*)Gt.begin(),r,r);
  arma::mat  Cor((const double*)CorMatrix.begin(),n,n);
  arma::mat  DescFFBS((const double*)DescMatrix.begin(),r,r);
  arma::cube mt(r,1,T+1); mt.zeros();
  arma::cube Ct(r,r,T+1);
  arma::cube ht(r,1,T+1);
  arma::cube Ht(r,r,T+1);
  arma::cube at(r,1,T+1);
  arma::cube Rt(r,r,T+1);
  arma::cube ft(n,1,T+1);
  arma::cube Qt(n,n,T+1);
  arma::cube At(r,n,T+1);
  arma::cube et(n,1,T+1);
  arma::vec nt(T+1);
  arma::vec St(T+1);
  arma::vec nt2(T+1);
  arma::vec St2(T+1);
  arma::mat Maux(r,r);
  arma::cube Xtmu(n,1,T+1);
  arma::mat Thetatchain(T+1,r);
  arma::vec sigma2tchain(T+1);
  
  // FORWARD FILTERING (PROGRESSIVO)
  
  // prioris em t=0
  sigma2tchain(0) = 0.0000;
  Ct.slice(0) = Ct.slice(0).eye();
  nt(0) = 1.00000000;
  St(0) = 1.00000000;
  
  for (int t = 1; t < (T+1); t++) 
  {
    // Priori em t
    at.slice(t) = GtFFBS * mt.slice(t-1);
    Rt.slice(t) = symmetricMatrix(DescFFBS * (GtFFBS * Ct.slice(t-1) * GtFFBS.t()) * DescFFBS);
    Rt.slice(t) = symmetricMatrix(as<arma::mat>(makeposdef(wrap(Rt.slice(t)))));
    
    // Previsao 1 passo a frente
    ft.slice(t) = FtFFBS.slice(t).t() * at.slice(t);
    Qt.slice(t) = symmetricMatrix((FtFFBS.slice(t).t() * Rt.slice(t) * FtFFBS.slice(t)) + arma::as_scalar(St(t-1))*Cor);
    Qt.slice(t) = symmetricMatrix(as<arma::mat>(makeposdef(wrap(Qt.slice(t)))));
    
    // Posteriori em t
    At.slice(t) = Rt.slice(t) * FtFFBS.slice(t) * inv_sympd(Qt.slice(t));
    et.slice(t) = XtFFBS.slice(t) - ft.slice(t);
    mt.slice(t) = at.slice(t) + (At.slice(t) * et.slice(t));
    nt(t) = delta*arma::as_scalar(nt(t-1)) + 1.0000000;
    St(t) = (1.0000/arma::as_scalar(nt(t)))*(delta*arma::as_scalar(St(t-1))*arma::as_scalar(nt(t-1)) +
      ((1.0000000/n)*arma::as_scalar(St(t-1))*arma::as_scalar(et.slice(t).t()*inv_sympd(Qt.slice(t))*et.slice(t))));
    Ct.slice(t) = symmetricMatrix(arma::as_scalar(arma::as_scalar(St(t))/arma::as_scalar(St(t-1)))*(Rt.slice(t) - (At.slice(t) * Qt.slice(t) * At.slice(t).t())));
    Ct.slice(t) = symmetricMatrix(as<arma::mat>(makeposdef(wrap(Ct.slice(t)))));
  }
  
  // BACKWARD SAMPLING (REGRESSIVO)
  for (int t = T; t > -1; t--) 
  {
    if (t == T)
    {
      ht.slice(t) = mt.slice(t);
      Ht.slice(t) = Ct.slice(t);
      nt2(t) = nt(t); 
      St2(t) = St(t);
    } else {
      Maux = Ct.slice(t) * GtFFBS.t() * inv_sympd(Rt.slice(t+1));
      ht.slice(t) = mt.slice(t) + (Maux * (Thetatchain.row(t+1).t() - at.slice(t+1)));
      Ht.slice(t) = symmetricMatrix(Ct.slice(t) - (Maux * GtFFBS * Ct.slice(t)));
      Ht.slice(t) = symmetricMatrix(as<arma::mat>(makeposdef(wrap(Ht.slice(t)))));
      nt2(t) = (1.000-delta)*arma::as_scalar(nt(t)) + delta*arma::as_scalar(nt2(t+1)); 
      St2(t) = 1.0000/(((1.0000 - delta)*(1.00000/arma::as_scalar(St(t)))) +
        (delta*(1.00000/arma::as_scalar(St2(t+1))))); 
    }
    // Amostrando Thetat(k)|Thetat(k-1) ~ N(ht,Ct)
    Thetatchain.row(t) = as<arma::rowvec>(RMVN(wrap(ht.slice(t)),wrap(Ht.slice(t))));
    
    // Amostrando sigma2t(k)|sigma2t(k-1) ~ GI(nt,ntSt)
    sigma2tchain(t) = 1.00000/arma::as_scalar(as<arma::rowvec>(RGama(wrap(nt2(t)),wrap(nt2(t)*St2(t)))));
  }
  
  for (int t = 1; t < (T+1); ++t)
  {Xtmu.slice(t) = FtFFBS.slice(t).t()*Thetatchain.row(t).t();}
  
  //RESULTADOS
  return List::create(Named("ht") = ht,
                      Named("Ht") = Ht,
                      Named("Xtmu") = Xtmu,
                      Named("Thetatchain") = Thetatchain,
                      Named("nt2") = nt2,
                      Named("St2") = St2,
                      Named("sigma2tchain") = sigma2tchain);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

// ---------------------------- LOG VEROSSIMILHANCA PROPORCIONAL ----------------------------
double loglikep(arma::cube Yt,
                arma::cube Ft,
                arma::mat Thetat,
                arma::colvec sigma2t,
                double lambda,
                double phi,
                arma::mat DistMatrix,
                double csrd)
  
{
  // DECLARACAO DE ARGUMENTOS
  int n = Yt.n_rows;
  int T = Yt.n_slices - 1;
  int r = Ft.n_rows;
  arma::cube YtML((const double*)Yt.begin(),n,1,T+1);
  arma::cube Xt(n,1,T+1); Xt = BC(YtML,lambda,csrd);
  arma::mat ThetatML((const double*)Thetat.begin(),T+1,r);
  arma::cube FtML((const double*)Ft.begin(),r,n,T+1);
  arma::colvec sigma2tML((const double*)sigma2t.begin(),T+1);
  arma::mat D((const double*)DistMatrix.begin(),n,n);
  double normkern;
  double jacob;
  arma::mat CorInv(n,n);
  
  // CALCULOS
  
  // jacobiano
  jacob = 0;
  for (int t = 1; t < (T+1); t++) 
  {
    for (int l = 0; l < n; l++)
    {
      if (arma::as_scalar(YtML(l,0,t)) > csrd)
      {
        jacob = jacob + std::log(arma::as_scalar(YtML(l,0,t)));
      }
    }
  }
  jacob = (lambda-1)*jacob;
  
  // matriz de correlacao inversa
  CorInv = arma::inv(MatCor(phi,D));
  
  // nucleo da Normal
  normkern = 0;
  for (int t = 1; t < (T+1); t++) 
  {
    normkern = normkern + arma::as_scalar(((Xt.slice(t)-(FtML.slice(t).t()*ThetatML.row(t).t())).t())*
      (arma::as_scalar(1.00000/arma::as_scalar(sigma2tML(t)))*CorInv)*
      (Xt.slice(t)-(FtML.slice(t).t()*ThetatML.row(t).t())));
  }
  
  //RESULTADOS
  return ((-n*arma::sum(sigma2tML) - 
          (T/2.0000)*std::log(arma::det(MatCor(phi,D)))
             -.50000*normkern + jacob));
}

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// 
// double RTN(NumericVector mean, 
//            NumericVector sd, 
//            NumericVector csrd)
// {
//   NumericVector n;
//   NumericVector li;
//   n = 1;
//   li = -10;
//   
//   // NOME DO PACOTE QUE CONTEM A FUNCAO
//   Rcpp::Environment package_env("package:truncnorm"); 
//   
//   // FAZENDO A FUNCAO SER CHAMADA
//   Rcpp::Function rfunction3 = package_env["rtruncnorm"];    
//   
//   // RESULTADO
//   return as<double>(rfunction3(n,li,csrd,mean,sd));
// }

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// 
// // FUNCAO: GERACAO DA VARIAVEL LATENTE Zt
// List VarLatZ(arma::uvec In,
//              arma::uvec Out,
//              arma::mat Vt,
//              arma::vec Xt,
//              arma::vec Xtmean,
//              NumericVector csrd){
//   
//   // DECLARACAO DE ARGUMENTOS
//   int n = Xt.n_rows;
//   arma::mat  VtZ((const double*)Vt.begin(),n,n);
//   arma::vec  XtZ((const double*)Xt.begin(),n);
//   arma::vec  Xtmu((const double*)Xtmean.begin(),n);
//   arma::uvec  inZ(In.begin(),1);
//   arma::uvec  outZ(Out.begin(),n-1);
//   arma::mat  parcZ(1,n-1);
//   NumericVector Zmu;
//   NumericVector Zsigma;
//   
//   // CALCULOS
//   parcZ = VtZ(inZ,outZ)*inv_sympd(VtZ(outZ,outZ)); 
//   Zmu = arma::as_scalar(Xtmu(inZ) + parcZ*(XtZ(outZ)-Xtmu(outZ)));
//   Zsigma = arma::as_scalar(std::sqrt(arma::as_scalar(VtZ(inZ,inZ) - parcZ*VtZ(outZ,inZ))));
//   
//   //RESULTADOS
//   return List::create(Named("Z.mu") = Zmu,
//                       Named("Z.sigma") = Zsigma,
//                       Named("Z.random") = RTN(Zmu, Zsigma,csrd));
// }