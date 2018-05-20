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
      if (arma::as_scalar(YtBC(l,0,t)) >= csrd)
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

// FUNCAO: GERACAO DA MATRIZ Vt
arma::cube genVt(double b0,
                 double b1,
                 arma::cube s2t,
                 double phi,
                 arma::mat DistMatrix)
  
{
  // DECLARACAO DE ARGUMENTOS
  int n = DistMatrix.n_rows;
  int T = s2t.n_slices - 1;
  arma::cube s2tVt((const double*)s2t.begin(),1,n,T+1);
  arma::mat D((const double*)DistMatrix.begin(),n,n);
  arma::mat Cor(n,n);
  arma::colvec s2aux(n);
  arma::cube s2tfinal(n,n,T+1);
  
  // CALCULOS
  
  // matriz de correlacao
  Cor = MatCor(phi,D);
  
  // matriz de covariancia NREMOS
  for (int t = 1; t < (T+1); t++) 
  {
    for (int c = 0; c < n; c++) 
    {
      s2aux(c) = std::sqrt(b0 + (b1*arma::as_scalar(s2tVt(0,c,t))));
    }
    s2tfinal.slice(t) = (arma::diagmat(s2aux))*Cor*(arma::diagmat(s2aux));
  }
  
  // RESULTADOS
  return s2tfinal;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

// FUNCAO: LOG VEROSSIMILHANCA PROPORCIONAL
double loglikep(arma::cube Yt,
                arma::cube Ft,
                arma::mat  Thetat,
                arma::cube Vt,
                double lambda,
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
  arma::cube VtML((const double*)Vt.begin(),n,n,T+1);
  double sumlogdet;
  double normkern;
  double jacob;
  
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
  
  // somatorio dos log determinantes
  sumlogdet = 0;
  for (int t = 1; t < (T+1); t++) {sumlogdet = sumlogdet + arma::as_scalar(std::log(det(VtML.slice(t))));}
  
  // nucleo da Normal
  normkern = 0;
  for (int t = 1; t < (T+1); t++) 
  {
    normkern = normkern + arma::as_scalar(((Xt.slice(t)-(FtML.slice(t).t()*ThetatML.row(t).t())).t())*
      inv(VtML.slice(t))*
      (Xt.slice(t)-(FtML.slice(t).t()*ThetatML.row(t).t())));
  }
  
  //RESULTADOS
  return (-.5*sumlogdet -.5*normkern + jacob);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

// FUNCAO: FORWARD FILTERING BACKWARD SAMPLING
List FFBS(arma::cube Xt,
          arma::cube Ft,
          arma::mat Gt,
          arma::cube Vt,
          arma::mat DescMatrix)
  
{
  // DECLARACAO DE ARGUMENTOS
  int n = Xt.n_rows;
  int T = Xt.n_slices - 1;
  int r = Ft.n_rows;
  arma::cube XtFFBS((const double*)Xt.begin(),n,1,T+1);
  arma::cube FtFFBS((const double*)Ft.begin(),r,n,T+1);
  arma::mat  GtFFBS((const double*)Gt.begin(),r,r);
  arma::cube VtFFBS((const double*)Vt.begin(),n,n,T+1);
  arma::mat Desc((const double*)DescMatrix.begin(),r,r);
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
  arma::mat Thetatchain(T+1,r);
  arma::mat Maux(r,r);
  arma::cube Xtmu(n,1,T+1);
  
  // FORWARD FILTERING (PROGRESSIVO)
  
  // prioris em t=0
  Ct.slice(0) = Ct.slice(0).eye();
  
  for (int t = 1; t < (T+1); t++) 
  {
    // Priori em t
    at.slice(t) = GtFFBS * mt.slice(t-1);
    Rt.slice(t) = (Desc * (GtFFBS * Ct.slice(t-1) * GtFFBS.t()) * Desc);
    Rt.slice(t) = symmetricMatrix(as<arma::mat>(makeposdef(wrap(Rt.slice(t)))));
    
    // Previsao 1 passo a frente
    ft.slice(t) = FtFFBS.slice(t).t() * at.slice(t);
    Qt.slice(t) = (FtFFBS.slice(t).t() * Rt.slice(t) * FtFFBS.slice(t)) + VtFFBS.slice(t);
    Qt.slice(t) = symmetricMatrix(as<arma::mat>(makeposdef(wrap(Qt.slice(t)))));
    
    // Posteriori em t
    At.slice(t) = Rt.slice(t) * FtFFBS.slice(t) * inv_sympd(Qt.slice(t));
    et.slice(t) = XtFFBS.slice(t) - ft.slice(t);
    mt.slice(t) = at.slice(t) + (At.slice(t) * et.slice(t));
    Ct.slice(t) = Rt.slice(t) - (At.slice(t) * Qt.slice(t) * At.slice(t).t());
    Ct.slice(t) = symmetricMatrix(as<arma::mat>(makeposdef(wrap(Ct.slice(t)))));
  }
  
  // BACKWARD SAMPLING (REGRESSIVO)
  for (int t = T; t > -1; t--) 
  {
    if (t == T)
    {
      ht.slice(t) = mt.slice(t);
      Ht.slice(t) = Ct.slice(t);
    } else {
      Maux = Ct.slice(t) * GtFFBS.t() * inv_sympd(Rt.slice(t+1));
      ht.slice(t) = mt.slice(t) + (Maux * (Thetatchain.row(t+1).t() - at.slice(t+1)));
      Ht.slice(t) = Ct.slice(t) - (Maux * GtFFBS * Ct.slice(t));
      Ht.slice(t) = symmetricMatrix(as<arma::mat>(makeposdef(wrap(Ht.slice(t)))));
    }
    // Amostrando Thetat(k)|Thetat(k-1) ~ N(ht,Ct)
    Thetatchain.row(t) = as<arma::rowvec>(RMVN(wrap(ht.slice(t)),wrap(Ht.slice(t))));
  }
  
  for (int t = 1; t < (T+1); ++t)
  {Xtmu.slice(t) = FtFFBS.slice(t).t()*Thetatchain.row(t).t();}
  
  //RESULTADOS
  return List::create(Named("mt") = mt,
                      Named("Ct") = Ct,
                      Named("at") = at,
                      Named("Rt") = Rt,
                      Named("ht") = ht,
                      Named("Ht") = Ht,
                      Named("Xtmu") = Xtmu,
                      Named("Thetatchain") = Thetatchain);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double RTN(NumericVector mean, 
           NumericVector sd, 
           NumericVector csrd)
{
  NumericVector n;
  NumericVector li;
  n = 1;
  li = -10;
  
  // NOME DO PACOTE QUE CONTEM A FUNCAO
  Rcpp::Environment package_env("package:truncnorm"); 
  
  // FAZENDO A FUNCAO SER CHAMADA
  Rcpp::Function rfunction3 = package_env["rtruncnorm"];    
  
  // RESULTADO
  return as<double>(rfunction3(n,li,csrd,mean,sd));
}

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