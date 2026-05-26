#include <TMB.hpp>
#include <math.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // Vastemuuttujan matriisi
  DATA_MATRIX(y);
  // Ympäristökovariaattien design-matriisi 
  DATA_MATRIX(x);
  // SPDE-objekti, jonka avulla muodostetaan keskittyneisyysmatriisi latenteille muuttujille
  DATA_STRUCT(spde, spde_t);
  // SPDE-verkko-projektio matriisi 
  DATA_SPARSE_MATRIX(projmat);
  
  // Regressiokertoimet
  PARAMETER_MATRIX(b);
  // Satunnaisefektit eli latentit muuttujat
  PARAMETER_ARRAY(u);
  // NB-jakauman lajikohtainen hajontaparametri
  PARAMETER_VECTOR(lg_phi);
  // Matern-kovarianssifunktion skaalaparametri
  PARAMETER(lg_range); 
  // Latausmatriisi
  PARAMETER_VECTOR(lataus);
  // Latenttien muuttujien kovarianssimatriisi
  PARAMETER_VECTOR(lv_cor);
  
  // Vastemuuttujamatriisin dimensiot
  int n = y.rows();
  int p = y.cols();
  // Latenttien muuttujien lukumäärä
  int d = u.cols();
  
  Type nu = 1.0;  // Matern-kovarianssifunktion muotoparametri nu
  Type pi = 3.1415926536;
  // Matern-kovarianssifunktion kappa-parametri
  Type lg_kappa = 0.5*log(8*nu) - lg_range;
  // Matern-kovarianssifunktion tau-parametri
  Type lg_tau = -0.5*log(4*pi) - lg_kappa;  
  
  
  // Alustetaan lineaarinen ennustin
  matrix<Type> eta(n,p);
  eta.fill(0.0);
  
  // Lisätään regressiokertoimet lineaariseen ennustimeen
  eta += x*b;
  
  // Määritetään latausmatriisi, joka on yläkolmiomatriisi, jonka diagonaali on positiivinen
  int k = 0;
  matrix<Type> lam(d,p); 
  for (int i=0; i < d; i++){
    for (int j=0; j < p; j++){
      if (i > j) {
        lam(i,j) = 0.0;
      }
      else if (i == j) {
        lam(i,j) = exp(lataus(k));
        k++;
      }
      else {
        lam(i,j) = lataus(k);
        k++;
      }
    }
  }
  
  // Funktio latausmatriisin rakenteen tarkastukseen
  REPORT(lam);
  
  // Lisätään satunnaisefektit lineaariseen ennustimeen
  eta += projmat*(u.matrix()*lam);
  
  // Log-uskottavuuden alkuarvo
  Type nll = 0.0;
  
  // Latenttien muuttujien kovarianssimatriisi
  UNSTRUCTURED_CORR_t<Type> C(lv_cor);
  
  // Muodostetaan harva keskittyneisyysmatriisi SPDE-objektin ja Matern-kovarianssin parametrin avulla
  SparseMatrix<Type> Q = Q_spde(spde, exp(lg_kappa)); 
  
  // Muodostetaan separoituva kovarianssirakenne satunnaisefekteille
  nll += SEPARABLE(C, SCALE(GMRF(Q), 1/exp(lg_tau)))(u);  
  
  // Mallinnetaan vastemuuttujaa negatiivisella binomijakaumalla
  for (int i=0; i<n; i++) {
    for (int j=0; j<p; j++) {
      nll -= dnbinom_robust(y(i,j), eta(i,j), 2*eta(i,j) - lg_phi(j), 1); 
    }
  }
  
  return nll;
}
