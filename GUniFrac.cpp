// [[Rcpp::depends(RcppArmadillo)]]

# include <RcppArmadillo.h>


using namespace Rcpp;


// [[Rcpp::export()]]
arma::cube GUniFracCpp2 (arma::mat cum,arma::mat brlen, arma::mat powV) {
  
  // inputs
  const int n = cum.n_cols ;
  const int powN = powV.n_cols;
  
  // containers
  arma::cube d(n,n,(powN+1)) ;
  d.zeros() ; // initialize distance cube to 0
  
  for (int i=1; i< n; i ++) {
    for (int j=0; j<=i-1; j++) {
      arma::vec  cum11 = cum.col(i);
      arma::vec cum21 = cum.col(j);
      arma::vec cum1 = cum11.elem(arma::find( (cum11+cum21)>0));
      arma::vec cum2 = cum21.elem(arma::find( (cum11+cum21)>0));
      arma::vec brlen2 = brlen.elem(arma::find((cum11+cum21)>0));
      
      arma::mat diff = abs(cum1 - cum2) / (cum1 + cum2);
      
      // calculate the distance for different weights
      for (int powi = 0; powi<powN; powi++) {
        arma::mat weight = pow(cum1+cum2,powV(0,powi));
        arma::mat  w = brlen2 % weight;
        d(i,j,powi) = accu(diff % w) /accu(w);
        d(j,i,powi) = d(i,j,powi);
      }
      arma::vec cum3 = arma::conv_to<arma::vec>::from(cum1>0);
      arma::vec cum4 = arma::conv_to<arma::vec>::from(cum2>0);
      arma::vec  diff2 = abs(cum3-cum4);
      d(i,j,powN) = accu(diff2 % brlen2)/accu(brlen2);
      d(j,i,powN) = d(i,j,powN);
   
    }
  }
  
  // returns
  return(d) ;
}

// [[Rcpp::export()]]
arma::mat cumall (arma::mat otu,arma::mat edge,arma::vec len) {
  
  // inputs
  const int ntip = otu.n_cols;
  const int n = otu.n_rows;
  const int nbr = edge.n_rows;
  
  // containers
  arma::vec edge2 = edge.col(1);
  arma::mat cum(nbr,n) ;
  cum.fill(0) ; // initialize d0 to 0
  
  for (int i=0; i< ntip; i ++) {
    arma::uvec  tiploc = arma::find(edge2 == (i+1));
    arma::rowvec tmpotu = otu.col(i).t();
    cum.rows(tiploc) =  cum.rows(tiploc) + tmpotu;
    arma::vec tmpnode = edge.col(0);
    arma::vec node = tmpnode.elem(tiploc);
    arma::uvec nodeloc = arma::find(edge2==node(0));
    while (nodeloc.n_elem>=1) {
      cum.rows(nodeloc) = cum.rows(nodeloc) + tmpotu;
      node = tmpnode.elem(nodeloc);
      nodeloc = find(edge2==node(0));
    }
  }
  
  // returns
  return(cum) ;
}



// [[Rcpp::export()]]
arma::mat BCdist (arma::mat X,arma::mat weight) {
  
  // inputs
  const int n = X.n_cols;
  // containers
  arma::mat bc(n,n) ;
  bc.fill(0) ; // initialize d0 to 0
  
  for (int i=1; i< n; i ++) {
    for (int j=0; j<=i-1; j++) {
      arma::mat  cum1 = X.col(i);
      arma::mat cum2 = X.col(j);
      arma::mat diff = abs(cum1 - cum2);
      bc(i,j) = accu(diff % weight) /accu(weight %(cum1 + cum2) );
      bc(j,i) = bc(i,j);
    }
  }
  
  // returns
  return(bc) ;
}
