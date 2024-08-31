#include <Rcpp.h>
using namespace Rcpp;

// This is a the SIR model coded in Rcpp
// Note that in Rcpp, you must 'declare' all variables 


// [[Rcpp::export]]
List SIR_cpp_model(NumericVector times,NumericVector state,
                            NumericVector parms) {
  
    // Define variables
    double S = state["S"];
    double I = state["I"];
    double R = state["R"];
    double N = S + I + R;
    
    // Extract parameters
    double beta = parms["beta"];
    double gamma = parms["gamma"];
    
    // Define differential equations
    double dS = - (beta * S * I) / N;
    double dI = (beta * S * I) / N - gamma * I;
    double dR = gamma * I;
    
    NumericVector res_vec = NumericVector::create(dS, dI, dR );
    
    List res = List::create(res_vec);
    
    return(res);
}


