// calculates function tau_n as used in asymptotic variance
// calculations for each observation

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fun_tau_n(NumericVector H, NumericVector Y) {
    int n = H.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++){
        out[i] = (double)(sum(Y > Y[i] & H > H[i]) + 
                sum(Y < Y[i] & H < H[i])) / n;    
    }
    return out;
}
