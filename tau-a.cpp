// calculates kendall's tau-a

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double taua(NumericVector H, NumericVector Y) {
	NumericVector H0 = H[Y == 0];
	NumericVector H1 = H[Y == 1];
	int n0 = H0.size();
	int n1 = H1.size();
	double num1 = 0;  // concordant pair count
	// int num2 = 0;  // discordant pair count

	for (int i = 0; i < n0; i++) {
		for (int j = 0; j < n1; j++) {
			if (H0[i] < H1[j]) {
				num1 += 1;
			} else if (H0[i] == H1[j]) {
				num1 += 0.5;
			} else {
				continue;
			}
			// if (H0[i] > H1[j]) num2 ++;
		}
	}

	double res = (double) num1 / (n0 * n1);
	return res;
}