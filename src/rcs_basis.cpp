#include <Rcpp.h>
using namespace Rcpp;

inline double pos3(double u) { return u > 0 ? u * u * u : 0.0; }
inline double pos2(double u) { return u > 0 ? u * u : 0.0; }

//' Restricted cubic spline basis, computed in C++
//'
//' Same basis as `rcs_basis()` (Durrleman-Simon / Royston-Parmar
//' parameterisation): a linear term plus one term per interior knot.
//'
//' @param x numeric vector (typically log time).
//' @param knots numeric vector of ALL knots (boundary knots first/last,
//'   interior knots in between), already sorted.
//' @param derivative if `TRUE`, return d(basis)/dx.
//' @return an `n x (length(knots) - 1)` matrix.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rcs_basis_cpp(NumericVector x, NumericVector knots, bool derivative) {
  int nknots = knots.size();
  double kmin = knots[0];
  double kmax = knots[nknots - 1];
  int nk = nknots - 2;
  int n = x.size();

  NumericMatrix out(n, nk + 1);

  if (!derivative) {
    for (int i = 0; i < n; i++) out(i, 0) = x[i];
    for (int j = 0; j < nk; j++) {
      double kj = knots[j + 1];
      double lambda_j = (kmax - kj) / (kmax - kmin);
      for (int i = 0; i < n; i++) {
        out(i, j + 1) = pos3(x[i] - kj) - lambda_j * pos3(x[i] - kmin) - (1.0 - lambda_j) * pos3(x[i] - kmax);
      }
    }
  } else {
    for (int i = 0; i < n; i++) out(i, 0) = 1.0;
    for (int j = 0; j < nk; j++) {
      double kj = knots[j + 1];
      double lambda_j = (kmax - kj) / (kmax - kmin);
      for (int i = 0; i < n; i++) {
        out(i, j + 1) = 3.0 * pos2(x[i] - kj) - lambda_j * 3.0 * pos2(x[i] - kmin) - (1.0 - lambda_j) * 3.0 * pos2(x[i] - kmax);
      }
    }
  }

  return out;
}
