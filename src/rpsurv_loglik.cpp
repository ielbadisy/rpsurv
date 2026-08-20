// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

// scale: 0 = PH (log cumulative hazard), 1 = PO (odds), 2 = probit
static const int SCALE_PH = 0;
static const int SCALE_PO = 1;
static const int SCALE_PROBIT = 2;

inline double softplus(double eta) {
  // numerically stable log(1 + exp(eta))
  if (eta > 0) return eta + std::log1p(std::exp(-eta));
  return std::log1p(std::exp(eta));
}

struct RPLogLikGrad : public Worker {
  const RMatrix<double> X;
  const RMatrix<double> dX;
  const RVector<double> status;
  const RVector<double> logtime;
  const RVector<double> beta;
  const int scale;

  double loglik;
  std::vector<double> grad;

  RPLogLikGrad(const NumericMatrix X_, const NumericMatrix dX_,
               const NumericVector status_, const NumericVector logtime_,
               const NumericVector beta_, int scale_)
    : X(X_), dX(dX_), status(status_), logtime(logtime_), beta(beta_), scale(scale_),
      loglik(0.0), grad(beta_.size(), 0.0) {}

  RPLogLikGrad(const RPLogLikGrad& other, Split)
    : X(other.X), dX(other.dX), status(other.status), logtime(other.logtime),
      beta(other.beta), scale(other.scale), loglik(0.0), grad(other.beta.size(), 0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t p = beta.size();
    for (std::size_t i = begin; i < end; i++) {
      double eta = 0.0, deta = 0.0;
      for (std::size_t k = 0; k < p; k++) {
        eta  += X(i, k)  * beta[k];
        deta += dX(i, k) * beta[k];
      }

      double logS, dlogS_deta, g, gprime;

      if (scale == SCALE_PH) {
        double H = std::exp(eta);
        logS = -H;
        dlogS_deta = -H;
        g = eta;
        gprime = 1.0;
      } else if (scale == SCALE_PO) {
        double sp = softplus(eta);
        double S = std::exp(-sp);
        logS = -sp;
        dlogS_deta = -(1.0 - S);
        g = eta - sp;
        gprime = S;
      } else { // probit
        double logS_ = R::pnorm(eta, 0.0, 1.0, /*lower*/0, /*log*/1);
        double logphi = R::dnorm(eta, 0.0, 1.0, /*log*/1);
        logS = logS_;
        dlogS_deta = -std::exp(logphi - logS_);
        g = logphi - logS_;
        gprime = -eta + std::exp(logphi - logS_);
      }

      double d = status[i];
      double li;
      if (d > 0.0) {
        if (deta <= 0.0) {
          // non-monotone hazard for these parameters: heavily penalise
          li = -1e10;
        } else {
          double logh = std::log(deta) - logtime[i] + g;
          li = d * logh + logS;
          for (std::size_t k = 0; k < p; k++) {
            double dli = d * (dX(i, k) / deta + gprime * X(i, k)) + dlogS_deta * X(i, k);
            grad[k] += dli;
          }
        }
      } else {
        li = logS;
        for (std::size_t k = 0; k < p; k++) {
          grad[k] += dlogS_deta * X(i, k);
        }
      }
      loglik += li;
    }
  }

  void join(const RPLogLikGrad& rhs) {
    loglik += rhs.loglik;
    for (std::size_t k = 0; k < grad.size(); k++) grad[k] += rhs.grad[k];
  }
};

//' Negative log-likelihood and gradient for a Royston-Parmar model
//'
//' @param beta parameter vector.
//' @param X design matrix for the linear predictor eta.
//' @param dX design matrix for d(eta)/d(log t).
//' @param logtime log of observed time.
//' @param status event indicator (1 = event, 0 = censored).
//' @param scale integer scale code: 0 = PH, 1 = PO, 2 = probit.
//' @return a list with `value` (negative log-likelihood) and `gradient`.
//' @keywords internal
// [[Rcpp::export]]
List rp_negloglik_grad_cpp(NumericVector beta, NumericMatrix X, NumericMatrix dX,
                            NumericVector logtime, NumericVector status, int scale) {
  RPLogLikGrad worker(X, dX, status, logtime, beta, scale);
  parallelReduce(0, X.nrow(), worker);

  NumericVector grad(worker.grad.size());
  for (std::size_t k = 0; k < worker.grad.size(); k++) grad[k] = -worker.grad[k];

  return List::create(
    Named("value") = -worker.loglik,
    Named("gradient") = grad
  );
}
