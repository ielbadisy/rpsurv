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

// Scale-specific log-survival, its derivative wrt eta, and the hazard
// building blocks g(eta) = log h - log(deta/dlogt) + log t, g'(eta).
struct ScaleTerms {
  double logS;
  double dlogS_deta;
  double g;
  double gprime;
};

inline ScaleTerms scale_terms(double eta, int scale) {
  ScaleTerms t;
  if (scale == SCALE_PH) {
    double H = std::exp(eta);
    t.logS = -H;
    t.dlogS_deta = -H;
    t.g = eta;
    t.gprime = 1.0;
  } else if (scale == SCALE_PO) {
    double sp = softplus(eta);
    double S = std::exp(-sp);
    t.logS = -sp;
    t.dlogS_deta = -(1.0 - S);
    t.g = eta - sp;
    t.gprime = S;
  } else { // probit
    double logS_ = R::pnorm(eta, 0.0, 1.0, /*lower*/0, /*log*/1);
    double logphi = R::dnorm(eta, 0.0, 1.0, /*log*/1);
    t.logS = logS_;
    t.dlogS_deta = -std::exp(logphi - logS_);
    t.g = logphi - logS_;
    t.gprime = -eta + std::exp(logphi - logS_);
  }
  return t;
}

struct RPLogLikGrad : public Worker {
  const RMatrix<double> X;
  const RMatrix<double> dX;
  const RMatrix<double> Xentry;
  const RVector<double> hasEntry;
  const RVector<double> status;
  const RVector<double> logtime;
  const RVector<double> beta;
  const int scale;

  double loglik;
  std::vector<double> grad;

  RPLogLikGrad(const NumericMatrix X_, const NumericMatrix dX_, const NumericMatrix Xentry_,
               const NumericVector hasEntry_, const NumericVector status_,
               const NumericVector logtime_, const NumericVector beta_, int scale_)
    : X(X_), dX(dX_), Xentry(Xentry_), hasEntry(hasEntry_), status(status_), logtime(logtime_),
      beta(beta_), scale(scale_), loglik(0.0), grad(beta_.size(), 0.0) {}

  RPLogLikGrad(const RPLogLikGrad& other, Split)
    : X(other.X), dX(other.dX), Xentry(other.Xentry), hasEntry(other.hasEntry),
      status(other.status), logtime(other.logtime), beta(other.beta), scale(other.scale),
      loglik(0.0), grad(other.beta.size(), 0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t p = beta.size();
    for (std::size_t i = begin; i < end; i++) {
      double eta = 0.0, deta = 0.0;
      for (std::size_t k = 0; k < p; k++) {
        eta  += X(i, k)  * beta[k];
        deta += dX(i, k) * beta[k];
      }
      ScaleTerms ts = scale_terms(eta, scale);

      bool entry = hasEntry[i] > 0.0;
      double eta_e = 0.0;
      ScaleTerms te;
      if (entry) {
        for (std::size_t k = 0; k < p; k++) eta_e += Xentry(i, k) * beta[k];
        te = scale_terms(eta_e, scale);
      }
      // survival must not increase from entry to exit; a violation means the
      // fitted eta(log t) is non-monotone somewhere in (entry, exit] even
      // though it passes the pointwise check at the exit time below
      bool non_monotone_interval = entry && (ts.logS - te.logS > 1e-8);

      double d = status[i];
      double li;
      if (non_monotone_interval) {
        li = -1e10;
      } else if (d > 0.0) {
        if (deta <= 0.0) {
          li = -1e10;
        } else {
          double logh = std::log(deta) - logtime[i] + ts.g;
          li = d * logh + ts.logS;
          if (entry) li -= te.logS;
          for (std::size_t k = 0; k < p; k++) {
            double dli = d * (dX(i, k) / deta + ts.gprime * X(i, k)) + ts.dlogS_deta * X(i, k);
            if (entry) dli -= te.dlogS_deta * Xentry(i, k);
            grad[k] += dli;
          }
        }
      } else {
        li = ts.logS;
        if (entry) li -= te.logS;
        for (std::size_t k = 0; k < p; k++) {
          double dli = ts.dlogS_deta * X(i, k);
          if (entry) dli -= te.dlogS_deta * Xentry(i, k);
          grad[k] += dli;
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
//' Supports left truncation / counting-process data (genuine time-varying
//' covariates): each row contributes `logS(stop) - logS(entry)` to the
//' log-likelihood, `Xentry` and `hasEntry` giving the design row and flag
//' for the entry (left-truncation) time. Pass `hasEntry` all-zero (and
//' `Xentry` a matching-shape dummy matrix) for standard right-censored data.
//'
//' @param beta parameter vector.
//' @param X design matrix for the linear predictor eta at the exit (stop) time.
//' @param dX design matrix for d(eta)/d(log t) at the exit time.
//' @param Xentry design matrix for eta at the entry (start) time.
//' @param hasEntry 1 if the row is left-truncated (entry > 0), else 0.
//' @param logtime log of the exit (stop) time.
//' @param status event indicator (1 = event, 0 = censored) at the exit time.
//' @param scale integer scale code: 0 = PH, 1 = PO, 2 = probit.
//' @return a list with `value` (negative log-likelihood) and `gradient`.
//' @keywords internal
// [[Rcpp::export]]
List rp_negloglik_grad_cpp(NumericVector beta, NumericMatrix X, NumericMatrix dX,
                            NumericMatrix Xentry, NumericVector hasEntry,
                            NumericVector logtime, NumericVector status, int scale) {
  RPLogLikGrad worker(X, dX, Xentry, hasEntry, status, logtime, beta, scale);
  parallelReduce(0, X.nrow(), worker);

  NumericVector grad(worker.grad.size());
  for (std::size_t k = 0; k < worker.grad.size(); k++) grad[k] = -worker.grad[k];

  return List::create(
    Named("value") = -worker.loglik,
    Named("gradient") = grad
  );
}
