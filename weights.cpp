#include <Eigen/Dense>

std::complex<double> semidisk(const double t) {
  const std::complex<double> iu(0, 1);

  if (0 <= t && t < 1)
    return iu * t;
  if (1 <= t && t < M_PI + 1)
    return std::exp(iu * (t + M_PI_2 - 1));
  if (M_PI + 1 <= t && t < M_PI + 2)
    return iu * (t - (M_PI + 2));

  return std::complex<double>(0, 0);
}

class WeightsBuilder {
public:
  WeightsBuilder(const int, const int, const double);

  Eigen::VectorXd predictors() { return compute_coeff(predictor_matrix()); }
  Eigen::VectorXd correctors() { return compute_coeff(corrector_matrix()); }

private:
  Eigen::VectorXcd lambdas;
  Eigen::VectorXd times;
  double timestep, future_time;

  Eigen::MatrixXcd predictor_matrix() const;
  Eigen::MatrixXcd corrector_matrix() const;
  Eigen::VectorXcd rhs_vector() const;

  Eigen::VectorXd compute_coeff(const Eigen::MatrixXcd &) const;
};

WeightsBuilder::WeightsBuilder(const int n_lambda, const int n_time,
                               const double radius)
    : lambdas(n_lambda), times(Eigen::VectorXd::LinSpaced(n_time, -1, 1)),
      timestep(2.0 / (n_time - 1)), future_time(1 + timestep) {
  Eigen::VectorXd xs(Eigen::VectorXd::LinSpaced(n_lambda + 1, 0, M_PI + 2));
  for (int i = 0; i < n_lambda; ++i)
    lambdas[i] = radius * semidisk(xs[i]);
}

Eigen::MatrixXcd WeightsBuilder::predictor_matrix() const {
  Eigen::MatrixXcd result(lambdas.size(), 2 * times.size());

  const Eigen::ArrayXXcd b((lambdas * times.transpose()).array().exp());

  result.block(0, 0, lambdas.size(), times.size()) = b;
  result.block(0, times.size(), lambdas.size(), times.size()) =
      b.colwise() * lambdas.array();

  return result;
}

Eigen::MatrixXcd WeightsBuilder::corrector_matrix() const {
  Eigen::MatrixXcd result(lambdas.size(), 2 * times.size() + 1);

  result.block(0, 0, lambdas.size(), 2 * times.size()) = predictor_matrix();
  result.block(0, 2 * times.size(), lambdas.size(), 1) =
      lambdas.array() * rhs_vector().array();

  return result;
}

Eigen::VectorXcd WeightsBuilder::rhs_vector() const {
  Eigen::ArrayXcd b(lambdas * future_time);
  return b.exp();
}

Eigen::VectorXd
WeightsBuilder::compute_coeff(const Eigen::MatrixXcd &mat) const {
  Eigen::JacobiSVD<Eigen::MatrixXcd> decomp =
      mat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

  Eigen::VectorXcd b(rhs_vector()), least_squares(decomp.solve(b));

  return least_squares.real();
}

class Weights {
public:
  Weights(const int, const int, const double);

  Eigen::ArrayXXd ps, cs;
  double future_coef;

  int width() const { return n_time; }

private:
  int n_time;
};

Weights::Weights(const int n_lambda, const int n_time, const double radius)
    : n_time(n_time) {
  const double step_factor = (n_time - 1) / 2.0;
  WeightsBuilder builder(n_lambda, n_time, radius);

  Eigen::VectorXd predictors(builder.predictors());
  Eigen::VectorXd correctors(builder.correctors());

  // Eigen defaults to column major, hence the (n_time x 2) ordering
  ps = Eigen::Map<Eigen::ArrayXXd>(predictors.data(), n_time, 2).transpose();
  cs = Eigen::Map<Eigen::ArrayXXd>(correctors.data(), n_time, 2).transpose();

  ps.row(1) *= step_factor;
  cs.row(1) *= step_factor;
  future_coef = correctors(2 * n_time) * step_factor;
}

extern "C" {
void linear_weights_array(const int *n_lambda, const int *n_time, const double *radius,
                   double *output, double *future_coeff) {
  Weights weights(*n_lambda, *n_time, *radius);

  *future_coeff = weights.future_coef;

  // Linearize the coefficients to the output array
  for (int i = 0; i < *n_time; ++i) {
    const int offset = *n_time;

    // Predictor coefficients
    output[i + 0 * offset] = weights.ps(0, i);
    output[i + 1 * offset] = weights.ps(1, i);

    // Corrector coefficients
    output[i + 2 * offset] = weights.cs(0, i);
    output[i + 3 * offset] = weights.cs(1, i);
  }
}
}