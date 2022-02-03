#ifndef PARAM_INDEXER_H
#define PARAM_INDEXER_H

#include <stddef.h>
#include <vector>
#include <string>

/**
 * class to provide indices for the parameters in the model. If there are K
 * competing risks, with V covariates for each risks and L covariates for each
 * trajectory, the parameters are stored as the
 *
 *  a. V * K parameters for the risk stored consecutively for each cause.
 *  b. L * K parameters for the trajectory stored consecutively fo each cause.
 *  c. the K times K covariance matrix for the random effects.
 */
class param_indexer {
  size_t v_n_cov_risk{0}, v_n_cov_traject{0}, v_n_causes{0};
  static constexpr size_t risk_param_start{0};
  size_t time_param_start{risk_param_start + n_cov_risk() * n_causes()},
         vcov_start{time_param_start + n_cov_traject() * n_causes()};

public:
  param_indexer
    (size_t const n_cov_risk, size_t const n_cov_traject, size_t const n_causes):
  v_n_cov_risk{n_cov_risk}, v_n_cov_traject{n_cov_traject},
  v_n_causes{n_causes} { }

  /// returns the number of causes
  size_t n_causes() const { return v_n_causes; }

  /// returns the number of covariates for each risk
  size_t n_cov_risk() const { return v_n_cov_risk; }

  /// returns the number of covariates for each trajectory
  size_t n_cov_traject() const { return v_n_cov_traject; }

  /// returns the index of covariates for the risks
  size_t risk() const {
    return risk_param_start;
  }

  /**
   * returns the index for the covariates for the given risk.  Cause is zero
   * indexed.
   */
  size_t risk(size_t const cause) const {
    return risk_param_start + cause * n_cov_risk();
  }

  /// returns the index of the covariates for the time-to-event outcomes
  size_t surv() const {
    return time_param_start;
  }

  /**
   * returns the index of the covariates for the time-to-event outcome for the
   * given risk. Cause is zero indexed.
   */
  size_t surv(size_t const cause) const {
    return time_param_start + cause * n_cov_traject();
  }

  /// returns the index for the covariance parameter
  size_t vcov() const {
    return vcov_start;
  }

  /// returns the number of parameters
  template<bool upper_tri>
  size_t n_par() const {
    return upper_tri ? vcov_start + (n_causes() * (n_causes() + 1)) / 2
                     : vcov_start + n_causes() * n_causes();
  }

  /// returns a string with the names for parameters
  template<bool upper_tri>
  std::vector<std::string> param_names() const {
    std::vector<std::string> out;
    out.reserve(n_par<upper_tri>());

    auto add_brackets = [](size_t const i, size_t const j){
      return "[" + std::to_string(i) + "," + std::to_string(j) + "]";
    };

    for(size_t k = 1; k <= n_causes(); ++k)
      for(size_t i = 1; i <= n_cov_risk(); ++i)
        out.emplace_back("beta" + add_brackets(i, k));
    for(size_t k = 1; k <= n_causes(); ++k)
      for(size_t i = 1; i <= n_cov_traject(); ++i)
        out.emplace_back("gamma" + add_brackets(i, k));

    if constexpr (upper_tri){
      size_t counter{};
      for(size_t counter = 0; counter < (n_causes() * (n_causes() + 1)) / 2;)
        out.emplace_back("vcov[" + std::to_string(++counter) + "]");

    } else {
      for(size_t j = 1; j <= n_causes(); ++j)
        for(size_t i = 1; i <= n_causes(); ++i)
          out.emplace_back("vcov" + add_brackets(i, j));

    }

    return out;
  }
};

#endif
