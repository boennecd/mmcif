#include <testthat.h>
#include "mmcif-logLik.h"

namespace {
constexpr double ghq_nodes[]{-3.43615911883774, -2.53273167423279, -1.75668364929988, -1.03661082978951, -0.342901327223705, 0.342901327223705, 1.03661082978951, 1.75668364929988, 2.53273167423279, 3.43615911883774},
               ghq_weights[]{7.6404328552326e-06, 0.00134364574678124, 0.0338743944554811, 0.240138611082314, 0.610862633735326, 0.610862633735326, 0.240138611082315, 0.033874394455481, 0.00134364574678124, 7.64043285523265e-06};

const ghqCpp::ghq_data ghq_dat_use{ghq_nodes, ghq_weights, 10};

/*
# util functions
 upper_to_full <- \(x){
 dim <- (sqrt(8 * length(x) + 1) - 1) / 2
 out <- matrix(0, dim, dim)
 out[upper.tri(out, TRUE)] <- x
 out[lower.tri(out)] <- t(out)[lower.tri(out)]
 out
 }
 get_n_remove <- \(x, n){
 out <- x[1:n]
 eval(substitute(out <- out[-(1:n)], list(out = substitute(x), n = n)),
 parent.frame())
 out
 }
 d_upper_to_full <- \(x){
 dim <- (sqrt(8 * length(x) + 1) - 1) / 2
 out <- matrix(0, dim, dim)
 out[upper.tri(out, TRUE)] <- x
 out[upper.tri(out)] <- out[upper.tri(out)] / 2
 out[lower.tri(out)] <- t(out)[lower.tri(out)]
 out
 }
 */

} // namespace

context("mmcif_logLik works as expected with singleton data") {
  test_that("mmcif_logLik works") {
    /*
     set.seed(321)
     n_cov_traject <- 4L
     n_cov_risk <- 2L
     n_causes <- 3L

     g <- \(x) cbind(x, log(x))
     d_g <- \(x) cbind(1, 1/x)

     dput(covs_risk <- runif(n_cov_risk, -1) |> round(3))
     covs_traject <- runif(n_cov_traject, -1) |> round(3)
     vcov <-
     drop(rWishart(1, 2 * n_causes, diag(1/2/n_causes, 2 * n_causes))) |>
     drop() |> round(3)

     coefs_risk <- runif(n_causes * n_cov_risk, -1) |>
     round(3) |> matrix(n_cov_risk)
     coefs_traject <- runif(n_causes * n_cov_traject, -1) |>
     round(3) |> matrix(n_cov_traject)
     coefs_traject <- rbind(-.5, -.25, coefs_traject)
     c(coefs_risk, coefs_traject, vcov) |> dput()

# sample a plausible value
     u <- mvtnorm::rmvnorm(1, sigma = vcov) |> drop()
     etas <- covs_risk %*% coefs_risk |> drop()
     etas <- etas + u[1:n_causes]
     phats <- c(exp(etas), 1) / (1 + sum(exp(etas)))
     dput(cause <- sample(n_causes + 1, 1, prob = phats))
     stopifnot(cause < n_causes + 1) # we did not implement this case

     g_w_covs <- \(ti)
     cbind(g(ti), rep(covs_traject, each = length(ti)) |> matrix(length(ti)))
     rng <- runif(1)
     obs_time <- uniroot(
     \(ti){
     lp <- -g_w_covs(ti) %*% coefs_traject[, cause] - u[cause + n_causes]
     pnorm(lp) - rng
     },
     c(1e-16, 1000))$root

     d_g_w_covs <- \(ti)
     cbind(d_g(ti), matrix(0., length(ti), length(covs_traject)))

# the matrices we will need
     dput(covs_traject_w_time <- g_w_covs(obs_time) |> drop())
     dput(d_covs_traject_w_time <- d_g_w_covs(obs_time) |> drop())

     dput(gl <- fastGHQuad::gaussHermiteData(10L))

# compute the log densities we need. Start with the censored case
# trajectory probability is one
     library(ghqCpp)
     f <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |> matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)

     etas <- covs_risk %*% coefs_risk |> drop()
     res <- mixed_mult_logit_term(
     eta = as.matrix(etas), Sigma = vcov[1:n_causes, 1:n_causes],
     which_category = 0L, weights = gl$w, nodes = gl$x)
     log(res)
     }

     par <- c(coefs_risk, coefs_traject, vcov[upper.tri(vcov, TRUE)])
     f(par) |> dput()

     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     gr <- c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov)))
     dput(gr + 1.5)

# then with the observed outcome
     f <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |> matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)

     idx_cause <- cause + n_causes

     res <- log(-d_covs_traject_w_time %*% coefs_traject[, cause]) |> drop()
     res <- res + dnorm(-covs_traject_w_time %*% coefs_traject[, cause] |> drop(),
     sd = sqrt(1 + vcov[idx_cause, idx_cause]), log = TRUE)

     M <- solve(vcov)
     M[idx_cause, idx_cause] <- M[idx_cause, idx_cause] + 1
     M <- solve(M)

     mean_cond <- M[, idx_cause] *
     sum(-covs_traject_w_time %*% coefs_traject[, cause])
     etas <- covs_risk %*% coefs_risk |> drop()
     etas <- etas + mean_cond[1:n_causes]

     res_integral <- mixed_mult_logit_term(
     eta = as.matrix(etas), Sigma = M[1:n_causes, 1:n_causes],
     which_category = cause, weights = gl$w, nodes = gl$x)

     res + log(res_integral)
     }

     f(par) |> dput()

     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     gr <- c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov)))
     dput(gr + .5)

# the censored case
     f <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |> matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)

     integrand <- 1
     etas <- covs_risk %*% coefs_risk |> drop()
     vcov_sub <- vcov[1:n_causes, 1:n_causes]

     for(cause in 1:n_causes){
     idx_cause <- cause + n_causes
     z <- -solve(vcov_sub, vcov[1:n_causes, idx_cause]) |> drop()
     s <- sqrt(1 + vcov[idx_cause, idx_cause] -
     sum(vcov[idx_cause, 1:n_causes] * (-z)))

     offset <- -covs_traject_w_time %*% coefs_traject[, cause] |> drop()

     integrand <- integrand - mixed_mult_logit_n_probit_term(
     eta = as.matrix(etas),
     which_category = cause, s = s, eta_probit = offset,
     Sigma = vcov_sub, z = z, weights = gl$w, nodes = gl$x)
     }
     log(integrand)
     }

     f(par) |> dput()

     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     gr <- c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov)))
     dput(gr + .5)
     */

    constexpr size_t n_causes{3}, n_cov_risk{2}, n_cov_traject{6};
    constexpr double covs_risk[]{0.912, 0.875},
                  covs_traject[]{4.1936419665673, 1.43356956082597, -0.524, -0.49, -0.219, -0.318},
                d_covs_traject[]{1, 0.238456217286129, 0, 0, 0, 0},
                           par[]{0.355, -0.737, -0.475, -0.365, 0.353, -0.874, -0.5, -0.25, -0.871, 0.322, -0.791, 0.474, -0.5, -0.25, -0.85, 0.782, 0.78, -0.123, -0.5, -0.25, -0.81, 0.091, -0.198, -0.974, 0.771, -0.126, -0.198, 0.068, 0.879, -0.145, -0.126, 1.158, 0.184, 1.053, -0.211, 0.595, -0.198, 0.184, 1.456, -0.43, -1.171, -0.332, 0.068, 1.053, -0.43, 2.234, 0.431, 0.462, 0.879, -0.211, -1.171, 0.431, 1.673, -0.001, -0.145, 0.595, -0.332, 0.462, -0.001, 0.832};

    ghqCpp::simple_mem_stack<double> mem;
    param_indexer indexer{n_cov_risk, n_cov_traject, n_causes};

    {
      // the censored case with one as the probability of the trajectory
      mmcif_data dat{covs_traject, d_covs_traject, covs_risk, false, n_causes};
      double res{mmcif_logLik(par, indexer, dat, mem, ghq_dat_use)};
      constexpr double truth{-1.21423993502317};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      double * gr{mem.get(indexer.n_par<false>())};
      constexpr double shift{1.5},
                   true_gr[]{1.26365754060744, 1.2732459955188, 1.34146548646436, 1.34789725971731, 1.30061198121401, 1.30870118784666, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.46383170621554, 1.53592172902375, 1.54244818397882, 1.5, 1.5, 1.5, 1.53592172902375, 1.4623397748844, 1.53217128260594, 1.5, 1.5, 1.5, 1.54244818397882, 1.53217128260594, 1.46703163289562, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

      std::fill(gr, gr + indexer.n_par<false>(), shift);
      res = mmcif_logLik_grad(par, gr, indexer, dat, mem, ghq_dat_use);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      for(size_t i = 0; i < indexer.n_par<false>(); ++i)
        expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-6);
    }

    {
      // the observed case works
      constexpr unsigned cause{2L};
      mmcif_data dat{covs_traject, d_covs_traject, covs_risk, true, cause};
      double res{mmcif_logLik(par, indexer, dat, mem, ghq_dat_use)};
      constexpr double truth{-4.26186714415087};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      double * gr{mem.get(indexer.n_par<false>())};
      constexpr double shift{.5},
                   true_gr[]{0.351552206482173, 0.357574759929686, 0.282709259050858, 0.2915247825613, 1.10517471974935, 1.0806226750331, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 3.43029581567022, 1.68645040335615, -0.0894243337010308, -0.0511792441427352, 0.253656624124913, 0.142295919574298, 0.461770301601013, 0.532975866589409, 0.45954237109379, 0.500000003534147, 0.500000000567034, 0.402421270586787, 0.532975866589409, 0.460900127868904, 0.447644291654297, 0.499999999797272, 0.499999999555276, 0.393827579249467, 0.45954237109379, 0.447644291654297, 0.648183312792179, 0.500000000362377, 0.499999999926135, 0.852567814444038, 0.500000003534147, 0.499999999797272, 0.500000000362377, 0.499999999996471, 0.500000000000939, 0.499999999532829, 0.500000000567034, 0.499999999555276, 0.499999999926135, 0.500000000000939, 0.499999999824738, 0.5, 0.402421270586787, 0.393827579249467, 0.852567814444038, 0.499999999532829, 0.5, 0.846467447468548};

      std::fill(gr, gr + indexer.n_par<false>(), shift);
      res = mmcif_logLik_grad(par, gr, indexer, dat, mem, ghq_dat_use);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      for(size_t i = 0; i < indexer.n_par<false>(); ++i)
        expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
    }

    // the censored case with a probability of the trajectory in (0, 1) works
    mmcif_data dat{covs_traject, d_covs_traject, covs_risk, true, n_causes};
    double res{mmcif_logLik(par, indexer, dat, mem, ghq_dat_use)};
    constexpr double truth{-1.04608338904514};
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    double * gr{mem.get(indexer.n_par<false>())};
    constexpr double shift{.5},
                 true_gr[]{0.330648438629487, 0.33751906083885, 0.350726364390927, 0.356782421714224, 0.352693743921712, 0.358669984812303, 0.833300440713479, 0.613936613988158, 0.458353757359891, 0.461055994306309, 0.482594413812261, 0.47472613513888, 0.638228487237129, 0.547252520431373, 0.482728204306005, 0.48384889336008, 0.492781444202289, 0.489518261901478, 0.800680721846452, 0.60278577298994, 0.462429625797045, 0.464867397951781, 0.484297878241694, 0.47719965841653, 0.477058826417935, 0.524688031342338, 0.532398053180518, 0.524647165965919, 0.493874870084439, 0.494164701133886, 0.524688031342338, 0.468439391298966, 0.525877057569041, 0.491270491704209, 0.511791124947451, 0.491458569479165, 0.532398053180518, 0.525877057569041, 0.477259089716002, 0.494835460522858, 0.498832652888146, 0.523788667244225, 0.524647165965919, 0.491270491704209, 0.494835460522858, 0.527860742985394, 0.5, 0.5, 0.493874870084439, 0.511791124947451, 0.498832652888146, 0.5, 0.517999932390265, 0.5, 0.494164701133886, 0.491458569479165, 0.523788667244225, 0.5, 0.5, 0.540325615606095};

    std::fill(gr, gr + indexer.n_par<false>(), shift);
    res = mmcif_logLik_grad(par, gr, indexer, dat, mem, ghq_dat_use);
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    for(size_t i = 0; i < indexer.n_par<false>(); ++i)
      expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
  }
}

namespace {
/*
 set.seed(321)
 n_cov_traject <- 4L
 n_cov_risk <- 2L
 n_causes <- 3L

 g <- \(x) cbind(x, log(x))
 d_g <- \(x) cbind(1, 1/x)

 vcov <-
 drop(rWishart(1, 2 * n_causes, diag(1/2/n_causes, 2 * n_causes))) |>
 drop() |> round(3)

 coefs_risk <- runif(n_causes * n_cov_risk, -1) |>
 round(3) |> matrix(n_cov_risk)
 coefs_traject <- runif(n_causes * n_cov_traject, -1) |>
 round(3) |> matrix(n_cov_traject)
 coefs_traject <- rbind(-.5, -.25, coefs_traject)
 c(coefs_risk, coefs_traject, vcov) |> dput()

# sample a plausible value
 u <- mvtnorm::rmvnorm(1, sigma = vcov) |> drop()

 g_w_covs <- \(ti, covs_traject)
 cbind(g(ti), rep(covs_traject, each = length(ti)) |> matrix(length(ti)))
 d_g_w_covs <- \(ti, covs_traject)
 cbind(d_g(ti), matrix(0., length(ti), length(covs_traject)))

 sample_obs <- \(){
 covs_risk <- runif(n_cov_risk, -1) |> round(3)
 covs_traject <- runif(n_cov_traject, -1) |> round(3)

 etas <- covs_risk %*% coefs_risk |> drop()
 etas <- etas + u[1:n_causes]
 phats <- c(exp(etas), 1) / (1 + sum(exp(etas)))
 cause <- sample(n_causes + 1, 1, prob = phats)
 stopifnot(cause < n_causes + 1) # we did not implement this case

 rng <- runif(1)
 obs_time <- uniroot(
 \(ti){
 lp <- -g_w_covs(ti, covs_traject) %*% coefs_traject[, cause] -
 u[cause + n_causes]
 pnorm(lp) - rng
 },
 c(1e-16, 1000))$root

 list(covs_risk = covs_risk,
 covs_traject_w_time = g_w_covs(obs_time, covs_traject) |> drop(),
 d_covs_traject_w_time = d_g_w_covs(obs_time, covs_traject) |> drop(),
 cause = cause)
 }

 dput(obs <- replicate(2, sample_obs(), simplify = FALSE))

 dput(gl <- fastGHQuad::gaussHermiteData(10L))
 */

constexpr size_t n_cov_risk{2},
                 n_cov_traject{6},
                 n_causes{3},
                 cause[]{0, 2};

param_indexer const indexer{n_cov_risk, n_cov_traject, n_causes};

constexpr double covs_risk1[]{-0.605, -0.876},
                 covs_risk2[]{0.584, 0.576},
              covs_traject1[]{0.85860256800553, -0.152449132274265, -0.734, 0.062, 0.957, 0.198},
              covs_traject2[]{4.78877056139744, 1.5662737107113, 0.304, 0.267, -0.669, 0.204},
            d_covs_traject1[]{1, 1.16468321580138, 0, 0, 0, 0},
            d_covs_traject2[]{1, 0.208821865065129, 0, 0, 0, 0},
                        par[]{-0.647, 0.076, -0.439, -0.98, -0.806, 0.324, -0.5, -0.25, 0.355, -0.737, -0.475, -0.365, -0.5, -0.25, 0.353, -0.874, -0.871, 0.322, -0.5, -0.25, -0.791, 0.474, -0.85, 0.782, 1.974, -0.235, 0.154, 0.195, 0.108, 0.059, -0.235, 0.401, 0.163, -0.161, 0.597, -0.155, 0.154, 0.163, 0.4, 0.026, 0.047, 0.179, 0.195, -0.161, 0.026, 0.506, -0.493, 0.018, 0.108, 0.597, 0.047, -0.493, 2.031, 0.052, 0.059, -0.155, 0.179, 0.018, 0.052, 0.578};
} // namespace

context("mmcif_logLik works as expected with bivariate data") {
  test_that("mmcif_logLik works when both individuals are observed") {
    /*
     library(ghqCpp)
     f <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)

     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     covs_traject <- sapply(
     obs,
     \(x, cause) -x$covs_traject_w_time %*% coefs_traject[, x$cause])
     d_covs_traject <- sapply(
     obs,
     \(x) -x$d_covs_traject_w_time %*% coefs_traject[, x$cause])

     V <- matrix(0, 2, 2 * n_causes)
     V[1, obs[[1]]$cause + n_causes] <- 1
     V[2, obs[[2]]$cause + n_causes] <- 1
     out <- sum(log(d_covs_traject)) +
     mvtnorm::dmvnorm(
     covs_traject, sigma = diag(2) + V %*% tcrossprod(vcov, V), log = TRUE)

     vcov_cond <- solve(crossprod(V) + solve(vcov))
     mean_cond <- vcov_cond %*% crossprod(V, covs_traject) |> drop()

     integral <- mixed_mult_logit_term(
     eta = etas + mean_cond[1:n_causes],
     Sigma = vcov_cond[1:n_causes, 1:n_causes],
     which_category = sapply(obs, `[[`, "cause"),
     weights = gl$w, nodes = gl$x)

     out + log(integral)
     }

     par <- c(coefs_risk, coefs_traject, vcov[upper.tri(vcov, TRUE)])
     f(par) |> dput()
     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     (c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov))) - 4) |> dput()
     */
    ghqCpp::simple_mem_stack<double> mem;

    mmcif_data
      obs1{covs_traject1, d_covs_traject1, covs_risk1, true, cause[0]},
      obs2{covs_traject2, d_covs_traject2, covs_risk2, true, cause[1]};

    double res{mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
    constexpr double truth{-7.77301433778719};
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    double * gr{mem.get(indexer.n_par<false>())};
    constexpr double shift{-4},
                 true_gr[]{-4.60044096320519, -4.76143069955199, -3.88925107716226, -3.81874883303431, -3.44226495101285, -3.38684934060633, -4.63995658752879, -5.58289389453021, -4.53343755123284, -3.95494124232041, -3.30449627176941, -3.85610267649002, -4, -4, -4, -4, -4, -4, 0.276190040738317, -2.38723521835543, -3.61357900292289, -3.66061050597717, -4.85038041763269, -3.7406911732045, -4.08660377973573, -4.01372576629899, -3.9001876680933, -3.91005332184739, -4.00000000811773, -3.8691233090845, -4.01372576629899, -4.04876332825872, -4.03306286597394, -4.13643808454184, -3.99999997350787, -4.22974199345051, -3.9001876680933, -4.03306286597394, -4.02777117647809, -3.81863162075935, -4.00000061380221, -3.66503522108475, -3.91005332184739, -4.13643808454184, -3.81863162075935, -4.0715034184455, -4.00000001334417, -3.53643481460751, -4.00000000811773, -3.99999997350787, -4.00000061380221, -4.00000001334417, -4.0000000095682, -3.99999984574, -3.8691233090845, -4.22974199345051, -3.66503522108475, -3.53643481460751, -3.99999984574, -3.51292697279014};

    std::fill(gr, gr + indexer.n_par<false>(), shift);
    res = mmcif_logLik_grad(par, gr, indexer, obs1, obs2, mem, ghq_dat_use);
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    for(size_t i = 0; i < indexer.n_par<false>(); ++i)
      expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
  }

  test_that("mmcif_logLik works when one individual is observed and one is censored") {
    /*
     library(ghqCpp)
     f <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)
     vcov <- (vcov + t(vcov)) / 2

     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     covs_traject <- sapply(
     obs,
     \(x, cause) -x$covs_traject_w_time %*% coefs_traject[, x$cause])
     d_covs_traject <- sapply(
     obs,
     \(x) -x$d_covs_traject_w_time %*% coefs_traject[, x$cause])

     v <- numeric(2 * n_causes)
     v[n_causes + obs[[1]]$cause] <- 1
     out <- log(d_covs_traject[1]) + dnorm(
     covs_traject[1], sd = sqrt(1 + v %*% vcov %*% v), log = TRUE)

     vcov_cond <- solve(tcrossprod(v) + solve(vcov))
     mean_cond <- drop(vcov_cond %*% v) * covs_traject[1]

     integral <- mixed_mult_logit_term(
     eta = etas + mean_cond[1:n_causes],
     Sigma = vcov_cond[1:n_causes, 1:n_causes],
     which_category = c(obs[[1]]$cause, 0L),
     weights = gl$w, nodes = gl$x)

     out + log(integral)
     }

     par <- c(coefs_risk, coefs_traject, vcov[upper.tri(vcov, TRUE)])
     f(par) |> dput()
     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     (c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov))) - 4) |> dput()
     */

    ghqCpp::simple_mem_stack<double> mem;

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, true, cause[0]},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, false, n_causes};

      double res{mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
      constexpr double truth{-4.5240025456571};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      double * gr{mem.get(indexer.n_par<false>())};
      constexpr double shift{-4},
                   true_gr[]{-4.60677825999459, -4.76876219981359, -3.86894022626248, -3.78371293952945, -4.02220331100695, -3.97886624482647, -4.62087217947148, -5.58628242493471, -4.54975238220723, -3.95356315005794, -3.28322475503409, -3.85170167375117, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4.08651843041873, -4.01607930345544, -3.99983074431344, -3.90823923485969, -3.99999986891621, -4.00000003552164, -4.01607930345544, -4.03105028895141, -3.88643373484883, -4.16983704586129, -4.00000001555981, -3.9999999372787, -3.99983074431344, -3.88643373484883, -4.05931545957518, -4.12982271375971, -3.99999956557467, -4.00000006439723, -3.90823923485969, -4.16983704586129, -4.12982271375971, -4.05536084376252, -4.00000000252952, -3.99999957702496, -3.99999986891621, -4.00000001555981, -3.99999956557467, -4.00000000252952, -3.99999999300296, -4.00000024008675, -4.00000003552164, -3.9999999372787, -4.00000006439723, -3.99999957702496, -4.00000024008675, -4.00000001062483};

      std::fill(gr, gr + indexer.n_par<false>(), shift);
      res = mmcif_logLik_grad(par, gr, indexer, obs1, obs2, mem, ghq_dat_use);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      for(size_t i = 0; i < indexer.n_par<false>(); ++i)
        expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
    }

    /*
     f <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)
     vcov <- (vcov + t(vcov)) / 2

     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     covs_traject <- sapply(
     obs,
     \(x, cause) -x$covs_traject_w_time %*% coefs_traject[, x$cause])
     d_covs_traject <- sapply(
     obs,
     \(x) -x$d_covs_traject_w_time %*% coefs_traject[, x$cause])

     v <- numeric(2 * n_causes)
     v[n_causes + obs[[1]]$cause] <- 1
     out <- log(d_covs_traject[1]) + dnorm(
     covs_traject[1], sd = sqrt(1 + v %*% vcov %*% v), log = TRUE)

     vcov_cond <- solve(tcrossprod(v) + solve(vcov))
     mean_cond <- drop(vcov_cond %*% v) * covs_traject[1]

     integral <- mixed_mult_logit_term(
     eta = as.matrix(etas[, 1]) + mean_cond[1:n_causes],
     Sigma = vcov_cond[1:n_causes, 1:n_causes],
     which_category = obs[[1]]$cause, weights = gl$w, nodes = gl$x)

     integrals <- sapply(1:n_causes, \(cause){
     idx_cause <- cause + n_causes
     rng_coefs <- solve(vcov_cond[1:n_causes, 1:n_causes],
     vcov_cond[1:n_causes, idx_cause])

     s <- sqrt(1 + vcov_cond[idx_cause, idx_cause] -
     vcov_cond[idx_cause, 1:n_causes] %*% rng_coefs)
     shift <- -obs[[2]]$covs_traject_w_time %*% coefs_traject[, cause] -
     mean_cond[idx_cause]

     rng_coefs <- -rng_coefs

     mixed_mult_logit_n_probit_term(
     eta = etas + mean_cond[1:n_causes],
     which_category = c(obs[[1]]$cause, cause), s = s, eta_probit = shift,
     Sigma = vcov_cond[1:n_causes, 1:n_causes], z = rng_coefs,
     weights = gl$w, nodes = gl$x)
     })
     integral <- integral - sum(integrals)

     out + log(integral)
     }

     par <- c(coefs_risk, coefs_traject, vcov[upper.tri(vcov, TRUE)])
     f(par) |> dput()
     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     (c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov))) - 4) |> dput()
     */

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, true, cause[0]},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, true, n_causes};

      constexpr double truth{-4.38551973954109},
                       shift{-4},
                   true_gr[]{-4.5671761432415, -4.72543744207062, -3.85884139924293, -3.77625459152193, -3.9998547663936, -3.95693855135666, -4.034641265227, -5.37854908528959, -4.48343114496485, -3.9213991068482, -3.40487582437673, -3.83268709528386, -3.87120856084756, -3.95787590094769, -3.99182408151401, -3.99281917690393, -4.01799239946641, -3.99451352851197, -3.70694231274736, -3.90414897823354, -3.98139615681146, -3.98366044018079, -4.04094069464688, -3.9875158420464, -4.06985265337591, -4.01580871407255, -3.99979185548028, -3.86305499015696, -3.99623191508629, -3.99371522788762, -4.01580871407255, -4.03410448588676, -3.89887558953613, -4.16294785134668, -3.99471714867857, -4.01026172632862, -3.99979185548028, -3.89887558953613, -4.05858870479936, -4.12763313213162, -4.0044183991572, -3.98484286950757, -3.86305499015696, -4.16294785134668, -4.12763313213162, -3.94364666714109, -3.98520159014463, -3.97776221043803, -3.99623191508629, -3.99471714867857, -4.0044183991572, -3.98520159014463, -3.98865659435786, -4.00000027578256, -3.99371522788762, -4.01026172632862, -3.98484286950757, -3.97776221043803, -4.00000027578256, -3.96110568888788};
      {
        double res
          {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
        expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

        double * gr{mem.get(indexer.n_par<false>())};
        std::fill(gr, gr + indexer.n_par<false>(), shift);
        res = mmcif_logLik_grad(par, gr, indexer, obs1, obs2, mem, ghq_dat_use);
        expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

        for(size_t i = 0; i < indexer.n_par<false>(); ++i)
          expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
      }

      double res
        {mmcif_logLik(par, indexer, obs2, obs1, mem, ghq_dat_use)};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      double * gr{mem.get(indexer.n_par<false>())};
      std::fill(gr, gr + indexer.n_par<false>(), shift);
      res = mmcif_logLik_grad(par, gr, indexer, obs2, obs1, mem, ghq_dat_use);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      for(size_t i = 0; i < indexer.n_par<false>(); ++i)
        expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
    }
  }

  test_that("mmcif_logLik works when both individuals are censored") {
    /*
     f <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)
     vcov <- (vcov + t(vcov)) / 2

     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     mixed_mult_logit_term(
     eta = etas, Sigma = vcov[1:n_causes, 1:n_causes],
     which_category = integer(2), weights = gl$w, nodes = gl$x) |>
     log()
     }

     par <- c(coefs_risk, coefs_traject, vcov[upper.tri(vcov, TRUE)])
     f(par) |> dput()
     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     (c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov))) - 4) |> dput()
     */

    ghqCpp::simple_mem_stack<double> mem;

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, false, n_causes},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, false, n_causes};

      double res
        {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
      constexpr double truth{-3.01479043790748};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      double * gr{mem.get(indexer.n_par<false>())};
      constexpr double shift{-4},
                   true_gr[]{-4.00600061167854, -3.95429338506848, -3.81463703103875, -3.68922104267058, -4.03296854671522, -3.98498043300022, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4.00651859808217, -3.86236765460961, -3.89687944480354, -4, -4, -4, -3.86236765460961, -3.95860579328351, -3.82086704477511, -4, -4, -4, -3.89687944480354, -3.82086704477511, -4.05656998194764, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4};

      std::fill(gr, gr + indexer.n_par<false>(), shift);
      res = mmcif_logLik_grad(par, gr, indexer, obs1, obs2, mem, ghq_dat_use);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      for(size_t i = 0; i < indexer.n_par<false>(); ++i)
        expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
    }

    /*
     f <- function(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)
     vcov <- (vcov + t(vcov)) / 2

     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     integral <- mixed_mult_logit_term(
     eta = as.matrix(etas[, 2]), Sigma = vcov[1:n_causes, 1:n_causes],
     which_category = 0L, weights = gl$w, nodes = gl$x)

     integrals <- sapply(1:n_causes, \(cause){
     lp_traject <- -obs[[1]]$covs_traject_w_time %*% coefs_traject[, cause] |>
     drop()
     idx_cause <- cause + n_causes
     rng_coefs <-
     solve(vcov[1:n_causes, 1:n_causes], vcov[1:n_causes, idx_cause])
     s <- sqrt(1 + vcov[idx_cause, idx_cause] -
     vcov[idx_cause, 1:n_causes] %*% rng_coefs)
     rng_coefs <- -rng_coefs
     mixed_mult_logit_n_probit_term(
     eta = etas, which_category = c(cause, 0L), s = s,
     eta_probit = lp_traject, Sigma = vcov[1:n_causes, 1:n_causes],
     z = rng_coefs, weights = gl$w, nodes = gl$x)
     })
     log(integral - sum(integrals))
     }

     par <- c(coefs_risk, coefs_traject, vcov[upper.tri(vcov, TRUE)])
     f(par) |> dput()
     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     (c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov))) + 1.5) |> dput()
     */

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, true, n_causes},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, false, n_causes};

      constexpr double truth{-2.18417107496657},
                       shift{1.5},
                   true_gr[]{1.42886432212848, 1.45568701354807, 1.52143756531699, 1.57829984536793, 1.34545167726186, 1.34168374815065, 1.60455529481686, 1.4814356903032, 1.41061803300473, 1.50754997550706, 1.61653752380935, 1.52411121171035, 1.68542942744349, 1.46707608763835, 1.34148055836095, 1.51338992552404, 1.70667998038838, 1.54276137510894, 1.62477679072493, 1.47784526361077, 1.39333113149789, 1.50901017690837, 1.63907644033356, 1.52877443596108, 1.48077385432297, 1.58360118439227, 1.5560064124667, 1.51232498904552, 1.46513088846104, 1.47203468864447, 1.58360118439227, 1.47288847995647, 1.58987619120269, 1.47503865108853, 1.52039508343806, 1.45664598170391, 1.5560064124667, 1.58987619120269, 1.43986715799987, 1.47858474718355, 1.4588547390977, 1.53606382712226, 1.51232498904552, 1.47503865108853, 1.47858474718355, 1.54560315711551, 1.5, 1.5, 1.46513088846104, 1.52039508343806, 1.4588547390977, 1.5, 1.55038611784376, 1.5, 1.47203468864447, 1.45664598170391, 1.53606382712226, 1.5, 1.5, 1.51294882404452};
      {
        double res
          {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
        expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

        double * gr{mem.get(indexer.n_par<false>())};
        std::fill(gr, gr + indexer.n_par<false>(), shift);
        res = mmcif_logLik_grad(par, gr, indexer, obs1, obs2, mem, ghq_dat_use);
        expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

        for(size_t i = 0; i < indexer.n_par<false>(); ++i)
          expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
      }

      double res
        {mmcif_logLik(par, indexer, obs2, obs1, mem, ghq_dat_use)};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      double * gr{mem.get(indexer.n_par<false>())};
      std::fill(gr, gr + indexer.n_par<false>(), shift);
      res = mmcif_logLik_grad(par, gr, indexer, obs2, obs1, mem, ghq_dat_use);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      for(size_t i = 0; i < indexer.n_par<false>(); ++i)
        expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
    }

    /*
     library(ghqCpp)
     f <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(NROW(coefs_risk))
     coefs_traject <- get_n_remove(x, length(coefs_traject)) |>
     matrix(NROW(coefs_traject))
     vcov <- upper_to_full(x)
     vcov <- (vcov + t(vcov)) / 2

     vcov_sub <- vcov[1:n_causes, 1:n_causes]

     integrals_substract <- sapply(obs, \(obs_i){
     sapply(1:n_causes, \(cause){
     etas <- obs_i$covs_risk %*% coefs_risk |> drop()

     idx_cause <- cause + n_causes
     z <- -solve(vcov_sub, vcov[1:n_causes, idx_cause]) |> drop()
     s <- sqrt(1 + vcov[idx_cause, idx_cause] -
     sum(vcov[idx_cause, 1:n_causes] * (-z)))

     offset <- -obs_i$covs_traject_w_time %*% coefs_traject[, cause] |>
     drop()

     mixed_mult_logit_n_probit_term(
     eta = as.matrix(etas),
     which_category = cause, s = s, eta_probit = offset,
     Sigma = vcov_sub, z = z, weights = gl$w, nodes = gl$x)
     })
     })

     vcov_sub <- vcov[1:n_causes, 1:n_causes]

     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     rng_coefs <- solve(vcov_sub, vcov[1:n_causes, -(1:n_causes)]) |> t()
     vcov_cond <- vcov[1:n_causes + n_causes, 1:n_causes + n_causes] -
     vcov[1:n_causes + n_causes, 1:n_causes] %*% t(rng_coefs)

     integrals_adds <- sapply(
     1:n_causes,
     \(cause1) sapply(
     1:n_causes, \(cause2){
     V <- matrix(0, 2, n_causes)
     V[1, cause1] <- 1
     V[2, cause2] <- 1

     shift <- -c(obs[[1]]$covs_traject_w_time %*% coefs_traject[, cause1],
     obs[[2]]$covs_traject_w_time %*% coefs_traject[, cause2])

     rng_coefs <- -V %*% rng_coefs
     vcov_cond_pbvn <- diag(2) + V %*% tcrossprod(vcov_cond, V)

     mixed_mult_logit_n_cond_pbvn(
     eta = etas, which_category = c(cause1, cause2),
     eta_pbvn = -shift, Psi = vcov_cond_pbvn, V = -rng_coefs,
     Sigma = vcov_sub, weights = gl$w, nodes = gl$x)
     }))

     log(1 - sum(integrals_substract) + sum(integrals_adds))
     }

     par <- c(coefs_risk, coefs_traject, vcov[upper.tri(vcov, TRUE)])
     f(par) |> dput()
     gr <- numDeriv::grad(f, par)
     dim_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
     (c(head(gr, -dim_vcov), d_upper_to_full(tail(gr, dim_vcov))) + 1.5) |> dput()
     */

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, true, n_causes},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, true, n_causes};

      double res
        {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
      constexpr double truth{-2.00247380641032};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      double * gr{mem.get(indexer.n_par<false>())};
      constexpr double shift{1.5},
                   true_gr[]{1.4392362543624, 1.4657164371737, 1.55664825291928, 1.60435190353244, 1.37028357120591, 1.36660932696277, 1.77073755216239, 1.53802079510158, 1.42522865934166, 1.51674277488534, 1.5877749109731, 1.53035913222007, 2.06304090224125, 1.59159588972652, 1.36729663448639, 1.5344107677658, 1.65140759199596, 1.5584698669516, 1.94148079104312, 1.58523847207786, 1.4203664235269, 1.52654397141941, 1.58536679406214, 1.54084713201152, 1.48091073070106, 1.56488498397801, 1.54783934446057, 1.52202302260823, 1.4560063687567, 1.46367963873689, 1.56488498397801, 1.49314850571281, 1.56749376537132, 1.47422099404467, 1.56126698068742, 1.45113517043605, 1.54783934446057, 1.56749376537132, 1.45156009955787, 1.47872440662632, 1.45122709058752, 1.56317435473337, 1.52202302260823, 1.47422099404467, 1.47872440662632, 1.58411611553288, 1.50274810184889, 1.50547469903672, 1.4560063687567, 1.56126698068742, 1.45122709058752, 1.50274810184889, 1.58144204720316, 1.50857653786166, 1.46367963873689, 1.45113517043605, 1.56317435473337, 1.50547469903672, 1.50857653786166, 1.5622719445094};

      std::fill(gr, gr + indexer.n_par<false>(), shift);
      res = mmcif_logLik_grad(par, gr, indexer, obs1, obs2, mem, ghq_dat_use);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      for(size_t i = 0; i < indexer.n_par<false>(); ++i)
        expect_true(std::abs(gr[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-3);
    }
  }
}
