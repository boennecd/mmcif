
We load and prepare the prostate cancer data set.

``` r
library(mets)
#> Loading required package: timereg
#> Loading required package: survival
#> Loading required package: lava
#> mets version 1.2.9
library(mmcif)

data(prt)
str(prt)
#> 'data.frame':    29222 obs. of  6 variables:
#>  $ country: Factor w/ 4 levels "Denmark","Finland",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ time   : num  97 80.9 68 61.5 78.8 ...
#>  $ status : int  1 1 1 1 1 1 1 1 1 2 ...
#>  $ zyg    : Factor w/ 2 levels "DZ","MZ": 1 1 1 1 1 1 2 2 1 1 ...
#>  $ id     : int  1 1 3 3 5 5 9 9 12 12 ...
#>  $ cancer : num  0 0 0 0 0 0 0 0 0 1 ...

with(prt, quantile(time[status > 0], probs = seq(0, 1, by = .1)))
#>       0%      10%      20%      30%      40%      50%      60%      70% 
#>   9.9212  53.6796  62.8195  68.1046  72.4418  76.4447  80.2888  83.9669 
#>      80%      90%     100% 
#>  87.8642  93.0579 117.6224
max_time <- 90
prt <- within(prt, {
  status[time >= max_time] <- 0
  time <- pmin(time, max_time)
})

# use another code for censoring
prt_dz <- subset(prt, zyg == "DZ") |> 
  transform(status = ifelse(status == 0, 3L, status))
prt_mz <- subset(prt, zyg == "MZ") |> 
  transform(status = ifelse(status == 0, 3L, status))

# check distribution of observed events disregarding time
table(prt_dz$status)
#> 
#>     1     2     3 
#>  3916   522 13553
xtabs(~ status + country, prt_dz)
#>       country
#> status Denmark Finland Norway Sweden
#>      1    1303     792    482   1339
#>      2      80     118     69    255
#>      3    4808    1923   1842   4980
prop.table(xtabs(~ status + country, prt_dz), margin = 2L)
#>       country
#> status  Denmark  Finland   Norway   Sweden
#>      1 0.210467 0.279562 0.201421 0.203681
#>      2 0.012922 0.041652 0.028834 0.038789
#>      3 0.776611 0.678786 0.769745 0.757530

table(prt_mz$status)
#> 
#>    1    2    3 
#> 1941  327 8963
xtabs(~ status + country, prt_mz)
#>       country
#> status Denmark Finland Norway Sweden
#>      1     586     308    300    747
#>      2      53      48     47    179
#>      3    2841     737   1367   4018
prop.table(xtabs(~ status + country, prt_mz), margin = 2L)
#>       country
#> status  Denmark  Finland   Norway   Sweden
#>      1 0.168391 0.281793 0.175029 0.151092
#>      2 0.015230 0.043916 0.027421 0.036206
#>      3 0.816379 0.674291 0.797550 0.812702
```

Then we fit separate models for DZ and MZ twins.

``` r
# fits a model given a data set
fit_model <- function(data, res_name){
  res_name <- file.path("cache", res_name)
  
  if(!file.exists(res_name)){
    get_comp_obj <- function(dat = data){
      ghq_data <- fastGHQuad::gaussHermiteData(10) |> 
        with(list(node = x, weight = w))
      
      mmcif_data(
        formula = ~ country - 1, data, cause = status, time = time, 
        cluster_id = id, max_time = max_time, spline_df = 5L, 
        strata = country, ghq_data = ghq_data)
    }
    comp_obj <- get_comp_obj()
    
    start_time <- 
      system.time(start_vals <- mmcif_start_values(comp_obj, n_threads = 4L))
    
    fit_time <- system.time(fit <- mmcif_fit(
      start_vals$upper, object = comp_obj, n_threads = 4L))
    
    ll_start <- mmcif_logLik(
      comp_obj, start_vals$upper, n_threads = 4L, is_log_chol = TRUE)
    
    gr_mle <- mmcif_logLik_grad(
      comp_obj, fit$par, n_threads = 4L, is_log_chol = TRUE)
    
    sandwich_time <- system.time(sandwich_est <- mmcif_sandwich(
      comp_obj, fit$par, n_threads = 4L, order = 0L))
    
    vcov_est <- log_chol_inv(fit$par[comp_obj$indices$vcov_upper])
    vcov_est[lower.tri(vcov_est)] <- NA_real_
    vcov_SE <- matrix(NA_real_, NROW(vcov_est), NCOL(vcov_est))
    vcov_SE[upper.tri(vcov_SE, TRUE)] <- 
      (attr(sandwich_est, "res vcov") |> diag() |> sqrt())[
        comp_obj$indices$vcov_upper]
    
    list(
      fit = fit, start = start_vals, ll_start = ll_start, 
      sandwich_est = sandwich_est, get_comp_obj = get_comp_obj, 
      time = list(start = start_time, fit = fit_time, 
                  sandwich = sandwich_time), 
      vcov_est = vcov_est, vcov_SE = vcov_SE, gr_mle = gr_mle) |> 
      saveRDS(res_name)
  }
  
  readRDS(res_name)
}

# fit models on the two types of twins
fit_dz <- fit_model(prt_dz, "fit_dz.RDS")
fit_mz <- fit_model(prt_mz, "fit_mz.RDS")
```

The results are summarized below.

``` r
# check convergence codes and counts
fit_dz$fit$convergence
#> [1] 0
fit_mz$fit$convergence
#> [1] 0
fit_dz$fit$counts
#> function gradient 
#>      659      657
fit_mz$fit$counts
#> function gradient 
#>      815      790

sqrt(sum(fit_dz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 0.087731
sqrt(sum(fit_mz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 0.26188

# the maximum log composite likelihood
print(-fit_dz$fit$value, digits = 10)
#> [1] -23833.42907
print(-fit_mz$fit$value, digits = 10)
#> [1] -12299.93411

# the time to estimate the model and to compute the sandwich estimator
with(fit_dz$time, start + fit + sandwich)
#>     user   system  elapsed 
#> 2309.747    0.091  579.794
with(fit_mz$time, start + fit + sandwich)
#>     user   system  elapsed 
#> 2041.195    0.248  511.015

# the parameter estimates along with standard errors
rbind(`DZ estimate` = fit_dz$fit$par,
      `DZ SE` = diag(fit_dz$sandwich_est) |> sqrt(),
      `MZ estimate` = fit_mz$fit$par,
      `MZ SE` = diag(fit_mz$sandwich_est) |> sqrt())
#>             cause1:risk:countryDenmark cause1:risk:countryFinland
#> DZ estimate                   0.724843                   0.637673
#> DZ SE                         0.068237                   0.091622
#> MZ estimate                   0.660979                   0.644025
#> MZ SE                         0.116715                   0.170183
#>             cause1:risk:countryNorway cause1:risk:countrySweden
#> DZ estimate                   0.42625                  0.433724
#> DZ SE                         0.10121                  0.057797
#> MZ estimate                   0.35726                  0.236986
#> MZ SE                         0.15380                  0.096953
#>             cause2:risk:countryDenmark cause2:risk:countryFinland
#> DZ estimate                   -2.47755                   -1.49759
#> DZ SE                          0.20252                    0.18404
#> MZ estimate                   -3.00505                   -2.29820
#> MZ SE                          0.36858                    0.39040
#>             cause2:risk:countryNorway cause2:risk:countrySweden
#> DZ estimate                  -1.89833                  -1.58886
#> DZ SE                         0.21030                   0.14789
#> MZ estimate                  -2.60269                  -2.17919
#> MZ SE                         0.37918                   0.27582
#>             cause1:strataDenmark:spline1 cause1:strataDenmark:spline2
#> DZ estimate                    -1.753065                    -2.343591
#> DZ SE                           0.067376                     0.076701
#> MZ estimate                    -1.757339                    -2.276793
#> MZ SE                           0.100200                     0.109660
#>             cause1:strataDenmark:spline3 cause1:strataDenmark:spline4
#> DZ estimate                    -3.109425                     -4.79099
#> DZ SE                           0.095186                      0.13491
#> MZ estimate                    -3.167256                     -4.61107
#> MZ SE                           0.135462                      0.18299
#>             cause1:strataDenmark:spline5 cause1:strataFinland:spline1
#> DZ estimate                     -3.96744                     -1.93999
#> DZ SE                            0.11904                      0.09887
#> MZ estimate                     -3.71391                     -2.08626
#> MZ SE                            0.14483                      0.15740
#>             cause1:strataFinland:spline2 cause1:strataFinland:spline3
#> DZ estimate                     -2.54388                     -3.11192
#> DZ SE                            0.11054                      0.12042
#> MZ estimate                     -2.60385                     -3.08892
#> MZ SE                            0.17197                      0.19312
#>             cause1:strataFinland:spline4 cause1:strataFinland:spline5
#> DZ estimate                     -5.04574                     -3.98129
#> DZ SE                            0.20826                      0.16389
#> MZ estimate                     -5.23771                     -4.10505
#> MZ SE                            0.31732                      0.19645
#>             cause1:strataNorway:spline1 cause1:strataNorway:spline2
#> DZ estimate                   -1.652542                    -2.26003
#> DZ SE                          0.096453                     0.10524
#> MZ estimate                   -1.977526                    -2.54615
#> MZ SE                          0.149576                     0.16522
#>             cause1:strataNorway:spline3 cause1:strataNorway:spline4
#> DZ estimate                    -2.92812                    -4.67719
#> DZ SE                           0.13291                     0.17465
#> MZ estimate                    -3.18134                    -5.37444
#> MZ SE                           0.18100                     0.30403
#>             cause1:strataNorway:spline5 cause1:strataSweden:spline1
#> DZ estimate                    -3.77936                   -1.923810
#> DZ SE                           0.14614                    0.086231
#> MZ estimate                    -4.39523                   -2.006642
#> MZ SE                           0.21357                    0.123061
#>             cause1:strataSweden:spline2 cause1:strataSweden:spline3
#> DZ estimate                   -2.579416                   -3.204515
#> DZ SE                          0.094813                    0.097127
#> MZ estimate                   -2.694445                   -3.354998
#> MZ SE                          0.133619                    0.133474
#>             cause1:strataSweden:spline4 cause1:strataSweden:spline5
#> DZ estimate                    -5.05609                    -4.02992
#> DZ SE                           0.17188                     0.10731
#> MZ estimate                    -5.16856                    -4.14860
#> MZ SE                           0.24646                     0.15127
#>             cause1:traject:countryDenmark cause1:traject:countryFinland
#> DZ estimate                      2.468410                       2.65896
#> DZ SE                            0.066602                       0.10017
#> MZ estimate                      2.650570                       2.78071
#> MZ SE                            0.107355                       0.16102
#>             cause1:traject:countryNorway cause1:traject:countrySweden
#> DZ estimate                     2.557636                     2.918105
#> DZ SE                           0.089394                     0.089729
#> MZ estimate                     3.040609                     3.122142
#> MZ SE                           0.158449                     0.138853
#>             cause2:strataDenmark:spline1 cause2:strataDenmark:spline2
#> DZ estimate                     -2.00570                     -2.83584
#> DZ SE                            0.26829                      0.31849
#> MZ estimate                     -3.10799                     -3.59102
#> MZ SE                            0.51976                      0.51078
#>             cause2:strataDenmark:spline3 cause2:strataDenmark:spline4
#> DZ estimate                     -3.04877                     -5.50706
#> DZ SE                            0.36258                      0.54860
#> MZ estimate                     -4.33119                     -6.79958
#> MZ SE                            0.64098                      0.82446
#>             cause2:strataDenmark:spline5 cause2:strataFinland:spline1
#> DZ estimate                     -4.17395                     -1.79497
#> DZ SE                            0.42959                      0.23083
#> MZ estimate                     -5.92990                     -2.72471
#> MZ SE                            0.82164                      0.47389
#>             cause2:strataFinland:spline2 cause2:strataFinland:spline3
#> DZ estimate                     -2.03362                     -2.93599
#> DZ SE                            0.24308                      0.34395
#> MZ estimate                     -3.18184                     -4.46812
#> MZ SE                            0.47399                      0.69657
#>             cause2:strataFinland:spline4 cause2:strataFinland:spline5
#> DZ estimate                     -4.81308                     -3.99673
#> DZ SE                            0.46961                      0.40678
#> MZ estimate                     -7.00528                     -4.67909
#> MZ SE                            0.95679                      0.65856
#>             cause2:strataNorway:spline1 cause2:strataNorway:spline2
#> DZ estimate                    -2.24983                    -2.53974
#> DZ SE                           0.29842                     0.30499
#> MZ estimate                    -3.21152                    -3.38075
#> MZ SE                           0.60099                     0.61389
#>             cause2:strataNorway:spline3 cause2:strataNorway:spline4
#> DZ estimate                    -3.38967                    -4.94046
#> DZ SE                           0.39273                     0.50166
#> MZ estimate                    -4.72648                    -7.16877
#> MZ SE                           0.66761                     1.18818
#>             cause2:strataNorway:spline5 cause2:strataSweden:spline1
#> DZ estimate                    -3.85490                    -1.93650
#> DZ SE                           0.42199                     0.19048
#> MZ estimate                    -5.77992                    -2.54150
#> MZ SE                           0.85752                     0.29768
#>             cause2:strataSweden:spline2 cause2:strataSweden:spline3
#> DZ estimate                    -2.52216                    -3.04438
#> DZ SE                           0.21940                     0.26701
#> MZ estimate                    -3.17232                    -3.81595
#> MZ SE                           0.34038                     0.40524
#>             cause2:strataSweden:spline4 cause2:strataSweden:spline5
#> DZ estimate                    -5.06469                    -3.64861
#> DZ SE                           0.42761                     0.31946
#> MZ estimate                    -6.40343                    -4.86401
#> MZ SE                           0.63474                     0.52031
#>             cause2:traject:countryDenmark cause2:traject:countryFinland
#> DZ estimate                       2.86485                       2.69866
#> DZ SE                             0.32375                       0.29799
#> MZ estimate                       4.85647                       4.62852
#> MZ SE                             0.60668                       0.55911
#>             cause2:traject:countryNorway cause2:traject:countrySweden
#> DZ estimate                      2.82290                      2.86566
#> DZ SE                            0.32192                      0.27108
#> MZ estimate                      5.01252                      4.35518
#> MZ SE                            0.66260                      0.44357
#>             vcov:risk1:risk1 vcov:risk1:risk2 vcov:risk2:risk2
#> DZ estimate         -0.35589          0.23443          0.12023
#> DZ SE                0.14527          0.28126          0.14895
#> MZ estimate          0.20185          0.67922          0.70662
#> MZ SE                0.11637          0.53492          0.12816
#>             vcov:risk1:traject1 vcov:risk2:traject1 vcov:traject1:traject1
#> DZ estimate           -0.249749           -0.049394                -0.9488
#> DZ SE                  0.086135            0.113239                 0.2252
#> MZ estimate           -0.527129           -0.067804                -4.3329
#> MZ SE                  0.071186            0.256263                 1.4000
#>             vcov:risk1:traject2 vcov:risk2:traject2 vcov:traject1:traject2
#> DZ estimate            -0.28278            -0.35555               -0.24542
#> DZ SE                   0.18897             0.16151                0.22180
#> MZ estimate            -0.14784            -0.91285               -0.80571
#> MZ SE                   0.33947             0.24349                0.28907
#>             vcov:traject2:traject2
#> DZ estimate               -2.52041
#> DZ SE                      0.68679
#> MZ estimate               -2.35689
#> MZ SE                      2.27580

# the estimated covariance matrices along with standard errors
fit_dz$vcov_est # the estimates for DZ
#>         [,1]    [,2]     [,3]       [,4]
#> [1,] 0.49077 0.16423 -0.17496 -0.1980980
#> [2,]      NA 1.32679 -0.11425 -0.4672686
#> [3,]      NA      NA  0.21474 -0.0068421
#> [4,]      NA      NA       NA  0.2730783
fit_dz$vcov_SE # the SEs
#>         [,1]    [,2]     [,3]     [,4]
#> [1,] 0.14258 0.20441 0.055565 0.137138
#> [2,]      NA 0.42385 0.105609 0.182604
#> [3,]      NA      NA 0.047271 0.084765
#> [4,]      NA      NA       NA 0.194869
fit_mz$vcov_est # the estimates for MZ
#>        [,1]    [,2]     [,3]     [,4]
#> [1,] 1.4974 0.83114 -0.64503 -0.18091
#> [2,]     NA 4.57061 -0.49549 -1.95089
#> [3,]     NA      NA  0.28264  0.12925
#> [4,]     NA      NA       NA  1.51330
fit_mz$vcov_SE # the SEs
#>        [,1]   [,2]     [,3]    [,4]
#> [1,] 0.3485 0.6695 0.088133 0.41438
#> [2,]     NA 1.3781 0.311675 0.60304
#> [3,]     NA     NA 0.073664 0.16610
#> [4,]     NA     NA       NA 0.48133
```

``` r
# plot the conditional cumulative incidence functions
plot_curves <- function(fit){
  comp_obj <- fit$get_comp_obj()
  
  # get the estimates
  n_causes <- 2L
  coef_risk_est <- fit$fit$par[comp_obj$indices$coef_risk] |> 
    matrix(ncol = n_causes)
  coef_traject_time_est <- fit$fit$par[comp_obj$indices$coef_trajectory_time] |> 
    matrix(ncol = n_causes)
  coef_traject_est <- fit$fit$par[comp_obj$indices$coef_trajectory] |> 
    matrix(ncol = n_causes)
  
  # plot the estimated cumulative incidence functions
  par(mar = c(5, 5, 1, 1), mfcol = c(1, 2))
  for(i in 1:2){
    pts <- seq(min(prt$time), max_time * (1 - 1e-8), length.out = 1000)
    curves <- sapply(1:4 - 1L, \(strata){
      # compute the risk probabilities  
      probs_est <- exp(coef_risk_est[1 + strata, ]) / 
        (1 + sum(exp(coef_risk_est[1 + strata, ])))
      
      probs_est[i] * pnorm(
        -comp_obj$time_expansion(pts, cause = i, which_strata = strata + 1L) %*% 
          coef_traject_time_est[, i] - 
          coef_traject_est[21 + strata, i]) |> drop()
    })
    
    matplot(
      pts, curves, ylim = c(0, 1), bty = "l",  xlab = "time", type = "l",
      ylab = sprintf("Cumulative incidence; cause %d", i), lty = 1:4, 
      col = "black", yaxs = "i", xaxs = "i")
    grid(ny = 10)
    legend("topleft", legend = levels(prt$country), bty = "n", 
           cex = par("cex") * .8, lty = 1:4)
  }
}

plot_curves(fit_dz)
```

<img src="figures/prostate/plots-1.png" width="100%" />

``` r
plot_curves(fit_mz)
```

<img src="figures/prostate/plots-2.png" width="100%" />

``` r
# compare to marginal estimates (NOT directly comparable; the above are 
# conditional)
par(mar = c(5, 5, 1, 1))
plot_cif <- \(data, cause){
  cif(Event(time, status) ~ strata(country), data = data,
              cause = cause, cens.code = 3) |>
    bplot(se = TRUE, ylim = c(0, 1), bty = "l", yaxs = "i", xaxs = "i")
  grid(ny = 10)
}

plot_cif(prt_dz, 1)
plot_cif(prt_dz, 2)
```

<img src="figures/prostate/plots-3.png" width="100%" />

``` r
plot_cif(prt_mz, 1)
plot_cif(prt_mz, 2)
```

<img src="figures/prostate/plots-4.png" width="100%" />
