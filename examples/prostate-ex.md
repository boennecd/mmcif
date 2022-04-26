
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
    ghq_list <- lapply(
      c(4, 10), \(n_nodes)
        fastGHQuad::gaussHermiteData(n_nodes) |> 
          with(list(node = x, weight = w)))
    
    get_comp_obj <- function(dat = data){
      mmcif_data(
        formula = ~ country - 1, data, cause = status, time = time, 
        cluster_id = id, max_time = max_time, spline_df = 5L, 
        strata = country, ghq_data = ghq_list[[length(ghq_list)]])
    }
    comp_obj <- get_comp_obj()
    
    start_time <- 
      system.time(start_vals <- mmcif_start_values(comp_obj, n_threads = 4L))
    
    fit_time <- system.time(fits <- mmcif_fit(
      start_vals$upper, object = comp_obj, n_threads = 4L, ghq_data = ghq_list))
    
    fit <- fits[[length(fits)]]
    
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
      fit = fit, fits = fits, start = start_vals, ll_start = ll_start, 
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
#> [1] 1
fit_mz$fit$convergence
#> [1] 1
fit_dz$fit$counts
#> function gradient 
#>      286      234
fit_mz$fit$counts
#> function gradient 
#>      194      114

sqrt(sum(fit_dz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 0.0034399
sqrt(sum(fit_mz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 0.34094

# the maximum log composite likelihood
print(-fit_dz$fit$value, digits = 10)
#> [1] -23833.41659
print(-fit_mz$fit$value, digits = 10)
#> [1] -12296.7432

# the time to estimate the model and to compute the sandwich estimator
with(fit_dz$time, start + fit + sandwich)
#>     user   system  elapsed 
#> 1257.739    0.293  315.381
with(fit_mz$time, start + fit + sandwich)
#>    user  system elapsed 
#>  779.87    0.04  195.86

# the parameter estimates along with standard errors
rbind(`DZ estimate` = fit_dz$fit$par,
      `DZ SE` = diag(fit_dz$sandwich_est) |> sqrt(),
      `MZ estimate` = fit_mz$fit$par,
      `MZ SE` = diag(fit_mz$sandwich_est) |> sqrt())
#>             cause1:risk:countryDenmark cause1:risk:countryFinland
#> DZ estimate                   0.724824                   0.637623
#> DZ SE                         0.068234                   0.091616
#> MZ estimate                   0.685456                   0.660546
#> MZ SE                         0.119975                   0.174573
#>             cause1:risk:countryNorway cause1:risk:countrySweden
#> DZ estimate                   0.42608                  0.433614
#> DZ SE                         0.10120                  0.057784
#> MZ estimate                   0.37892                  0.249915
#> MZ SE                         0.15766                  0.092680
#>             cause2:risk:countryDenmark cause2:risk:countryFinland
#> DZ estimate                   -2.47811                   -1.49776
#> DZ SE                          0.20242                    0.18394
#> MZ estimate                   -2.95504                   -2.25473
#> MZ SE                          0.36929                    0.39135
#>             cause2:risk:countryNorway cause2:risk:countrySweden
#> DZ estimate                  -1.89892                  -1.58944
#> DZ SE                         0.21023                   0.14783
#> MZ estimate                  -2.55645                  -2.14015
#> MZ SE                         0.38218                   0.27561
#>             cause1:strataDenmark:spline1 cause1:strataDenmark:spline2
#> DZ estimate                    -1.752890                    -2.343360
#> DZ SE                           0.067369                     0.076693
#> MZ estimate                    -1.829579                    -2.372409
#> MZ SE                           0.107167                     0.118632
#>             cause1:strataDenmark:spline3 cause1:strataDenmark:spline4
#> DZ estimate                     -3.10914                     -4.79053
#> DZ SE                            0.09518                      0.13490
#> MZ estimate                     -3.30410                     -4.80048
#> MZ SE                            0.14819                      0.20233
#>             cause1:strataDenmark:spline5 cause1:strataFinland:spline1
#> DZ estimate                     -3.96712                    -1.939739
#> DZ SE                            0.11904                     0.098856
#> MZ estimate                     -3.87288                    -2.170988
#> MZ SE                            0.16104                     0.165674
#>             cause1:strataFinland:spline2 cause1:strataFinland:spline3
#> DZ estimate                     -2.54363                     -3.11157
#> DZ SE                            0.11053                      0.12042
#> MZ estimate                     -2.71170                     -3.21849
#> MZ SE                            0.18180                      0.20469
#>             cause1:strataFinland:spline4 cause1:strataFinland:spline5
#> DZ estimate                     -5.04531                     -3.98108
#> DZ SE                            0.20825                      0.16393
#> MZ estimate                     -5.45758                     -4.28318
#> MZ SE                            0.33643                      0.20879
#>             cause1:strataNorway:spline1 cause1:strataNorway:spline2
#> DZ estimate                   -1.652212                    -2.25964
#> DZ SE                          0.096419                     0.10521
#> MZ estimate                   -2.057418                    -2.64852
#> MZ SE                          0.156670                     0.17345
#>             cause1:strataNorway:spline3 cause1:strataNorway:spline4
#> DZ estimate                    -2.92781                    -4.67646
#> DZ SE                           0.13290                     0.17459
#> MZ estimate                    -3.32339                    -5.59961
#> MZ SE                           0.19651                     0.32480
#>             cause1:strataNorway:spline5 cause1:strataSweden:spline1
#> DZ estimate                    -3.77917                   -1.923488
#> DZ SE                           0.14616                    0.086202
#> MZ estimate                    -4.58154                   -2.089501
#> MZ SE                           0.23443                    0.131644
#>             cause1:strataSweden:spline2 cause1:strataSweden:spline3
#> DZ estimate                   -2.579064                   -3.204194
#> DZ SE                          0.094784                    0.097116
#> MZ estimate                   -2.810152                   -3.501061
#> MZ SE                          0.145351                    0.148013
#>             cause1:strataSweden:spline4 cause1:strataSweden:spline5
#> DZ estimate                    -5.05543                    -4.02964
#> DZ SE                           0.17183                     0.10731
#> MZ estimate                    -5.38148                    -4.32715
#> MZ SE                           0.26738                     0.16819
#>             cause1:traject:countryDenmark cause1:traject:countryFinland
#> DZ estimate                      2.468221                       2.65870
#> DZ SE                            0.066592                       0.10016
#> MZ estimate                      2.689367                       2.82308
#> MZ SE                            0.105992                       0.16374
#>             cause1:traject:countryNorway cause1:traject:countrySweden
#> DZ estimate                      2.55725                     2.917746
#> DZ SE                            0.08935                     0.089696
#> MZ estimate                      3.08996                     3.171599
#> MZ SE                            0.15961                     0.138517
#>             cause2:strataDenmark:spline1 cause2:strataDenmark:spline2
#> DZ estimate                     -1.99878                     -2.82643
#> DZ SE                            0.26596                      0.31505
#> MZ estimate                     -3.16236                     -3.65558
#> MZ SE                            0.53495                      0.52683
#>             cause2:strataDenmark:spline3 cause2:strataDenmark:spline4
#> DZ estimate                     -3.03970                     -5.48905
#> DZ SE                            0.35963                      0.54185
#> MZ estimate                     -4.39835                     -6.92122
#> MZ SE                            0.65915                      0.85487
#>             cause2:strataDenmark:spline5 cause2:strataFinland:spline1
#> DZ estimate                     -4.16122                     -1.78931
#> DZ SE                            0.42491                      0.22873
#> MZ estimate                     -6.02251                     -2.77921
#> MZ SE                            0.84458                      0.48433
#>             cause2:strataFinland:spline2 cause2:strataFinland:spline3
#> DZ estimate                     -2.02752                     -2.92586
#> DZ SE                            0.24072                      0.34021
#> MZ estimate                     -3.24612                     -4.55765
#> MZ SE                            0.48323                      0.71519
#>             cause2:strataFinland:spline4 cause2:strataFinland:spline5
#> DZ estimate                     -4.79712                     -3.98414
#> DZ SE                            0.46346                      0.40213
#> MZ estimate                     -7.15681                     -4.76385
#> MZ SE                            0.98427                      0.67660
#>             cause2:strataNorway:spline1 cause2:strataNorway:spline2
#> DZ estimate                    -2.24325                    -2.53184
#> DZ SE                           0.29634                     0.30242
#> MZ estimate                    -3.28357                    -3.45221
#> MZ SE                           0.63140                     0.64435
#>             cause2:strataNorway:spline3 cause2:strataNorway:spline4
#> DZ estimate                    -3.37990                    -4.92537
#> DZ SE                           0.38960                     0.49685
#> MZ estimate                    -4.79748                    -7.33208
#> MZ SE                           0.69253                     1.25061
#>             cause2:strataNorway:spline5 cause2:strataSweden:spline1
#> DZ estimate                    -3.84369                    -1.92948
#> DZ SE                           0.41809                     0.18808
#> MZ estimate                    -5.88214                    -2.57661
#> MZ SE                           0.88241                     0.30388
#>             cause2:strataSweden:spline2 cause2:strataSweden:spline3
#> DZ estimate                    -2.51348                    -3.03399
#> DZ SE                           0.21645                     0.26321
#> MZ estimate                    -3.21703                    -3.87027
#> MZ SE                           0.34822                     0.41438
#>             cause2:strataSweden:spline4 cause2:strataSweden:spline5
#> DZ estimate                    -5.04659                    -3.63682
#> DZ SE                           0.42090                     0.31498
#> MZ estimate                    -6.49379                    -4.93436
#> MZ SE                           0.64843                     0.53163
#>             cause2:traject:countryDenmark cause2:traject:countryFinland
#> DZ estimate                       2.85592                       2.69060
#> DZ SE                             0.32090                       0.29513
#> MZ estimate                       4.93565                       4.70827
#> MZ SE                             0.62366                       0.57432
#>             cause2:traject:countryNorway cause2:traject:countrySweden
#> DZ estimate                      2.81451                      2.85627
#> DZ SE                            0.31941                      0.26803
#> MZ estimate                      5.11589                      4.41663
#> MZ SE                            0.69723                      0.45377
#>             vcov:risk1:risk1 vcov:risk1:risk2 vcov:risk2:risk2
#> DZ estimate        -0.356432          0.23367          0.12078
#> DZ SE               0.145505          0.28122          0.14863
#> MZ estimate         0.290451          0.81292          0.70174
#> MZ SE               0.099517          0.34984          0.12748
#>             vcov:risk1:traject1 vcov:risk2:traject1 vcov:traject1:traject1
#> DZ estimate           -0.249975           -0.049630               -0.94979
#> DZ SE                  0.086241            0.113268                0.22593
#> MZ estimate           -0.333847           -0.033447               -0.69992
#> MZ SE                  0.096794            0.118830                0.19452
#>             vcov:risk1:traject2 vcov:risk2:traject2 vcov:traject1:traject2
#> DZ estimate            -0.28212            -0.35485               -0.24491
#> DZ SE                   0.18844             0.16085                0.22137
#> MZ estimate            -0.23030            -0.92243               -0.08474
#> MZ SE                   0.32900             0.24830                0.50169
#>             vcov:traject2:traject2
#> DZ estimate               -7.35447
#> DZ SE                      0.69221
#> MZ estimate               -0.19197
#> MZ SE                      0.37014

# the estimated covariance matrices along with standard errors
fit_dz$vcov_est # the estimates for DZ
#>         [,1]    [,2]     [,3]       [,4]
#> [1,] 0.49024 0.16361 -0.17503 -0.1975319
#> [2,]      NA 1.32784 -0.11441 -0.4663260
#> [3,]      NA      NA  0.21458 -0.0066027
#> [4,]      NA      NA       NA  0.2654897
fit_dz$vcov_SE # the SEs
#>         [,1]    [,2]     [,3]     [,4]
#> [1,] 0.14266 0.20428 0.055575 0.136673
#> [2,]      NA 0.42322 0.105646 0.181994
#> [3,]      NA      NA 0.047268 0.084483
#> [4,]      NA      NA       NA 0.190917
fit_mz$vcov_est # the estimates for MZ
#>        [,1]   [,2]     [,3]      [,4]
#> [1,] 1.7876 1.0869 -0.44636 -0.307924
#> [2,]     NA 4.7302 -0.33886 -2.048004
#> [3,]     NA     NA  0.35921  0.065655
#> [4,]     NA     NA       NA  1.592272
fit_mz$vcov_SE # the SEs
#>        [,1]    [,2]     [,3]    [,4]
#> [1,] 0.3558 0.50537 0.122935 0.43945
#> [2,]     NA 1.22903 0.236608 0.59248
#> [3,]     NA      NA 0.076067 0.19166
#> [4,]     NA      NA       NA 0.50428
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
