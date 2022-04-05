
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
max_time <- 100
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
#>  4524   574 12893
xtabs(~ status + country, prt_dz)
#>       country
#> status Denmark Finland Norway Sweden
#>      1    1487     860    534   1643
#>      2      89     129     75    281
#>      3    4615    1844   1784   4650
prop.table(xtabs(~ status + country, prt_dz), margin = 2L)
#>       country
#> status  Denmark  Finland   Norway   Sweden
#>      1 0.240187 0.303565 0.223151 0.249924
#>      2 0.014376 0.045535 0.031341 0.042744
#>      3 0.745437 0.650900 0.745508 0.707332

table(prt_mz$status)
#> 
#>    1    2    3 
#> 2295  357 8579
xtabs(~ status + country, prt_mz)
#>       country
#> status Denmark Finland Norway Sweden
#>      1     690     335    331    939
#>      2      56      54     52    195
#>      3    2734     704   1331   3810
prop.table(xtabs(~ status + country, prt_mz), margin = 2L)
#>       country
#> status  Denmark  Finland   Norway   Sweden
#>      1 0.198276 0.306496 0.193116 0.189927
#>      2 0.016092 0.049405 0.030338 0.039442
#>      3 0.785632 0.644099 0.776546 0.770631
```

Then we fit separate models for DZ and MZ twins.

``` r
# fits a model given a data set
fit_model <- function(data, res_name){
  res_name <- file.path("cache", res_name)
  
  if(!file.exists(res_name)){
    get_comp_obj <- function(dat = data){
      ghq_data <- fastGHQuad::gaussHermiteData(8) |> 
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
#>      201      199
fit_mz$fit$counts
#> function gradient 
#>      444      442

sqrt(sum(fit_dz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 1.2629
sqrt(sum(fit_mz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 0.63822

# the maximum log composite likelihood
print(-fit_dz$fit$value, digits = 10)
#> [1] -26076.95688
print(-fit_mz$fit$value, digits = 10)
#> [1] -13596.91875

# the time to estimate the model and to compute the sandwich estimator
with(fit_dz$time, start + fit + sandwich)
#>    user  system elapsed 
#> 3248.51    0.12  817.58
with(fit_mz$time, start + fit + sandwich)
#>     user   system  elapsed 
#> 4559.047    0.132 1147.912

# the parameter estimates along with standard errors
rbind(`DZ estimate` = fit_dz$fit$par,
      `DZ SE` = diag(fit_dz$sandwich_est) |> sqrt(),
      `MZ estimate` = fit_mz$fit$par,
      `MZ SE` = diag(fit_mz$sandwich_est) |> sqrt())
#>             cause1:risk:countryDenmark cause1:risk:countryFinland
#> DZ estimate                    2.61361                    2.42749
#> DZ SE                          0.17676                    0.25049
#> MZ estimate                    3.19929                    3.02735
#> MZ SE                          0.38008                    0.53446
#>             cause1:risk:countryNorway cause1:risk:countrySweden
#> DZ estimate                   1.96351                   2.40944
#> DZ SE                         0.23792                   0.15337
#> MZ estimate                   2.35121                   2.54771
#> MZ SE                         0.44806                   0.34863
#>             cause2:risk:countryDenmark cause2:risk:countryFinland
#> DZ estimate                   -0.63073                   0.283836
#> DZ SE                          0.24895                   0.285217
#> MZ estimate                   -1.23237                  -0.091496
#> MZ SE                          0.56369                   0.689388
#>             cause2:risk:countryNorway cause2:risk:countrySweden
#> DZ estimate                  -0.40363                   0.21402
#> DZ SE                         0.28902                   0.19986
#> MZ estimate                  -1.01869                  -0.54239
#> MZ SE                         0.59682                   0.43709
#>             cause1:strata1:spline1 cause1:strata1:spline2
#> DZ estimate              -1.663919              -2.245965
#> DZ SE                     0.057753               0.062916
#> MZ estimate              -1.804480              -2.471069
#> MZ SE                     0.096531               0.106469
#>             cause1:strata1:spline3 cause1:strata1:spline4
#> DZ estimate              -2.922424               -4.42857
#> DZ SE                     0.077088                0.10938
#> MZ estimate              -3.203297               -4.85757
#> MZ SE                     0.128754                0.17844
#>             cause1:strata1:spline5 cause1:strata2:spline1
#> DZ estimate              -3.615946               -1.85906
#> DZ SE                     0.096514                0.08622
#> MZ estimate              -4.081678               -2.16424
#> MZ SE                     0.136181                0.15350
#>             cause1:strata2:spline2 cause1:strata2:spline3
#> DZ estimate               -2.35595               -3.06925
#> DZ SE                      0.09478                0.10606
#> MZ estimate               -2.66464               -3.45803
#> MZ SE                      0.16972                0.18942
#>             cause1:strata2:spline4 cause1:strata2:spline5
#> DZ estimate               -4.74154               -3.71601
#> DZ SE                      0.17936                0.13435
#> MZ estimate               -5.43727               -4.13591
#> MZ SE                      0.32127                0.26233
#>             cause1:strata3:spline1 cause1:strata3:spline2
#> DZ estimate              -1.578850               -2.14677
#> DZ SE                     0.085618                0.09156
#> MZ estimate              -2.055968               -2.73611
#> MZ SE                     0.146073                0.16381
#>             cause1:strata3:spline3 cause1:strata3:spline4
#> DZ estimate               -2.99273               -4.47494
#> DZ SE                      0.11897                0.16600
#> MZ estimate               -3.57439               -5.57758
#> MZ SE                      0.19897                0.31662
#>             cause1:strata3:spline5 cause1:strata4:spline1
#> DZ estimate               -3.57662              -1.857043
#> DZ SE                      0.18003               0.073118
#> MZ estimate               -4.37820              -2.131585
#> MZ SE                      0.23461               0.116497
#>             cause1:strata4:spline2 cause1:strata4:spline3
#> DZ estimate              -2.406163               -3.16528
#> DZ SE                     0.076889                0.07607
#> MZ estimate              -2.778665               -3.56065
#> MZ SE                     0.124596                0.12467
#>             cause1:strata4:spline4 cause1:strata4:spline5
#> DZ estimate               -4.82022              -3.877115
#> DZ SE                      0.13707               0.083305
#> MZ estimate               -5.39588              -4.431196
#> MZ SE                      0.22809               0.139562
#>             cause1:traject:countryDenmark cause1:traject:countryFinland
#> DZ estimate                      2.470098                       2.63480
#> DZ SE                            0.055718                       0.08563
#> MZ estimate                      2.807306                       2.91035
#> MZ SE                            0.093041                       0.14969
#>             cause1:traject:countryNorway cause1:traject:countrySweden
#> DZ estimate                      2.57009                     2.910965
#> DZ SE                            0.07755                     0.073601
#> MZ estimate                      3.19889                     3.293960
#> MZ SE                            0.14680                     0.119748
#>             cause2:strata1:spline1 cause2:strata1:spline2
#> DZ estimate               -1.90069               -2.52779
#> DZ SE                      0.24712                0.28366
#> MZ estimate               -3.11175               -3.25448
#> MZ SE                      0.46993                0.43189
#>             cause2:strata1:spline3 cause2:strata1:spline4
#> DZ estimate               -2.91933               -4.74403
#> DZ SE                      0.32628                0.48780
#> MZ estimate               -4.61283               -6.41469
#> MZ SE                      0.65648                0.79033
#>             cause2:strata1:spline5 cause2:strata2:spline1
#> DZ estimate               -3.24548               -1.66577
#> DZ SE                      0.35864                0.21686
#> MZ estimate               -5.37722               -2.30169
#> MZ SE                      0.65331                0.42617
#>             cause2:strata2:spline2 cause2:strata2:spline3
#> DZ estimate               -1.95941               -2.91686
#> DZ SE                      0.24507                0.31346
#> MZ estimate               -2.99028               -3.47593
#> MZ SE                      0.46501                0.55819
#>             cause2:strata2:spline4 cause2:strata2:spline5
#> DZ estimate               -4.64161               -3.72661
#> DZ SE                      0.47590                0.43232
#> MZ estimate               -5.93081               -3.96342
#> MZ SE                      0.86665                0.58298
#>             cause2:strata3:spline1 cause2:strata3:spline2
#> DZ estimate               -2.10015               -2.33791
#> DZ SE                      0.27341                0.27093
#> MZ estimate               -2.95651               -3.03318
#> MZ SE                      0.53555                0.52712
#>             cause2:strata3:spline3 cause2:strata3:spline4
#> DZ estimate               -3.00421               -4.56059
#> DZ SE                      0.40052                0.51688
#> MZ estimate               -4.34365               -6.35231
#> MZ SE                      0.76397                1.09850
#>             cause2:strata3:spline5 cause2:strata4:spline1
#> DZ estimate               -3.56983               -1.87074
#> DZ SE                      0.45654                0.18021
#> MZ estimate               -4.95764               -2.37276
#> MZ SE                      0.72656                0.27260
#>             cause2:strata4:spline2 cause2:strata4:spline3
#> DZ estimate               -2.36296               -3.06051
#> DZ SE                      0.20806                0.27206
#> MZ estimate               -2.90484               -3.83225
#> MZ SE                      0.29595                0.43376
#>             cause2:strata4:spline4 cause2:strata4:spline5
#> DZ estimate               -4.90909               -3.52812
#> DZ SE                      0.41638                0.30888
#> MZ estimate               -6.19232               -4.89533
#> MZ SE                      0.63485                0.54244
#>             cause2:traject:countryDenmark cause2:traject:countryFinland
#> DZ estimate                       2.81590                       2.71443
#> DZ SE                             0.28742                       0.27129
#> MZ estimate                       3.84908                       3.79009
#> MZ SE                             0.58777                       0.54203
#>             cause2:traject:countryNorway cause2:traject:countrySweden
#> DZ estimate                      2.79397                      2.82131
#> DZ SE                            0.29095                      0.24475
#> MZ estimate                      4.13282                      3.48920
#> MZ SE                            0.60209                      0.42026
#>             vcov:risk1:risk1 vcov:risk1:risk2 vcov:risk2:risk2
#> DZ estimate       -0.0076402          0.46215         0.028833
#> DZ SE              0.1787029          0.31469         0.174708
#> MZ estimate        0.3479495          2.91706         0.597730
#> MZ SE              0.2608980          0.99217         0.345149
#>             vcov:risk1:traject1 vcov:risk2:traject1 vcov:traject1:traject1
#> DZ estimate           -0.393229           -0.068479               -3.32560
#> DZ SE                  0.046394            0.117828                7.50821
#> MZ estimate           -0.499229            0.457140               -4.30561
#> MZ SE                  0.186123            0.211756                0.96085
#>             vcov:risk1:traject2 vcov:risk2:traject2 vcov:traject1:traject2
#> DZ estimate            -0.23016           -0.406638              -0.055195
#> DZ SE                   0.17754            0.170609               1.145687
#> MZ estimate            -0.40880           -0.048982               0.919617
#> MZ SE                   0.31638            0.271958               0.210227
#>             vcov:traject2:traject2
#> DZ estimate               -2.22752
#> DZ SE                      1.74194
#> MZ estimate               -2.39998
#> MZ SE                      0.53358

# the estimated covariance matrices along with standard errors
fit_dz$vcov_est # the estimates for DZ
#>         [,1]    [,2]     [,3]     [,4]
#> [1,] 0.98484 0.45863 -0.39024 -0.22840
#> [2,]      NA 1.27294 -0.25221 -0.52490
#> [3,]      NA      NA  0.16061  0.11637
#> [4,]      NA      NA       NA  0.23299
fit_dz$vcov_SE # the SEs
#>         [,1]    [,2]     [,3]    [,4]
#> [1,] 0.35199 0.71348 0.146185 0.36592
#> [2,]      NA 0.53660 0.197034 0.37581
#> [3,]      NA      NA 0.030811 0.13221
#> [4,]      NA      NA       NA 0.18604
fit_mz$vcov_est # the estimates for MZ
#>        [,1]   [,2]     [,3]     [,4]
#> [1,] 2.0055  4.131 -0.70699 -0.57893
#> [2,]     NA 11.814 -0.62521 -1.28155
#> [3,]     NA     NA  0.45839  0.19410
#> [4,]     NA     NA       NA  1.02344
fit_mz$vcov_SE # the SEs
#>        [,1]   [,2]     [,3]    [,4]
#> [1,] 1.0465 4.6533 0.419097 0.94750
#> [2,]     NA 4.1705 0.755531 2.03162
#> [3,]     NA     NA 0.067427 0.35209
#> [4,]     NA     NA       NA 0.34643
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
