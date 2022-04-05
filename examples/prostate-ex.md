
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
#>      201      199
fit_mz$fit$counts
#> function gradient 
#>      409      400

sqrt(sum(fit_dz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 1.2599
sqrt(sum(fit_mz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 1.308

# the maximum log composite likelihood
print(-fit_dz$fit$value, digits = 10)
#> [1] -26076.95589
print(-fit_mz$fit$value, digits = 10)
#> [1] -13596.94963

# the time to estimate the model and to compute the sandwich estimator
with(fit_dz$time, start + fit + sandwich)
#>     user   system  elapsed 
#> 4720.476    0.196 1200.550
with(fit_mz$time, start + fit + sandwich)
#>     user   system  elapsed 
#> 6307.351    0.188 1611.355

# the parameter estimates along with standard errors
rbind(`DZ estimate` = fit_dz$fit$par,
      `DZ SE` = diag(fit_dz$sandwich_est) |> sqrt(),
      `MZ estimate` = fit_mz$fit$par,
      `MZ SE` = diag(fit_mz$sandwich_est) |> sqrt())
#>             cause1:risk:countryDenmark cause1:risk:countryFinland
#> DZ estimate                    2.61360                    2.42751
#> DZ SE                          0.17681                    0.25053
#> MZ estimate                    3.20997                    3.01315
#> MZ SE                          0.38910                    0.53142
#>             cause1:risk:countryNorway cause1:risk:countrySweden
#> DZ estimate                   1.96356                   2.40942
#> DZ SE                         0.23796                   0.15341
#> MZ estimate                   2.33990                   2.55787
#> MZ SE                         0.45071                   0.36083
#>             cause2:risk:countryDenmark cause2:risk:countryFinland
#> DZ estimate                   -0.63053                    0.28404
#> DZ SE                          0.24892                    0.28524
#> MZ estimate                   -1.23188                   -0.12168
#> MZ SE                          0.55980                    0.68112
#>             cause2:risk:countryNorway cause2:risk:countrySweden
#> DZ estimate                  -0.40335                   0.21419
#> DZ SE                         0.28900                   0.19992
#> MZ estimate                  -1.04540                  -0.53727
#> MZ SE                         0.59321                   0.43812
#>             cause1:strata1:spline1 cause1:strata1:spline2
#> DZ estimate              -1.663905              -2.245957
#> DZ SE                     0.057758               0.062926
#> MZ estimate              -1.803013              -2.471495
#> MZ SE                     0.096523               0.106586
#>             cause1:strata1:spline3 cause1:strata1:spline4
#> DZ estimate              -2.922416               -4.42855
#> DZ SE                     0.077103                0.10940
#> MZ estimate              -3.203099               -4.85792
#> MZ SE                     0.129035                0.17889
#>             cause1:strata1:spline5 cause1:strata2:spline1
#> DZ estimate              -3.615930              -1.859049
#> DZ SE                     0.096534               0.086223
#> MZ estimate              -4.083444              -2.161915
#> MZ SE                     0.136992               0.152840
#>             cause1:strata2:spline2 cause1:strata2:spline3
#> DZ estimate              -2.355930               -3.06923
#> DZ SE                     0.094785                0.10607
#> MZ estimate              -2.663535               -3.45699
#> MZ SE                     0.169078                0.18921
#>             cause1:strata2:spline4 cause1:strata2:spline5
#> DZ estimate               -4.74150               -3.71597
#> DZ SE                      0.17938                0.13437
#> MZ estimate               -5.42963               -4.13529
#> MZ SE                      0.31955                0.26187
#>             cause1:strata3:spline1 cause1:strata3:spline2
#> DZ estimate               -1.57884              -2.146746
#> DZ SE                      0.08562               0.091565
#> MZ estimate               -2.06277              -2.744967
#> MZ SE                      0.14760               0.165525
#>             cause1:strata3:spline3 cause1:strata3:spline4
#> DZ estimate               -2.99270               -4.47489
#> DZ SE                      0.11897                0.16601
#> MZ estimate               -3.58202               -5.59317
#> MZ SE                      0.19985                0.31951
#>             cause1:strata3:spline5 cause1:strata4:spline1
#> DZ estimate               -3.57659              -1.857027
#> DZ SE                      0.18004               0.073123
#> MZ estimate               -4.37903              -2.131680
#> MZ SE                      0.23321               0.116702
#>             cause1:strata4:spline2 cause1:strata4:spline3
#> DZ estimate              -2.406147              -3.165261
#> DZ SE                     0.076897               0.076088
#> MZ estimate              -2.779404              -3.559744
#> MZ SE                     0.124933               0.125078
#>             cause1:strata4:spline4 cause1:strata4:spline5
#> DZ estimate               -4.82018               -3.87710
#> DZ SE                      0.13709                0.08333
#> MZ estimate               -5.39591               -4.43014
#> MZ SE                      0.22881                0.14030
#>             cause1:traject:countryDenmark cause1:traject:countryFinland
#> DZ estimate                      2.470087                      2.634783
#> DZ SE                            0.055722                      0.085629
#> MZ estimate                      2.805731                      2.906389
#> MZ SE                            0.092889                      0.148814
#>             cause1:traject:countryNorway cause1:traject:countrySweden
#> DZ estimate                     2.570084                     2.910946
#> DZ SE                           0.077549                     0.073602
#> MZ estimate                     3.202949                     3.292289
#> MZ SE                           0.148409                     0.119783
#>             cause2:strata1:spline1 cause2:strata1:spline2
#> DZ estimate                -1.9009               -2.52802
#> DZ SE                       0.2480                0.28503
#> MZ estimate                -3.1012               -3.24816
#> MZ SE                       0.4658                0.42813
#>             cause2:strata1:spline3 cause2:strata1:spline4
#> DZ estimate               -2.91962               -4.74451
#> DZ SE                      0.32786                0.49056
#> MZ estimate               -4.59427               -6.39337
#> MZ SE                      0.64848                0.77990
#>             cause2:strata1:spline5 cause2:strata2:spline1
#> DZ estimate               -3.24579               -1.66595
#> DZ SE                      0.36049                0.21764
#> MZ estimate               -5.35583               -2.30082
#> MZ SE                      0.64279                0.42532
#>             cause2:strata2:spline2 cause2:strata2:spline3
#> DZ estimate               -1.95962               -2.91715
#> DZ SE                      0.24609                0.31500
#> MZ estimate               -2.99934               -3.42842
#> MZ SE                      0.46611                0.53532
#>             cause2:strata2:spline4 cause2:strata2:spline5
#> DZ estimate               -4.64212               -3.72697
#> DZ SE                      0.47831                0.43403
#> MZ estimate               -5.91248               -3.92354
#> MZ SE                      0.85588                0.56325
#>             cause2:strata3:spline1 cause2:strata3:spline2
#> DZ estimate               -2.10036               -2.33815
#> DZ SE                      0.27427                0.27205
#> MZ estimate               -2.92932               -3.01743
#> MZ SE                      0.52023                0.51301
#>             cause2:strata3:spline3 cause2:strata3:spline4
#> DZ estimate               -3.00447                -4.5611
#> DZ SE                      0.40179                 0.5193
#> MZ estimate               -4.32479                -6.3118
#> MZ SE                      0.75126                 1.0680
#>             cause2:strata3:spline5 cause2:strata4:spline1
#> DZ estimate               -3.57013               -1.87097
#> DZ SE                      0.45869                0.18126
#> MZ estimate               -4.93592               -2.37939
#> MZ SE                      0.71230                0.27319
#>             cause2:strata4:spline2 cause2:strata4:spline3
#> DZ estimate               -2.36322               -3.06084
#> DZ SE                      0.20950                0.27425
#> MZ estimate               -2.91045               -3.84355
#> MZ SE                      0.29584                0.43432
#>             cause2:strata4:spline4 cause2:strata4:spline5
#> DZ estimate               -4.90966               -3.52848
#> DZ SE                      0.41990                0.31142
#> MZ estimate               -6.20353               -4.89571
#> MZ SE                      0.63417                0.53922
#>             cause2:traject:countryDenmark cause2:traject:countryFinland
#> DZ estimate                       2.81623                       2.71478
#> DZ SE                             0.28770                       0.27170
#> MZ estimate                       3.83946                       3.78564
#> MZ SE                             0.58437                       0.54292
#>             cause2:traject:countryNorway cause2:traject:countrySweden
#> DZ estimate                      2.79434                      2.82169
#> DZ SE                            0.29131                      0.24524
#> MZ estimate                      4.11384                      3.49492
#> MZ SE                            0.58977                      0.41991
#>             vcov:risk1:risk1 vcov:risk1:risk2 vcov:risk2:risk2
#> DZ estimate       -0.0075784          0.46189         0.028412
#> DZ SE              0.1788157          0.31467         0.174085
#> MZ estimate        0.3542122          2.96485         0.585914
#> MZ SE              0.2655416          1.06490         0.377501
#>             vcov:risk1:traject1 vcov:risk2:traject1 vcov:traject1:traject1
#> DZ estimate           -0.393206           -0.068591                -3.3268
#> DZ SE                  0.046595            0.117840                 7.6790
#> MZ estimate           -0.495516            0.460557                -4.2756
#> MZ SE                  0.202764            0.227532                 1.7806
#>             vcov:risk1:traject2 vcov:risk2:traject2 vcov:traject1:traject2
#> DZ estimate            -0.23010           -0.406945              -0.055455
#> DZ SE                   0.17825            0.170872               1.172040
#> MZ estimate            -0.38416           -0.060491               0.929293
#> MZ SE                   0.30470            0.278583               0.204267
#>             vcov:traject2:traject2
#> DZ estimate              -2.228813
#> DZ SE                     1.742192
#> MZ estimate              -2.390810
#> MZ SE                     0.098484

# the estimated covariance matrices along with standard errors
fit_dz$vcov_est # the estimates for DZ
#>         [,1]   [,2]     [,3]     [,4]
#> [1,] 0.98496 0.4584 -0.39024 -0.22837
#> [2,]      NA 1.2718 -0.25218 -0.52495
#> [3,]      NA     NA  0.16060  0.11640
#> [4,]      NA     NA       NA  0.23322
fit_dz$vcov_SE # the SEs
#>         [,1]    [,2]     [,3]    [,4]
#> [1,] 0.35225 0.71359 0.146376 0.36745
#> [2,]      NA 0.53446 0.197064 0.37579
#> [3,]      NA      NA 0.030827 0.13285
#> [4,]      NA      NA       NA 0.18775
fit_mz$vcov_est # the estimates for MZ
#>        [,1]    [,2]     [,3]     [,4]
#> [1,] 2.0308  4.2251 -0.70614 -0.54745
#> [2,]     NA 12.0182 -0.64168 -1.24765
#> [3,]     NA      NA  0.45784  0.17542
#> [4,]     NA      NA       NA  1.02320
fit_mz$vcov_SE # the SEs
#>        [,1]   [,2]     [,3]    [,4]
#> [1,] 1.0785 4.9645 0.446623 0.91742
#> [2,]     NA 4.4888 0.779556 2.00501
#> [3,]     NA     NA 0.067354 0.34584
#> [4,]     NA     NA       NA 0.34386
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
