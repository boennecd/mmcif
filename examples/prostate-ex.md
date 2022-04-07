
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
#>      359      357
fit_mz$fit$counts
#> function gradient 
#>      475      473

sqrt(sum(fit_dz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 1.1649
sqrt(sum(fit_mz$gr_mle^2)) # L2 norm of the gradient at the MLE
#> [1] 0.6469

# the maximum log composite likelihood
print(-fit_dz$fit$value, digits = 10)
#> [1] -23833.7303
print(-fit_mz$fit$value, digits = 10)
#> [1] -12300.04458

# the time to estimate the model and to compute the sandwich estimator
with(fit_dz$time, start + fit + sandwich)
#>     user   system  elapsed 
#> 7063.236    0.264 1779.904
with(fit_mz$time, start + fit + sandwich)
#>     user   system  elapsed 
#> 6759.209    0.196 1708.102

# the parameter estimates along with standard errors
rbind(`DZ estimate` = fit_dz$fit$par,
      `DZ SE` = diag(fit_dz$sandwich_est) |> sqrt(),
      `MZ estimate` = fit_mz$fit$par,
      `MZ SE` = diag(fit_mz$sandwich_est) |> sqrt())
#>             cause1:risk:countryDenmark cause1:risk:countryFinland
#> DZ estimate                   0.725230                   0.638287
#> DZ SE                         0.068232                   0.091785
#> MZ estimate                   0.657855                   0.645005
#> MZ SE                         0.116359                   0.169977
#>             cause1:risk:countryNorway cause1:risk:countrySweden
#> DZ estimate                   0.42635                  0.433080
#> DZ SE                         0.10127                  0.057892
#> MZ estimate                   0.35831                  0.236878
#> MZ SE                         0.15352                  0.096164
#>             cause2:risk:countryDenmark cause2:risk:countryFinland
#> DZ estimate                   -2.45909                   -1.47770
#> DZ SE                          0.20479                    0.18491
#> MZ estimate                   -3.01983                   -2.30931
#> MZ SE                          0.37274                    0.39365
#>             cause2:risk:countryNorway cause2:risk:countrySweden
#> DZ estimate                  -1.88167                  -1.57149
#> DZ SE                         0.21183                   0.14942
#> MZ estimate                  -2.60405                  -2.19200
#> MZ SE                         0.38103                   0.27966
#>             cause1:strata1:spline1 cause1:strata1:spline2
#> DZ estimate              -1.752914              -2.343784
#> DZ SE                     0.067256               0.076519
#> MZ estimate              -1.756715              -2.276234
#> MZ SE                     0.100087               0.109530
#>             cause1:strata1:spline3 cause1:strata1:spline4
#> DZ estimate              -3.109840               -4.79175
#> DZ SE                     0.094965                0.13455
#> MZ estimate              -3.167480               -4.61039
#> MZ SE                     0.135408                0.18271
#>             cause1:strata1:spline5 cause1:strata2:spline1
#> DZ estimate               -3.96851              -1.938246
#> DZ SE                      0.11883               0.098653
#> MZ estimate               -3.71461              -2.088730
#> MZ SE                      0.14481               0.157625
#>             cause1:strata2:spline2 cause1:strata2:spline3
#> DZ estimate               -2.54238               -3.10682
#> DZ SE                      0.11026                0.11985
#> MZ estimate               -2.60701               -3.09055
#> MZ SE                      0.17220                0.19316
#>             cause1:strata2:spline4 cause1:strata2:spline5
#> DZ estimate               -5.03750               -3.97240
#> DZ SE                      0.20747                0.16270
#> MZ estimate               -5.24198               -4.10931
#> MZ SE                      0.31773                0.19692
#>             cause1:strata3:spline1 cause1:strata3:spline2
#> DZ estimate              -1.649655               -2.25624
#> DZ SE                     0.096054                0.10475
#> MZ estimate              -1.978733               -2.54715
#> MZ SE                     0.149497                0.16511
#>             cause1:strata3:spline3 cause1:strata3:spline4
#> DZ estimate               -2.92399               -4.66918
#> DZ SE                      0.13249                0.17358
#> MZ estimate               -3.18371               -5.37515
#> MZ SE                      0.18110                0.30375
#>             cause1:strata3:spline5 cause1:strata4:spline1
#> DZ estimate               -3.77441              -1.920475
#> DZ SE                      0.14551               0.085667
#> MZ estimate               -4.39703              -2.007745
#> MZ SE                      0.21358               0.123171
#>             cause1:strata4:spline2 cause1:strata4:spline3
#> DZ estimate              -2.575321               -3.20190
#> DZ SE                     0.094136                0.09658
#> MZ estimate              -2.695833               -3.35564
#> MZ SE                     0.133724                0.13346
#>             cause1:strata4:spline4 cause1:strata4:spline5
#> DZ estimate               -5.04805               -4.02642
#> DZ SE                      0.17045                0.10663
#> MZ estimate               -5.17019               -4.14927
#> MZ SE                      0.24662                0.15125
#>             cause1:traject:countryDenmark cause1:traject:countryFinland
#> DZ estimate                      2.468395                      2.657489
#> DZ SE                            0.066514                      0.099995
#> MZ estimate                      2.650195                      2.782664
#> MZ SE                            0.107235                      0.161250
#>             cause1:traject:countryNorway cause1:traject:countrySweden
#> DZ estimate                     2.554338                     2.914434
#> DZ SE                           0.088952                     0.089143
#> MZ estimate                     3.041019                     3.123474
#> MZ SE                           0.158320                     0.138948
#>             cause2:strata1:spline1 cause2:strata1:spline2
#> DZ estimate               -2.16289               -3.04518
#> DZ SE                      0.35398                0.44359
#> MZ estimate               -3.00694               -3.48983
#> MZ SE                      0.49142                0.48347
#>             cause2:strata1:spline3 cause2:strata1:spline4
#> DZ estimate               -3.26427               -5.93010
#> DZ SE                      0.48061                0.80541
#> MZ estimate               -4.21176               -6.58762
#> MZ SE                      0.60667                0.77045
#>             cause2:strata1:spline5 cause2:strata2:spline1
#> DZ estimate               -4.47311               -1.90904
#> DZ SE                      0.61067                0.29858
#> MZ estimate               -5.77241               -2.62809
#> MZ SE                      0.78719                0.45102
#>             cause2:strata2:spline2 cause2:strata2:spline3
#> DZ estimate               -2.16150               -3.10913
#> DZ SE                      0.32135                0.45979
#> MZ estimate               -3.06565               -4.33862
#> MZ SE                      0.45114                0.67316
#>             cause2:strata2:spline4 cause2:strata2:spline5
#> DZ estimate               -5.10930                -4.2497
#> DZ SE                      0.67031                 0.5755
#> MZ estimate               -6.76606                -4.5538
#> MZ SE                      0.90730                 0.6351
#>             cause2:strata3:spline1 cause2:strata3:spline2
#> DZ estimate               -2.42546               -2.73355
#> DZ SE                      0.38676                0.40958
#> MZ estimate               -3.04777               -3.20868
#> MZ SE                      0.52732                0.53719
#>             cause2:strata3:spline3 cause2:strata3:spline4
#> DZ estimate               -3.64299               -5.33707
#> DZ SE                      0.52176                0.70684
#> MZ estimate               -4.55254               -6.81696
#> MZ SE                      0.62283                1.03563
#>             cause2:strata3:spline5 cause2:strata4:spline1
#> DZ estimate               -4.14271               -2.06304
#> DZ SE                      0.57723                0.27266
#> MZ estimate               -5.59379               -2.48165
#> MZ SE                      0.81264                0.28531
#>             cause2:strata4:spline2 cause2:strata4:spline3
#> DZ estimate               -2.68398               -3.23349
#> DZ SE                      0.32621                0.40012
#> MZ estimate               -3.09836               -3.72756
#> MZ SE                      0.32534                0.38780
#>             cause2:strata4:spline4 cause2:strata4:spline5
#> DZ estimate               -5.39106               -3.86360
#> DZ SE                      0.65828                0.47619
#> MZ estimate               -6.24782               -4.74432
#> MZ SE                      0.60421                0.49547
#>             cause2:traject:countryDenmark cause2:traject:countryFinland
#> DZ estimate                        3.0640                       2.86208
#> DZ SE                              0.4255                       0.38934
#> MZ estimate                        4.7352                       4.49721
#> MZ SE                              0.5778                       0.53223
#>             cause2:traject:countryNorway cause2:traject:countrySweden
#> DZ estimate                      3.02877                      3.04929
#> DZ SE                            0.42107                      0.37468
#> MZ estimate                      4.82667                      4.26301
#> MZ SE                            0.58886                      0.42694
#>             vcov:risk1:risk1 vcov:risk1:risk2 vcov:risk2:risk2
#> DZ estimate         -0.35748          0.22323         0.095034
#> DZ SE                0.14465          0.28445         0.162007
#> MZ estimate          0.19849          0.67629         0.711584
#> MZ SE                0.11687          0.52666         0.129090
#>             vcov:risk1:traject1 vcov:risk2:traject1 vcov:traject1:traject1
#> DZ estimate           -0.249503           -0.050647               -0.95200
#> DZ SE                  0.084676            0.116358                0.22319
#> MZ estimate           -0.527867           -0.068836               -4.44259
#> MZ SE                  0.071126            0.249240                1.44251
#>             vcov:risk1:traject2 vcov:risk2:traject2 vcov:traject1:traject2
#> DZ estimate            -0.37174            -0.40048               -0.35432
#> DZ SE                   0.22991             0.19691                0.28944
#> MZ estimate            -0.14517            -0.90046               -0.74747
#> MZ SE                   0.32477             0.23768                0.28825
#>             vcov:traject2:traject2
#> DZ estimate               -2.45898
#> DZ SE                      1.75689
#> MZ estimate               -2.35657
#> MZ SE                      0.53165

# the estimated covariance matrices along with standard errors
fit_dz$vcov_est # the estimates for DZ
#>         [,1]    [,2]     [,3]      [,4]
#> [1,] 0.48922 0.15613 -0.17451 -0.260008
#> [2,]      NA 1.25916 -0.11139 -0.523387
#> [3,]      NA      NA  0.21379 -0.023724
#> [4,]      NA      NA       NA  0.431429
fit_dz$vcov_SE # the SEs
#>         [,1]    [,2]     [,3]    [,4]
#> [1,] 0.14153 0.41184 0.110076 0.32911
#> [2,]      NA 0.43827 0.211346 0.40909
#> [3,]      NA      NA 0.046865 0.19367
#> [4,]      NA      NA       NA 0.33544
fit_mz$vcov_est # the estimates for MZ
#>        [,1]    [,2]     [,3]     [,4]
#> [1,] 1.4873 0.82477 -0.64377 -0.17704
#> [2,]     NA 4.60761 -0.49722 -1.93260
#> [3,]     NA      NA  0.28352  0.12982
#> [4,]     NA      NA       NA  1.39959
fit_mz$vcov_SE # the SEs
#>         [,1]   [,2]     [,3]    [,4]
#> [1,] 0.34765 1.3157 0.175995 0.79028
#> [2,]      NA 1.3855 0.613378 1.17683
#> [3,]      NA     NA 0.073669 0.32462
#> [4,]      NA     NA       NA 0.44785
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
