``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS  10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.0.2  magrittr_2.0.2  fastmap_1.1.0   cli_3.1.1      
    ##  [5] tools_4.0.2     htmltools_0.5.2 rstudioapi_0.13 yaml_2.2.2     
    ##  [9] stringi_1.7.6   rmarkdown_2.10  knitr_1.37      stringr_1.4.0  
    ## [13] xfun_0.29       digest_0.6.29   rlang_1.0.1     evaluate_0.14

``` r
library(tidyverse)
library(latex2exp)
source("functions.R")
```

## Best ITH association study design

This document reproduces some parts of numerical studies we performed
for the best ITH association study design paper. Here is the outline:

1.  Generate multiregion genomic datasets.
2.  Estimate parameters using generated datasets
3.  Given estimated parameters, find the optimal design for
    - pre-collection scenario
    - post-collection scenario

### Generate multiregion genomic datasets

We generate multiregion genomic datasets to estimate three parameters,
$\tau^2, \sigma^2$ and $\rho$ which will be used in objective functions
later. The function `simulate_tumors()` which is in `functions.R`
generates a list of matrix where each of the matrix represent a
multiregion genomic dataset of a tumor. Important inputs the function
require are:

- number of tumors (`nTumors`)
- number of samples for each tumor (`nSamp`)
- number of genes or probes or copy number variations or others (`nSeg`)
- underlying mutation rate of a patient (`mutationRates`) which will be
  used to generate underlying mutation status vector
- lower bound of the uniform distribution from which $\theta$ will be
  drawn
- upper bound of the uniform distribution from which ITH controlling
  parameter, $\theta$, will be drawn

$\theta$ controls ITH with $\theta = 1$ rendering no ITH, and
$\theta = 0.5$ rendering maximum ITH. The following code chunk generates
multiregion genomic datasets of 50 tumors. The columns represent samples
and each row represents a mutation status of a gene across the samples.

``` r
nTumors <- 50
nSamp <- 10
nSeg <- 40
theta_lb <- 0.5
theta_ub <- 0.6
mutationRates <- rbeta(nSeg,2,3)

set.seed(1)
tumor_mat_list <- simulate_tumors(nTumors=nTumors,nSamp=nSamp,nSeg=nSeg,
                                  mutationRates=mutationRates,
                                  theta_lb=theta_lb,theta_ub=theta_ub)
tumor_mat_list[[1]]
```

    ##       s1 s2 s3 s4 s5 s6 s7 s8 s9 s10
    ##  [1,]  0  0  0  0  0  0  0  0  0   0
    ##  [2,]  0  0  0  0  0  0  0  0  0   0
    ##  [3,]  0  0  0  0  0  0  0  0  0   0
    ##  [4,]  0  1  0  0  1  1  0  0  1   0
    ##  [5,]  0  0  0  0  0  0  0  0  0   0
    ##  [6,]  0  1  0  1  0  1  0  1  1   0
    ##  [7,]  1  1  0  1  1  0  1  1  1   0
    ##  [8,]  1  1  1  0  0  1  0  1  1   1
    ##  [9,]  1  1  1  0  1  0  1  1  1   1
    ## [10,]  0  0  0  0  0  0  0  0  0   0
    ## [11,]  0  0  0  0  0  0  0  0  0   0
    ## [12,]  0  0  0  0  0  0  0  0  0   0
    ## [13,]  0  0  0  0  0  0  0  0  0   0
    ## [14,]  1  1  1  0  0  1  1  0  1   1
    ## [15,]  0  0  1  1  0  0  1  1  0   1
    ## [16,]  0  0  0  0  0  0  0  0  0   0
    ## [17,]  0  0  1  0  1  1  0  0  0   1
    ## [18,]  0  1  0  1  0  0  0  1  1   1
    ## [19,]  0  0  0  0  0  0  0  0  0   0
    ## [20,]  1  0  0  0  1  0  1  0  1   1
    ## [21,]  0  1  0  0  0  1  1  1  1   0
    ## [22,]  0  0  0  0  0  0  0  0  0   0
    ## [23,]  1  1  1  0  1  1  1  0  1   0
    ## [24,]  0  0  0  0  0  0  0  0  0   0
    ## [25,]  0  0  0  0  0  0  0  0  0   0
    ## [26,]  0  0  0  0  0  0  0  0  0   0
    ## [27,]  0  0  0  0  0  0  0  0  0   0
    ## [28,]  0  0  1  1  0  0  1  1  1   1
    ## [29,]  1  1  1  0  0  1  1  0  1   0
    ## [30,]  0  0  0  0  0  0  0  0  0   0
    ## [31,]  0  0  0  0  0  0  0  0  0   0
    ## [32,]  0  0  0  0  0  0  0  0  0   0
    ## [33,]  0  1  1  1  1  0  0  1  1   0
    ## [34,]  0  0  0  0  0  0  0  0  0   0
    ## [35,]  1  1  0  1  1  1  0  0  1   0
    ## [36,]  0  0  0  0  0  0  0  0  0   0
    ## [37,]  1  1  0  1  0  1  1  0  1   0
    ## [38,]  0  0  0  0  0  0  0  0  0   0
    ## [39,]  1  1  1  0  1  0  1  0  0   1
    ## [40,]  0  0  0  0  0  0  0  0  0   0

### Estimate parameters

We use the generated dataset above to illustrate how the parameters
$\tau^2, \sigma^2$ and $\rho$ can be estimated. We define APITH,
$\hat{D}_n$, and its conditional variance,
$\textrm{Var}(\hat{D}_n|\mu_n)$ where $\mu_n$ is the true expected ITH
of a subject $n$:

The function, `estimate_parameters()` which is defined in `functions.R`
takes a list of multiregion genomic profile matrices, representing a
single study, and estimates the parameters..

``` r
estimate_parameters(tumor_mat_list)
```

    ## sigma_square          rho   tau_square 
    ##    4.1155556    0.1244444    2.0723978

### Find the optimal design

The objective function we want to maximize is

$$
\varphi := \sum_{n=1}^N \frac{1}{\tau^2 + f(K_n;\sigma^2,\rho)}
$$

where $N$ is the total number of subjects, $K_n$ is the number of
samples for subject $n$ and

$$
f(K_n;\sigma^2,\rho) := \frac{2}{K_n(K_n - 1)}\sigma^2 + \frac{4(K_n - 2)}{K_n(K_n-1)}\rho.
$$ Given the budget to profile $M$ tumor samples at most, we want to
find $(N,K_1,...,K_n)$ that maximize $\varphi$ under the constraints
$K_n \geq 2$ and $\sum_{i=1}^N K_i \leq M.$

#### Pre-collection scenario

Pre-collection scenario assumes no samples have been collected. Thus,
one has freedom to select any $N$ and $(K_1,...,K_N)$ to maximize
$\varphi$. For a patient with $K$ tumor samples, we quantify the
contribution of each tumor sample to $\varphi$ as

$$
\xi(K) = \frac{1}{K}\frac{1}{\tau^2 + f(K;\sigma^2,\rho)}.
$$

Assume $K_{max}$ is the integer that maximizes $\xi(K)$ so that
$\varphi_{max} = M\cdot\xi(K_{max})$. To find $K_{max}$, we examine
$\varphi$ at some range of values of $K$.
In the following example, we work with different sets of estimated
parameters $\tau^2, \sigma^2$ and $\rho$ which were computed using tumor
datasets comprised of 200 tumors, each of them with 30 samples. We
compute the optimal number of samples $K_{max}$ as well as
$\varphi_{max} = M\cdot\xi(K_{max})$.
Import parameter data

``` r
estimated_parameters_tab <- openxlsx::read.xlsx("estimated_parameters.xlsx")
estimated_parameters_tab
```

    ##    a b    l    u tau_sq sigma_sq   rho
    ## 1  5 2 0.55 0.60  2.009    6.800 0.147
    ## 2  5 2 0.50 0.55  2.093    6.880 0.012
    ## 3  5 2 0.60 0.65  1.939    6.754 0.393
    ## 4  5 2 0.65 0.70  1.752    6.696 0.707
    ## 5  5 2 0.75 0.80  1.419    6.186 1.409
    ## 6  5 2 0.60 0.70  2.031    6.698 0.537
    ## 7  5 2 0.50 0.70  2.355    6.749 0.301
    ## 8  5 2 0.60 0.80  3.123    6.527 0.877
    ## 9  5 2 0.50 0.80  3.788    6.609 0.600
    ## 10 5 2 0.80 0.90  2.036    5.242 1.685
    ## 11 5 2 0.70 0.90  4.718    5.831 1.463
    ## 12 5 2 0.60 0.90  6.935    6.133 1.158
    ## 13 5 2 0.50 0.90  7.933    6.244 0.864

We set $M = 100$ and search $K_{max} \in \{2,3,...,10\}$. `phi1()`
computes $\varphi$ given the parameters, a total budget $N$ and number
of samples $K$ as inputs.

``` r
M <- 100
nSamp_max <- 10
res <- apply(estimated_parameters_tab,1,function(x){
    tau_sq <- x[5]
    sigma_sq <- x[6]
    rho <- x[7]
    phi <- lapply(2:nSamp_max,function(nSamp) phi1(sigma_sq,rho,tau_sq,nSamp,M)) %>% unlist()
    K_max <- which.max(phi) + 1
    phi_max <- phi[which.max(phi)]
    c(K_max, phi_max)
}) %>% t()
colnames(res) <- c("K_max","phi_max")
res <- as_tibble(cbind(estimated_parameters_tab,res))
res
```

    ## # A tibble: 13 × 9
    ##        a     b     l     u tau_sq sigma_sq   rho K_max phi_max
    ##    <dbl> <dbl> <dbl> <dbl>  <dbl>    <dbl> <dbl> <dbl>   <dbl>
    ##  1     5     2  0.55  0.6    2.01     6.8  0.147     4    7.72
    ##  2     5     2  0.5   0.55   2.09     6.88 0.012     4    7.70
    ##  3     5     2  0.6   0.65   1.94     6.75 0.393     4    7.52
    ##  4     5     2  0.65  0.7    1.75     6.70 0.707     4    7.49
    ##  5     5     2  0.75  0.8    1.42     6.19 1.41      3    7.54
    ##  6     5     2  0.6   0.7    2.03     6.70 0.537     3    7.21
    ##  7     5     2  0.5   0.7    2.36     6.75 0.301     3    6.94
    ##  8     5     2  0.6   0.8    3.12     6.53 0.877     3    5.67
    ##  9     5     2  0.5   0.8    3.79     6.61 0.6       3    5.22
    ## 10     5     2  0.8   0.9    2.04     5.24 1.68      2    6.87
    ## 11     5     2  0.7   0.9    4.72     5.83 1.46      2    4.74
    ## 12     5     2  0.6   0.9    6.94     6.13 1.16      2    3.83
    ## 13     5     2  0.5   0.9    7.93     6.24 0.864     2    3.53

We plot $\varphi$ as a function of $K$ for selected parameters:

- $(\hat{\tau}^2,\hat{\sigma}^2,\hat{\rho}) = (2.04,5.24,1.68)$
- $(\hat{\tau}^2,\hat{\sigma}^2,\hat{\rho}) = (1.42,6.19,1.41)$
- $(\hat{\tau}^2,\hat{\sigma}^2,\hat{\rho}) = (2.01,6.80,0.15)$

``` r
res_selected <- res[c(10,5,1),]
par(mfrow=c(1,3))
par(oma = c(3,3,0,0))
par(mar = c(2,2,2,1))
apply(res_selected,1,function(x){
    tau_sq <- x[5]
    sigma_sq <- x[6]
    rho <- x[7]
    K <- 2:10
    phi <- lapply(2:nSamp_max,function(nSamp) phi1(sigma_sq,rho,tau_sq,nSamp,M)) %>% unlist()
    fit <- lm(phi~poly(K,6,raw=F))
    main_str <- TeX(sprintf("$(\\tau^2,\\sigma^2,\\rho) = (%0.2f,%0.2f,%0.2f)$",tau_sq,sigma_sq,rho))
    plot(K,phi,ylab="", xlab="",main = main_str,col=ifelse(phi==max(phi),"red","black"),
         pch=20,ylim = c(3,8))
    lines(K, predict(fit,data.frame(x=K)), col="blue", lwd = 0.5)
})
```

    ## NULL

``` r
mtext(TeX("$K$"), side = 1, outer = T, line = 1)
mtext(TeX("$\\varphi$"), side = 2, outer = T, line = 0.75,las=1)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

#### Post-collection scenario

For the post-collection scenario, we assume tumor samples have been
collected for $L$ patients, each of which has $m_l$ tumor samples
($l = 1,...,L$). We represent this sample collection scheme in the
following way:

| number of tumor samples collected | number of patients |
|:---------------------------------:|:------------------:|
|                 2                 |       $n_2$        |
|                 3                 |       $n_3$        |
|                 4                 |       $n_4$        |
|             $\cdots$              |      $\cdots$      |
|                $A$                |       $n_A$        |

$n_a$ is the number of patients with $a (2 \leq a \leq A)$ tumor samples
available. We let $\omega_a \geq 0$ be the number of patients who
contribute $a$ tumor samples to the study. We want to find
$(\omega_2,...,\omega_A$) meeting appropriate constraints and that
achieves $$
\varphi_{max} := \max_{\omega_2,...,\omega_A} \sum_{a=2}^A \frac{\omega_a}{\tau^2 + f(a,\sigma^2,\rho)}
$$ where $f(a,\sigma^2,\rho)$ is the same as what we have defined in the
pre-collection scenario. The details of constraints can be found in the
paper.
For illustration, we explore the optimal design for a study with the
following estimated parameters:

- $(\hat{\tau}^2,\hat{\sigma}^2,\hat{\rho}) = (2.04,5.24,1.68)$

We assume we have already collected 372 tumor samples from 84 patients
as the following:

| number of tumor samples collected | number of patients |
|:---------------------------------:|:------------------:|
|                 2                 |         20         |
|                 3                 |         16         |
|                 4                 |         14         |
|                 5                 |         10         |
|                 6                 |         8          |
|                 7                 |         6          |
|                 8                 |         4          |
|                 9                 |         4          |
|                10                 |         2          |

The function `recursive_search` that computes the optimal design takes
the following as inputs:

- `A`: the largest number of tumor samples collected among all subjects
- `M`: a total number of tumor samples budgeted
- `R`: a number of available subjects with more than `A` samples
- `Om`: a vector of “candidate” number of subjects
  ($(\omega_1,\omega_2,...,\omega_A)$)
- `N`: a vector of number of subjects collected for each number of
  samples ($(n_1,n_2,...,n_A)$)
- `tau_sq`: $\tau^2$
- `sigma_sq`: $\sigma^2$
- `rho`: $\rho$

We obtain optimal design solution for $M = 200.$

``` r
params <- c(2.04,5.24,1.68)
A <- 10
N <- c(0,20,16,14,10,8,6,4,4,2)
Om <- rep(0,10)
R <- 0
M <- 200

res <- recursive_search(A = A,M = M,R = R, Om = Om, N = N,
                        tau_sq = params[1], sigma_sq = params[2],rho = params[3])
```

The output of the function is a list of two things:

- The optimal design solution $(\omega_1,\omega_2,...,\omega_A)$ that
  achieves $\varphi_{max}$
- $\varphi_{max}$

``` r
res
```

    ## [[1]]
    ##  [1]  0 52 32  0  0  0  0  0  0  0
    ## 
    ## [[2]]
    ## [1] 13.6646
