# NNRPanel: Nonlinear Panel Estimation with Interactive Fixed Effects

## Description

The **NNRPanel** package implements nonlinear panel estimators with interactive fixed effects.  
The proposed method involves two steps: In the first step, we convexify the optimization problem using nuclear norm regularization (NNR) and obtain preliminary NNR estimators of the parameters, including the fixed effects. Then, we find the global solution of the original optimization problem using a standard gradient descent method initialized at these preliminary estimates. 

The proposed estimator avoids solving a high-dimensional non-convex optimization problem and can be feasibly computed in large nonlinear panels.   We also provide analytical and sample-splitting bias corrections.

---

## Details

NNRPanel is designed for panel data, and support for network data will be added soon. The package supports both logistic and Poisson models.   It also provides a data-driven tuning-parameter selection method and a simple procedure for determining the number of factors.

### Key Functionality

- Data-driven tuning-parameter selection used in the first-step estimation.  
- Nuclear-norm–regularized first-step estimation of nonlinear panel models.  
- Data-driven selection of the number of latent factors.  
- Second-step local MLE estimation using first-step estimators as initial values.  
- Analytical and sample-splitting bias corrections for high-dimensional panel models.

### Main Exported Functions

- `NNRPanel_estimate()` — Main estimator pipeline (NNR → local MLE → bias correction).  
- `calculate_tuning_parameter()` — Data-driven tuning parameter φ.  
- `determine_num_factor()` — Estimate the factor dimension.

This package is intended for econometric applications involving nonlinear panel models with interactive fixed effects when the number of individuals and/or time periods is large.

---

## Author

**Weisheng Zhang**  
Maintainer: *Weisheng Zhang* <wszhang.econ@gmail.com>

---

To install this package in R, run the following commands:

```R
install.packages("devtools")
library(devtools) 
install_github("wszhang-econ/NNRPanel")
```



## Examples

```r
library(NNRPanel)
# Simulated example (Logit)
data_list <- gen_data_logit(N = 200, T = 200, num_factor = 2, num_Z = 1, beta = c(0.2))
data_frame <- list_to_data_frame(data_list)

# Main estimator with data-driven factor selection
res <- NNRPanel_estimate(data_frame, func = "logit", R_max = 5, R_true = 2, delta = 0.05)
```
#### References
Andrei Zeleneev, Weisheng Zhang (2025). <b>Tractable Estimation of Nonlinear Panels with Interactive Fixed Effects</b> [<a href="https://arxiv.org/abs/2511.15427">link</a>]
