# Methodology Notes
## Project 2: Multisource SAE for SDG Monitoring in Niger

---

## 1. The Multisource SAE Problem

Standard Fay-Herriot models use a single source of auxiliary information 
(typically census data). This project extends the framework to incorporate 
**multiple heterogeneous data sources** as auxiliary variables, motivated by 
the increasing availability of:
- Administrative registers (INS Niger, ministries)
- Earth observation data (satellite imagery)
- Web-scraped price indices (FAO GIEWS)
- Mobile phone connectivity data

---

## 2. Non-Probability Data Bias Correction

Non-probability data (web-scraped, satellite, administrative) cannot be 
treated as representative samples. We apply:

### 2a. Propensity Score Weighting (Chen et al., 2020)
Estimate the probability of a unit being included in the non-probability 
source using a logistic regression model:
```
logit(P(included_i)) = alpha + gamma * X_i
```
Then reweight non-probability estimates by the inverse propensity scores.

### 2b. Doubly Robust Estimation
Combine propensity weighting with outcome regression for robustness:
```
theta_DR = theta_OR + (1/n) * sum(w_i * (Y_i - mu_i(X_i)))
```

---

## 3. Extended Fay-Herriot Model with Multisource Auxiliaries

**Standard FH model:**
```
y_i = x_i'β + v_i + e_i
```

**Extended multisource FH model:**
```
y_i = β_0 + β_1*x_census_i + β_2*x_admin_i + β_3*x_satellite_i + v_i + e_i
```

Model selection uses **AIC** and **cross-validation MSE** to identify the 
optimal set of auxiliary variables.

---

## 4. Differential Privacy for Official Statistics

### 4a. (ε, δ)-Differential Privacy
A mechanism M satisfies (ε, δ)-DP if for all adjacent databases D, D':
```
P(M(D) ∈ S) ≤ exp(ε) * P(M(D') ∈ S) + δ
```

### 4b. Gaussian Mechanism
For a query f with L2-sensitivity Δf:
```
M(D) = f(D) + N(0, σ²I)
```
where σ = Δf * sqrt(2*ln(1.25/δ)) / ε

### 4c. Calibration for SAE
For regional poverty rates with n_i households:
- Global sensitivity Δf = 1/n_i
- Recommended: ε = 2.0, δ = 1e-5
- This achieves < 2pp mean absolute error for regions with n ≥ 160

---

## References

1. Fay, R.E. & Herriot, R.A. (1979). Estimates of income for small places. *JASA* 74(366), 269-277.
2. Molina, I. (2022). Disaggregating data in household surveys. *ECLAC Statistics Series No. 97*.
3. Rao, J.N.K. & Molina, I. (2015). *Small Area Estimation* (2nd ed.). Wiley.
4. Chen, Y. et al. (2020). Doubly robust inference for non-probability samples. *JRSS-B* 82(2).
5. Dwork, C. & Roth, A. (2014). *Algorithmic Foundations of Differential Privacy*. Foundations & Trends.
6. Enamorado, T. et al. (2019). Using a probabilistic model to assist merging of large-scale datasets. *Political Analysis* 27(2).
7. UNSD (2020). *Handbook on Small Area Estimation*. United Nations Statistics Division.
