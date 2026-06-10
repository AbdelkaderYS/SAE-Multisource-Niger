# Methodology Notes
## Project 2: Multisource SAE for SDG Monitoring in Niger

---

## 1. The Multisource SAE Problem

Standard Fay-Herriot models rely on a single source of auxiliary information,
typically census data. This project extends the framework to incorporate
multiple heterogeneous data sources as auxiliary variables:

- Census data (INS Niger RGPH 2012)
- Administrative registers (INS Niger: health facilities, education)
- Earth observation data (MODIS NDVI, VIIRS nighttime lights)
- Open national statistics (World Bank Open Data API)

---

## 2. Probabilistic Record Linkage

When two data sources share no common identifier, probabilistic record
linkage estimates match probabilities from observable quasi-identifiers.

### Fellegi-Sunter Model (Enamorado et al., 2019)

For a pair (a, b), the posterior match probability is:
```
P(match | gamma) = P(gamma | match) * m / [P(gamma | match)*m + P(gamma | non-match)*(1-m)]
```

where gamma is the agreement pattern vector and m is the prior match rate.
Parameters are estimated via the EM algorithm. Pairs with posterior
probability above a threshold are declared matches.

**Implementation:** NIHR61FL (Household Recode) x NIIR61FL (Individual
Women Recode), linked on region, urban/rural, and education level.
Ground truth available via shared cluster-household identifiers,
allowing computation of precision, recall, and F1-score.

---

## 3. Non-Probability Bias Correction

Satellite and administrative data provide full territorial coverage but
are not probability samples. Units in accessible areas (near urban centres,
roads) are systematically over-represented.

### Propensity Score Weighting — Chen et al. (2020)

Estimate the probability of inclusion in the non-probability source:
```
logit(P(S_i = 1)) = alpha + gamma * X_i
```

Reweight non-probability estimates by inverse propensity scores:
```
theta_IPW = sum(w_i * Y_i) / sum(w_i)    where w_i = 1 / ps_i
```

**Implementation:** GPS cluster coordinates (NIGE61FL.shp, 476 clusters)
are used to compute haversine distance to Niamey. A logistic decay function
maps distance to inclusion probability, reflecting the satellite coverage
gradient across Niger's territory.

---

## 4. Extended Fay-Herriot Model with Multisource Auxiliaries

**Standard FH model (Fay & Herriot, 1979):**
```
y_i = x_i'β + v_i + e_i

  v_i ~ N(0, σ²_v)   area random effect
  e_i ~ N(0, ψ_i)    sampling error  (ψ_i = SE_i², known)
```

**EBLUP estimator:**
```
EBLUP_i = γ_i * y_i + (1 - γ_i) * x_i'β_hat

  γ_i = σ²_v / (σ²_v + ψ_i)    shrinkage factor
```

**Extended multisource specification:**
```
y_i = β_0 + β_1*x_census_i + β_2*x_admin_i + β_3*x_satellite_i + v_i + e_i
```

Four model variants are compared (M1 through M4), adding one data source
at a time. Model selection uses AIC, BIC, and the between-area variance
component σ²_v. MSE is estimated via parametric bootstrap (B = 200),
with Prasad-Rao analytical approximation as fallback.

---

## 5. Differential Privacy for Official Statistics

### (ε, δ)-Differential Privacy

A mechanism M satisfies (ε, δ)-DP if for all adjacent databases D, D':
```
P(M(D) ∈ S) ≤ exp(ε) * P(M(D') ∈ S) + δ
```

### Gaussian Mechanism

For a statistic f with L2-sensitivity Δf:
```
M(D) = f(D) + N(0, σ²)

  σ = Δf * sqrt(2 * ln(1.25/δ)) / ε
```

### Calibration for SAE Poverty Rates

For regional poverty rates bounded in [0,1]:

- Global sensitivity: Δf = 1/n_i
- Recommended configuration: ε = 2.0, δ = 1e-5
- Mean absolute error < 0.002 for all Niger regions (n_i ≥ 722)

---

## References

1. Fay, R.E. & Herriot, R.A. (1979). Estimates of income for small places.
   *JASA* 74(366), 269-277.
2. Rao, J.N.K. & Molina, I. (2015). *Small Area Estimation* (2nd ed.). Wiley.
3. Molina, I. (2022). Disaggregating data in household surveys.
   *ECLAC Statistics Series No. 97*.
4. Chen, Y. et al. (2020). Doubly robust inference for non-probability samples.
   *JRSS-B* 82(2), 391-411.
5. Dwork, C. & Roth, A. (2014). *Algorithmic Foundations of Differential Privacy*.
   Foundations & Trends in TCS 9(3-4).
6. Enamorado, T. et al. (2019). Using a probabilistic model to assist merging
   of large-scale datasets. *Political Analysis* 27(2), 222-229.
7. Edochie, I. et al. (2025). Small area estimation of poverty in four West
   African countries. *Journal of Official Statistics* 41(1), 96-124.
