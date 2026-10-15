# Cure-Model Approaches Using Piecewise Exponential Distributions in STEMI Patients Undergoing Primary PCI

This repository contains the R code implementation for parametric cure models (Bernoulli, Poisson, and Negative Binomial frameworks) using Weibull and Piecewise Exponential (PE / Power PE) baseline hazard distributions.

These functions were developed and applied in the study:
**"Cure-Model Approaches Using Piecewise Exponential Distributions in STEMI Patients Undergoing Primary PCI: Prognostic Implications of Apelin-12 and Renal Function"**

## Authors
- **Yolanda Gómez** (Department of Statistics, University of Bío-Bío, Concepción, Chile)
- **Allan Burboa** (Department of Statistics, University of Bío-Bío, Concepción, Chile)

---

## Overview

Classical survival models assume that all individuals will eventually experience the event of interest. Cure models extend this framework by allowing a proportion of long-term event-free (non-susceptible) individuals. 

This repository provides R functions to fit:
- **Mixture Cure Models**: Bernoulli formulation (`BerPE`, `BerPPE`, `BerW`).
- **Non-Mixture (Promotion Time) Cure Models**: Poisson (`PoPE`, `PoPPE`, `PoW`) and Negative Binomial (`NBPE`, `NBPPE`, `NBW`).
- **Baseline Distributions**: Weibull, Piecewise Exponential (PE), and Power Piecewise Exponential (PPE).

Key statistical features implemented in the code:
1. Optimization via Nelder-Mead with BFGS refinement on log-scale parameters to guarantee convergence at the global minimum.
2. Delta-method standard errors evaluated at the maximum likelihood estimates.
3. Computation of individual cure probabilities and average population cure fractions with 95% confidence intervals (via asymptotic normal simulation).
4. Residual diagnostic tools: Cox-Snell residuals and Randomized Quantile residuals.

---

## Prerequisites & Dependencies

To run the R code, ensure you have R installed along with the following required packages:

```R
install.packages(c("msm", "pracma", "MASS"))
