# Matrix Function Approximation with Chebyshev Polynomials

## Overview

This project implements Chebyshev polynomial methods for efficiently approximating matrix functions of the form

[
f(A)v
]

where (A) is a large sparse matrix and (v) is a vector.

Rather than computing the full matrix function (f(A)) through eigendecomposition, the implementation uses Chebyshev polynomial expansions and matrix-vector multiplications to obtain accurate approximations at a significantly lower computational cost.

The project was completed as part of Numerical Analysis coursework at UC Berkeley.

---

## Motivation

Computing matrix functions directly can be computationally expensive for large matrices.

A standard eigendecomposition approach requires approximately (O(N^3)) operations, making it impractical for large-scale problems.

This project uses:

* Chebyshev polynomial approximation
* Gauss-Chebyshev quadrature
* Three-term recurrence relations
* Sparse matrix-vector products

to efficiently approximate (f(A)v) without explicitly forming (f(A)).

---

## Files

| File            | Description                                      |
| --------------- | ------------------------------------------------ |
| `main.m`        | Runs all numerical experiments                   |
| `cheb_eval.m`   | Evaluates Chebyshev polynomials using recurrence |
| `cheb_coeff.m`  | Computes Chebyshev coefficients using quadrature |
| `cheb_matvec.m` | Approximates (f(A)v) using a Chebyshev expansion |

---

## Numerical Experiments

### Relative Error Analysis

The approximation was tested for:

* (N = 100, 200, \ldots, 1000)
* (C = 0.01, 0.1, 0.25, 0.5)

where the polynomial degree is chosen as

[
n = CN.
]

Relative errors were computed against the exact solution obtained using MATLAB's `expm` function.

Results showed that larger values of (C) produced significantly improved accuracy, with (C = 0.5) consistently achieving relative errors below (10^{-8}).

---

### Runtime Analysis

Runtime experiments were conducted for

[
N = 500, 1000, \ldots, 3000.
]

A quadratic curve was fit to the observed computational times.

The measured runtimes closely followed quadratic growth, confirming the expected scaling behavior of the implementation.

---

## Skills Demonstrated

* Numerical Linear Algebra
* Scientific Computing
* MATLAB Programming
* Sparse Matrix Computations
* Algorithm Design
* Computational Performance Analysis
* Approximation Theory
