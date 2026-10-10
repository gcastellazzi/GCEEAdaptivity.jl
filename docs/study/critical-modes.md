# Geometric critical modes and same-field comparisons

This extends the displacement-only controls by identifying which fields attain the normalized Gram maximum. It adds no PDE solves and uses no recovered data.

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/run_critical_modes.py
.venv/bin/python scripts/export_critical_modes.py --paper /path/to/PaperA
.venv/bin/python -m pytest -q
```

## Extraction and direct check

For Q4/Q8/Q9 and delta=0/.4/.8, the same interior cell (zero-based cell 5 on the n=4 mesh) is used. It attains the mesh maximum in these cases. Solve G a = lambda H a and normalize a.T H a=1. The corresponding exact Navier polynomial is defined modulo affine fields; subtract the mean displacement and the mean gradient about the physical centroid for plotting. This removes arbitrary affine content without changing the interpolation residual or the quotient norm.

The nine extracted maxima are verified by integrating the physical polynomial and its FE interpolant directly at quadrature order 21, independently of the order-13 operator contraction. Maximum relative difference: 2.496e-9. Maximum normalization defect: 2.443e-15. All relative gaps between the largest and second eigenvalues exceed 0.369, so each displayed largest mode is unique up to sign in these cases. No uniqueness claim is made for other geometries or repeated eigenvalues.

## A fixed field on a fixed cell

At delta=.8, each row fixes the source space's worst physical field and applies all three interpolation operators. Normalization is the same within the row.

| Source of worst field | Q4 response | Q8 response | Q9 response |
| --- | ---: | ---: | ---: |
| Q4 | 3.273 | 1.900 | 0.802 |
| Q8 | 2.582 | 2.444 | 0.925 |
| Q9 | 1.540 | 1.941 | 2.003 |

The diagonal attains each space's maximum. Every column is bounded by its diagonal maximum in all 27 cross-space evaluations (including delta=0/.4). A lower maximum does not imply a smaller interpolation error on every field. Nesting of approximation spaces orders their best-approximation errors, not necessarily the errors of their nodal interpolation operators. This comparison is local interpolation, not Galerkin error or locking.

## Covariance

Decompose the polynomial by total degree about the specified physical centre. For displacement interpolation residuals R2 and R3, the full squared error is ||R2||² + ||R3||² + 2(R2,R3). In the Q8 worst field at delta=.8, the cross term contributes 27.664% of the squared total. Dropping only this cross-degree term gives 85.05% of the full error norm. This is not the same ablation as dropping every off-diagonal channel entry. The degree decomposition depends on the expansion centre, while the full spectral response is unchanged by the tested basis and rigid-coordinate changes.

## Files and checks

- `gceeadaptivity/critical_modes.py`: eigenmode extraction, canonical affine representative and direct response.
- `results/critical_modes/modes.json`: nine mode diagnostics.
- `results/critical_modes/mode_*.npz`: polynomial coefficients, centres, scales, affine corrections and eigenvectors.
- `results/critical_modes/cross_space.json`: 27 fixed-field evaluations.
- `results/critical_modes/validation.json`: tolerances and maximum defects.
- `tests/test_critical_modes.py`: direct spectral attainment, normalization, finite-difference displacement gradients, homogeneous Navier equilibrium and the degree-covariance identity.
- `results/pytest_critical_modes.log`: all 30 tests passed, including the preceding scalar and elastic controls.

The manuscript adds one six-panel figure and one table. Exact field arrows and interpolation-residual arrows have separate magnifications; magnification is common across spaces within each row. Error-density colours share one scale. This geometric characterization still requires a fixed material and polynomial class; it does not supply a universal admissible distortion threshold.
