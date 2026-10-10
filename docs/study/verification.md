# Verification record — 3 October 2026

The full campaigns and their independent checks completed successfully. Exact environments and source hashes are in `results/validation.json`; console logs and raw results are retained.

## Checks

- 21 unit tests passed (`results/pytest.log`). They cover affine patches, nodal reproduction, cubic recovery, Q9 quadratic reproduction, curved mapping certificates, representable-field residuals, Gram identities and budgets, AD directional gradients, validity-preserving descent, affine anisotropy and nonzero cubic response at a zero-Hessian layer centre.
- All 93 recovery meshes have a positive whole-cell determinant certificate.
- All 93 reverse-triangle recovery budgets and their fitting/data decompositions pass.
- All 48 frozen optimization runs decrease their specified objective and maintain the corner determinant floor. These checks do not assert monotonic decrease of the true FEM error.
- Maximum direct Pythagoras relative defect: 3.25608065e-06 (tolerance 1e-5).
- Maximum relative change in six initial-mesh quadrature comparisons: 0.000151409596 (tolerance 2e-4), in the curved Q9 residual control.
- Maximum relative change in either true error on all 12 final relocated meshes: 1.27019286e-08 (tolerance 2e-8).

Recovery quadrature checks compare orders 7/11; initial relocation checks compare 9/13. Final meshes are independently reassembled at order 13 and errors integrated at order 15, compared with stored orders 9/11. These selected quadrature checks do not bound every quantity on every possible distorted mesh.

## Final GCEE error ratios

| Field | n | Final / initial Galerkin error |
| --- | ---: | ---: |
| layer | 8 | 0.843928 |
| bump | 8 | 0.906179 |
| layer | 12 | 0.746875 |
| bump | 12 | 0.751589 |

The n=8 layer improves to about 0.838 after three cycles, then rises to about 0.844 after the fourth. The uniform target outperforms GCEE for the n=8 bump. Both observations are retained in the manuscript.

No continuum reliability theorem, global optimizer convergence, unknown-solution recovery bound, 3D result or general TMOP superiority is established by these tests.

## Displacement-only elastic extension

The additional campaign has 150 solves on 10 geometric configurations, each represented with Q4/Q8/Q9. Its validation is in `results/elastic_patch/validation.json`; the complete 27-test suite passed (`results/pytest_elastic_patch.log`). New tests cover affine elastic patches, exact equilibrated quadratic/cubic fields, oracle identities, recovery budgets, Q9 quadratic reproduction, generalized Rayleigh bounds and translation/scale/rotation/basis invariance of the normalized geometric spectrum. The original records above refer to the initial scalar campaigns.

## Critical-mode extension

Nine largest modes and 27 cross-space responses pass direct higher-order integration checks. The maximum relative spectral/direct difference is 2.496e-9; all cross-space responses lie below the corresponding target-space maximum. The complete 30-test suite passed (`results/pytest_critical_modes.log`). New tests independently check displacement gradients by finite differences, Navier equilibrium, affine normalization and the degree-covariance decomposition. See `results/critical_modes/validation.json`.
