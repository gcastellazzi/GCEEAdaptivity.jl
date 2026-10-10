# Displacement-only elasticity: separate the three questions

## Definition and reproduction

This control adds **150 plane-strain solves on 30 straight-sided meshes (10 geometric configurations, each with three spaces)**. It does not perform elastic mesh relocation. E=1, nu=0.3; no body force is assembled; all boundary displacements are prescribed nodally; interior degrees of freedom solve equilibrium. Boundary reactions are recorded. The same physical boundary field is used for every distortion and element space.

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/run_elastic_patch.py
.venv/bin/python scripts/export_elastic_patch.py --paper /path/to/PaperA
.venv/bin/python -m pytest -q
```

The campaign uses Q4/Q8/Q9, n=4/8, alternating distortion delta=0/.2/.4/.6/.8. Midside nodes remain at straight-side midpoints; centre nodes remain at vertex means. The physical cell maps are identical across spaces. The original scalar and relocation campaigns are unchanged.

## Three distinct measurements

1. **Map response before selecting a field.** Form the full elastic interpolation-residual Gram matrix on the eight-dimensional space of homogeneous Navier quadratic/cubic fields modulo affine fields. Normalize by the exact strain energy after subtracting mean strain. The largest generalized eigenvalue gives q_K squared. This construction uses no recovery, boundary loading or FEM solve. Its meaning is conditional on the fixed material, interpolation operator and polynomial class. It is not an error estimator for an unspecified solution or a universal quality threshold.
2. **Equilibrium with prescribed displacements.** The affine patch verifies constant-strain reproduction. Four non-affine fields are gradients of harmonic cubic/quartic potentials; their divergence and vector Laplacian vanish identically, so they require no manufactured body forces. Only boundary nodal values are imposed. These probes test nonuniform strain, not volumetric locking. Their boundary interpolation error is retained and can differ between FE spaces.
3. **Recovery from computed data.** Exact cubic fitting to exact nodal data is exact, and cubic Taylor truncation is zero. Replacing those data by the FEM solution isolates data-induced recovery loss. A separate recovered-polynomial equilibrium defect is stored, since componentwise least squares need not satisfy Navier even when body forces are zero.

## Normalization

For a basis p_i of the Navier polynomial quotient, G_ij is the energy inner product of p_i-I_K p_i and p_j-I_K p_j. H_ij is the inner product of their exact strains after subtracting mean strain. H is positive definite on the quotient. Both forms are independent of the affine representative. q_K=sqrt(lambda_max(G,H)); the mesh statistic is max_K q_K. Full local G, H and spectra are stored, alongside the RMS of the modal eigenvalues.

The normalized spectrum is invariant under polynomial basis changes, translation, uniform scaling and rigid rotation for this isotropic material. These properties are tested. A raw coefficient-space eigenvalue without H would not have this basis-independent meaning. q_K can exceed one because nodal interpolation is not an energy projection. It is not a stiffness condition number and not a bound on every nonpolynomial field.

## First results

For the fixed cubic cosine probe at n=8:

| Space | q_max at delta=0 | q_max at delta=.8 | FEM error growth factor |
| --- | ---: | ---: | ---: |
| Q4 | 2.121 | 3.273 | 1.547 |
| Q8 | 0.849 | 2.444 | 4.334 |
| Q9 | 0.849 | 2.003 | 1.509 |

At delta=.8, Q8 has recovered/interpolation effectivity 0.866 and recovered/FEM effectivity 1.408. Removing applied forces does not remove the interpolation--Galerkin gap or the recovery error. A worst-case modal statistic and a particular imposed field need not deteriorate by the same factor.

All affine patches have energy error below 5.7e-14. Q9 reproduces the quadratic probes on these bilinear maps. The oracle/interpolation identity has maximum absolute defect 3.2e-16; the exact-nodal coefficient recovery error is below 1.8e-13. The largest recorded higher-order quadrature change is 3.3e-5, using a 1e-10 denominator floor for reproduced fields. The geometric maximum changes by less than 2.5e-9. There are 27 passing tests including the original suite. Full records and thresholds are in `results/elastic_patch/validation.json`.

The classical affine patch and homogeneous elasticity operator are described respectively in [Abaqus verification](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEVERRefMap/simaver-c-basicelemover.htm) and [MFEM examples](https://mfem.org/examples/). The normalized spectral response here is an explicit construction, not a reliability theorem attributed to those sources.

## Files and next questions

`gceeadaptivity/elastic_patch.py` implements the stiffness solve, analytic probes, vector recovery diagnostics and spectral response. `tests/test_elastic_patch.py` checks affine patches, exact-polynomial identities, geometric invariance and Rayleigh bounds. All local fields and spectra are under `results/elastic_patch`. The manuscript contains two additional vector figures and one additional table.

Next: use a selected modal-response tolerance to define a conditional distortion limit; compare complete spectra and worst modes, not only maxima; extend to 3D and transition spaces. Equilibrium-constrained recovery can then be compared against the present unconstrained fit. Displacement control still excites the structure through reactions; it isolates the loading mechanism without eliminating equilibrium.

## Critical-mode extension

The geometric maximizers and 27 cross-space evaluations are now available; see [CRITICAL_MODES.md](CRITICAL_MODES.md). The full suite now has 30 passing tests.
