# Provenance and methodological decisions

The starting notes are the author's working notes for Paper A; an unmodified dated copy is archived in `PaperA/archive`. The original remains unchanged.

The inherited Gram definition and projection construction belong to Paper 1 and Paper B (`../GCEEAlgebra`). This repository evaluates the full residual Gram form by quadrature. It does not regenerate the symbolic channel catalogue or use a truncated small-distortion expansion as an exact finite-distortion formula.

Mesh definitions were checked against the scripting directory of the first paper:

- `03_q8_convergence.py`: alternating and seeded random vertex perturbations, manufactured sin(2x)cosh(y) field.
- `05_q8_curved_convergence.py`: internal horizontal midside shifts and curved families.
- `06_q9_convergence.py`: Q9 geometry; its centre node stays at the mean of the four vertices. Q8 and Q9 have identical curved boundaries but need not have identical interior parameterizations.
- `09_radaptivity_example.py`: diagonal tanh layer motivation.

Kernels here are independent implementations, with no import-time dependency on those experimental scripts. Physical-node cubic fitting is a new reference implementation. The author's thesis or thesis code has not been supplied or integrated; it is not claimed as the source of this recovery implementation.

## Corrections to the preliminary framing

1. An exact-derivative oracle is a reference value, not automatically an upper bound.
2. Interpolation error is not local best-approximation error. The broken local best error cannot exceed the conforming FEM error.
3. With identical nodal boundary trace, Poisson Galerkin orthogonality yields EI² = Eh² + |uh-Ihu|². A large interpolation/FEM ratio is not evidence of supercloseness to this interpolant.
4. Random distortion need not eliminate the gap; the present Q8 data provide a counterexample to that general statement.
5. TMOP can use solution-driven targets. Distortion potentials can encode localized requirements. Comparisons must name the actual objectives and target fields.
6. The full physical residual energy sees affine anisotropy; a frozen cubic sees Hessian variation through third derivatives, not arbitrary within-cell variation.
7. Minimization does not automatically imply equidistribution on a fixed anisotropic mesh.

## Conventions

Reference square [-1,1]²; corner order SW, SE, NE, NW; midsides S, E, N, W; Q9 centre last. Geometry is isoparametric. The Jacobian has physical rows and reference columns. The seven Gram channels are x², xy, y², x³, x²y, xy², y³ in physical coordinates about F(0), with Taylor factorials included in the coefficients.

Straight-sided geometry is bilinear, including the initial higher-order meshes with straight sides. Corner determinant minima are exact. For curved Q8/Q9 the determinant lies in tensor degree (3,3); positive Bernstein coefficients give a sufficient whole-cell certificate. A negative coefficient alone would not prove an invalid cell; no such case is silently accepted.

Relocation is a block iteration: recovered physical polynomials and targets stay fixed inside each optimization, then the PDE is reassembled, solved and recovered again. The original connectivity and boundary nodes stay fixed. True errors are diagnostics, never optimization or acceptance data.
