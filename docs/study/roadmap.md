# Ordered research programme

The first five experiments follow the order in the source notes. The initial implementation and campaigns are complete; the research questions below remain open.

1. **Recovery alone — initial campaign complete.** Full-rank cubic preservation, exact-nodal/computed-nodal separation, quadratic/cubic coefficient errors and Gram-weighted budgets. Next: integrate the author's thesis recovery when its source is available; vary patch radius and polynomial order independently; test boundaries and more severe distortions. Develop a computable recovery uncertainty estimate and compare it with the exact diagnostic budget.
2. **Oracle survival — initial campaign complete.** Full Gram form, both effectivities, random/checkerboard/curved families. Next: isolate Taylor truncation from recovery error, investigate inactive and nearly cancelling channels, and test an H1 projection for internal-mode elements without conflating it with interpolation.
3. **Estimator comparison — initial campaign complete.** ZZ-type recovered gradient and volume/flux residual control. Next: order-aware residual scaling and boundary-data oscillation; compare an equilibrated flux estimator as a computable acceptance diagnostic. The current uncalibrated residual indicator is not a certified bound.
4. **Frozen relocation — initial campaign complete.** AD gradient checks, corner-validity constraint, affine anisotropy and cubic variation retained. Next: quantify surrogate mismatch against step size and recovery uncertainty; compare full and diagonal Gram objectives only as a controlled ablation. Add curved relocation only with a whole-cell validity certificate.
5. **Repeated relocation — initial campaign complete.** Four re-solves and two explicit TMOP shape controls. Next: design an implementable acceptance test after the re-solve, without analytical error. Add target-size adaptation to the solution-driven baseline and test more seeds/levels before any general performance claim.

## Subsequent extensions

- A constrained equidistribution analysis: the isotropic density argument in the paper is a model, not a theorem for fixed quadrilateral connectivity.
- Q8/Q9 relocation and boundary-conforming geometry updates.
- Displacement-only plane-strain diagnostics are now implemented (150 solves; see `DISPLACEMENT_CONTROLS.md`). Elastic relocation and 3D transition elements remain to be developed using Paper B's obstruction and covariance analysis. Keep displacement recovery, rigid modes and geometric validity separate.
- A wall/opening transition benchmark with identical initial geometry for every method; state distortion limits conditional on the field, space, error tolerance and validity criterion.
- Performance engineering only after the diagnostics and acceptance logic are fixed.

The first numerical results do not justify calling the method a reliable FEM-error estimator or a monotonically improving adaptive algorithm.

## Geometric critical modes — completed initial extension

The worst normalized Navier-class fields are extracted on three distortions and directly verified. Same-field comparisons across Q4/Q8/Q9 and quadratic–cubic covariance are documented in `CRITICAL_MODES.md`. Remaining steps are a tolerance-dependent distortion criterion, additional map families, and 3D/transition-space extension.
