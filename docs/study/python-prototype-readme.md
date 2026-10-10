# GCEEAdaptivity

Research code accompanying **Recovered Gram indicators for mesh relocation** (Paper A, 3 October 2026). The five initial experiments in `notes_A_recovery_radaptivity.tex` are implemented and their first campaigns are included. The initial relocation is scalar 2D; a separate displacement-only plane-strain control now separates geometry, equilibrium and recovery. This is a research prototype, not a certified adaptive solver.

## Reproduce

From this directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
export OPENBLAS_NUM_THREADS=1
python -m pytest -q
python scripts/run_recovery.py
python scripts/run_relocation.py
python scripts/check_benchmarks.py
python scripts/export_paper.py --paper '/Users/gcastellazzi/Dropbox/Papers/2026_a-priori_paper#2/PaperA'
```

Alternatively install `requirements-lock.txt` in a fresh environment to reproduce the recorded versions (Python 3.14.2). The general package supports Python >=3.10; the pinned dependencies may require a newer Python. The full two campaigns take approximately one minute on the development machine; timings are hardware-dependent. The `--quick` recovery option uses only levels 4 and 8 and overwrites the standard recovery outputs; restore the full campaign before paper export or benchmark checks, which require all 93 cases.

Compile the manuscript from `PaperA` with `latexmk -pdf -interaction=nonstopmode -halt-on-error paperA.tex`.

## Implemented experiments

| Step | Experiment | Implementation / data |
| --- | --- | --- |
| 1 | Recovery error separated from FE data error | Cubic weighted physical least squares; exact-nodal and computed-nodal fits on the same patches |
| 2 | Recovered versus oracle full Gram indicator | Seven quadratic/cubic channels, full covariance, diagnostic Gram error budgets |
| 3 | Both effectivities and estimator controls | Interpolation and Galerkin errors; ZZ-type gradient recovery and uncalibrated residual indicator |
| 4 | Frozen recovered-field relocation | Differentiated physical polynomial interpolation energy, fixed boundary and topology, feasible Armijo steps |
| 5 | Re-solved relocation and target controls | Four outer cycles; GCEE, uniform TMOP shape target, recovered-Hessian TMOP shape target |

The recovery campaign has **93 Poisson solves** on Q4/Q8/Q9 meshes, levels 4/8/16, regular, checkerboard, random (two seeds) and curved geometries. The relocation campaign has **12 sequences, 60 PDE solves and 48 frozen optimization runs**, for Q4 at levels 8/12 on a diagonal layer and an anisotropic Gaussian bump. All three objectives start from exactly the same mesh for each field and level.

## Displacement-only elastic control

The additional [control campaign](DISPLACEMENT_CONTROLS.md) contains 150 solves on 30 straight-sided Q4/Q8/Q9 meshes, with affine and exact non-affine equilibrated boundary probes, zero body forces and no prescribed tractions. A normalized eight-mode Gram spectrum characterizes map response without a selected load or recovered field. Run `scripts/run_elastic_patch.py` and `scripts/export_elastic_patch.py --paper /path/to/PaperA` to reproduce it. The full test suite now has 30 passing tests.

The [critical-mode extension](CRITICAL_MODES.md) extracts nine geometric maximizers and checks 27 same-field cross-space responses, without additional PDE solves.

## Layout

- `gceeadaptivity/fem.py`: shapes, mappings, manufactured fields, sparse Poisson solve, independent errors and validity certificates.
- `gceeadaptivity/recovery.py`: scaled fitting, residual Gram matrices and recovery diagnostics.
- `gceeadaptivity/estimators.py`: recovered-gradient and volume/flux residual controls.
- `gceeadaptivity/relocation.py`: automatic differentiation, objectives, safeguards and re-solve loop.
- `scripts/`: campaign runners, independent numerical checks and data-driven LaTeX/figure export.
- `results/recovery.csv`, `recovery.json`: global records; `local_*.npz`: local squared errors, indicators and budgets.
- `results/relocation.json`, `frozen.json`: true errors after each solve and complete inner histories.
- `results/meshes_*.npz`: initial plus four relocated meshes per sequence.
- `VERIFICATION.md`, `PROVENANCE.md`, `ROADMAP.md`, `literature/README.md`: checks, definitions, open work and bibliography.

## First findings and limits

At n=16 on the Q8 checkerboard mesh, the recovered indicator is 0.924 times interpolation error but 2.101 times Galerkin error. Recovery accuracy and the interpolation--Galerkin gap are separate issues. The Gram recovery budget is diagnostic: exact derivatives enter it, so it is not yet an operational uncertainty estimate.

GCEE reduces the final Galerkin error in all four relocation cases (about 9--25%). A single accepted outer update raises that error despite reducing the frozen objective. No true-error acceptance test conceals this result. Uniform-target smoothing outperforms GCEE for the coarse bump. The Hessian control adapts shape but not target area; its results do not establish superiority over general solution-driven TMOP. We implement the documented metric independently and do not execute MFEM's optimizer.

Current relocation is Q4 and scalar Poisson. Curved Q8/Q9 are covered by recovery experiments, not by relocation. Elastic relocation, 3D, internal-mode projections, thesis recovery integration, computable recovery uncertainty and reliable post-solve acceptance remain future work. The diagonal Gram approximation is never used.
