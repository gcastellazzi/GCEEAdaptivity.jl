# GCEEAdaptivity

Two studies in one package, both behind Paper A, *Recovered Gram indicators for
mesh relocation*.

**The recovered-indicator study** (Experiments 1--6) measures how a
geometry-induced interpolation indicator survives being rebuilt from a computed
solution, and what happens when it is used to move nodes. It carries its own
Q4/Q8/Q9 elements, mesh families, manufactured fields, patch recovery, Gram
budget, ZZ and residual estimators, a displacement-only plane-strain control,
the critical-mode extraction, and the quarter-annulus curved beam. Field
derivatives and relocation gradients come from forward-mode AD; mesh
perturbations use a stateless hash, not a random number generator, so a mesh
depends only on its arguments.

**The TMOP shortlist** compares that indicator against published target
constructions, with the quality metric and the solver held fixed.

Everything the manuscript needs:

```
julia --project=. test/runtests.jl
julia --project=. scripts/export_paper.jl --paper /path/to/PaperA
julia --project=. scripts/plot_benchmarks.jl
```

The first Python prototype of the recovered-indicator study has been retired;
its methodological notes are kept under [`docs/study/`](docs/study/).

The TMOP side uses the local
[FinEtools](https://github.com/PetrKryslUCSD/FinEtools.jl) checkout for the grid
and the element basis.

The note in [`docs/tmop-benchmarks.md`](docs/tmop-benchmarks.md) is the source
of the six cases. This package builds those targets, integrates the TMOP
objective, moves the nodes, and, where the case has an exact solution, solves
Poisson with bilinear elements before and after the relocation.

The node step is a scaled gradient with an Armijo line search. Corners stay
fixed. Nodes on the edges of the unit square slide along the edge. MFEM's
Newton solver, high-order geometry and published objective values are not
reproduced. The numbers in the note stay the reference for those runs.

## Setup

From this directory. `dev_finetools.jl` looks for the FinEtools checkouts
beside this repository; set `FINETOOLS_DEV_ROOT` if they live elsewhere:

```bash
julia --project=. scripts/dev_finetools.jl
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. scripts/run_shortlist.jl
```

`run_shortlist.jl` takes an optional element count per side and an iteration
count: `julia --project=. scripts/run_shortlist.jl 8 10`.

Mesh figures and the uniform-versus-TMOP convergence of the circular front:

```bash
julia --project=. scripts/plot_benchmarks.jl
```

The script writes the Poisson and shock solution figures, the TMOP / GCEE
meshes, the same meshes colored by the metric on the final nodes
(`*_metric.pdf`: metric id, final functional and maximum in the title), the
convergence plots, and the matching CSV files under `figures/` and `data/`.

## Cases

| Runner | Benchmark |
| --- | --- |
| `run_poisson_front(sharpness=20)` | Kolev et al. 2026, §8.1, circular front |
| `run_poisson_front(sharpness=200)` | Dobrev et al. 2022, §4.1, r-branch, metric μ₉ |
| `run_analytic_ring(kind=:paper)` | hr paper §4.2, radii 0.2 and 0.3 |
| `run_analytic_ring(kind=:tid9)` or `:tid10` | MFEM `mesh-optimizer -tid 9` and `-tid 10` |
| `run_sinusoidal_band(variant=:size\|:shape\|:both)` | Dobrev et al., IMR 2019, §6.1 |
| `run_oblique_shock()` | Kolev et al. 2026, §8.2, orientation target, metric 107 |
| `run_discrete_size()` | MFEM `-tid 5`, `-mid 94` |
| `run_material_interface()` | MFEM `material_indicator_2d`, metric μ₇ |

Every runner accepts `driver`. The published target is `:gradient` on the
Poisson front and `:analytic` on the other cases. `driver = :gcee` keeps that
case's metric and replaces the size and aspect field with the scalar Q4 Gram
indicator of the same feature (the exact solution, or the unclamped indicator
when the benchmark has no PDE). The Gram target is frozen on the initial mesh.
`run_shortlist` runs both drivers. The plot script writes the Poisson
convergence at both sharpness values, the shock convergence, and the TMOP / GCEE
meshes for the ring, the three band variants, the interface and the discrete
size, plus a copy of each mesh colored by the metric measured on the final nodes.

Uniform squares already give μ₂ = 0 for `W = I`. A pure shape target does not
move the Cartesian meshes used here. The 2026 paper reduces the `L²` error by
putting a PDE objective next to the metric. That objective is not part of
this relocation.
