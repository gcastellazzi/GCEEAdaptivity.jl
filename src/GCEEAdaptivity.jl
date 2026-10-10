"""
    GCEEAdaptivity

TMOP benchmark harness on quadrilateral meshes from FinEtools.

The six cases are the shortlist in `docs/tmop-benchmarks.md`: the circular-front
Poisson problem, the analytic ring, the IMR sinusoidal band, the oriented
oblique shock, and the MFEM size and material indicators. A caller-supplied
positive weight can replace `‖∇u‖` in the equidistributed size target; that is
the slot an a-priori estimator occupies. [`gcee_target`](@ref) fills that slot
with the scalar Q4 Gram indicator.
"""
module GCEEAdaptivity

using FinEtools
using ForwardDiff
using LinearAlgebra
using Printf

include("elements.jl")
include("meshes.jl")
include("fields.jl")
include("fesolve.jl")
include("recovery.jl")
include("indicators.jl")
include("relocation.jl")
include("part.jl")
include("elasticity.jl")
include("metrics.jl")
include("quadrature.jl")
include("targets.jl")
include("gcee.jl")
include("functional.jl")
include("optimize.jl")
include("poisson.jl")
include("benchmarks.jl")
include("studies.jl")

export unit_square_mesh, QuadPatch, integrate_scalar, minimum_area, nelems
export mu1, mu2, mu7, mu9, mu14, mu55, mu56, mu58, mu77, mu94, mu107, metric_value
export tmop_functional, element_metrics, field_target, relocate
export sinusoidal_eta, sinusoidal_target, ring_feature, ring_size_target, ring_aspect_target
export size_feature, material_feature
export circular_front, circular_front_grad, oblique_shock, oblique_shock_grad
export equidistributed_isotropic, gradient_shape_size, orientation_target
export q4_gram_indicator, gcee_indicators, gcee_target
export size_indicator, discrete_size_target, material_indicator, material_indicator_grad
export material_interface_target
export solve_poisson, poisson_errors
export BenchmarkRun, run_sinusoidal_band, run_analytic_ring, run_poisson_front
export run_oblique_shock, run_discrete_size, run_material_interface
export run_shortlist, print_runs

# Recovered-indicator study: elements, meshes, fields, recovery, relocation.
export CORN, MIDS, EXP, HIGH, monomials, basis, shapes, reference, gauss_legendre
export Mesh, structured_mesh, curved_beam, mesh_geometry, corner_determinants, certify
export with_points, ncells, npoints, node_layout, jitter, annulus_map
export manufactured, field_gradient, field_hessian, field_force, field_coefficients
export coefficients_at, solve_field, field_errors
export Recovery, fit, physical_coefficients, element_centers, element_recovery
export gram, contraction, diagnostics
export zz_indicator, residual_indicator
export target_inverse, relocation_objective, relocate_nodes, recovered_at
export recovery_cases, run_recovery_study, relocation_experiment
export run_relocation_study, run_part_study
export material, strain, elastic_probe, navier_residual, probe_values, probe_coefficients
export PatchSystem, patch_solve, assess_patch, navier_basis, geometric_response
export strain_energy, discrete_strain
export CriticalMode, critical_mode, direct_response, mode_value, mode_gradient, strain_of

end
