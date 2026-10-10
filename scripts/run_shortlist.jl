# Shortlist at a modest resolution. Published Newton / high-order numbers are
# printed in docs/tmop-benchmarks.md; this script does not try to match them.
#
#   julia --project=. scripts/run_shortlist.jl
#   julia --project=. scripts/run_shortlist.jl 6 4

using GCEEAdaptivity

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 6
iterations = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
runs = run_shortlist(; n=n, iterations=iterations, order=3)
println()
println("Adapted meshes keep a positive Jacobian. L2 is reported for the Poisson front and the oblique shock.")
for r in runs
    println(r.name, ": min area ", r.min_area)
end
