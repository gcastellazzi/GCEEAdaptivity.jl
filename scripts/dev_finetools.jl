using Pkg

const ROOT = get(ENV, "FINETOOLS_DEV_ROOT", "/Users/gcastellazzi/Dropbox/GitHub")

Pkg.activate(joinpath(@__DIR__, ".."))

for pkg in ("FinEtools.jl", "FinEtoolsDeforLinear.jl")
    path = joinpath(ROOT, pkg)
    if isdir(path)
        @info "Developing local package" pkg path
        Pkg.develop(path=path)
    else
        @warn "Local package checkout not found" pkg path
    end
end

Pkg.add("ForwardDiff")
Pkg.instantiate()
Pkg.precompile()
