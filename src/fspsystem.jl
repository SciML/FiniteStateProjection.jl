"""
Thin wrapper around `Catalyst.ReactionSystem` for use with this package.

Constructor: `FSPSystem(rs::ReactionSystem[, ih=DefaultIndexHandler(); combinatoric_ratelaw::Bool=true])`
"""
struct FSPSystem{IHT <: AbstractIndexHandler, RT}
    rs::ReactionSystem
    ih::IHT
    rfs::RT
end

function FSPSystem(rs::ReactionSystem, ih=DefaultIndexHandler(); combinatoric_ratelaw::Bool=true)
    isempty(Catalyst.get_systems(rs)) ||
        error("Supported Catalyst models can not contain subsystems. Please use `rs = Catalyst.flatten(rs::ReactionSystem)` to generate a single system with no subsystems from you Catalyst model.")
    any(eq -> !(eq isa Reaction), equations(rs)) &&
        error("Catalyst models that include constraint ODEs or algebraic equations are not supported.")

    rfs = create_ratefuncs(rs, ih; combinatoric_ratelaw=combinatoric_ratelaw)
    FSPSystem(rs, ih, rfs)
end

"""
    build_ratefuncs(rs, ih; state_sym::Symbol, combinatoric_ratelaw::Bool)::Vector

Return the rate functions converted to Julia expressions in the state variable
`state_sym`. Abundances of the species are computed using `getsubstitutions`.

See also: [`getsubstitutions`](@ref), [`build_rhs`](@ref)
"""
function build_ratefuncs(rs::ReactionSystem, ih::AbstractIndexHandler; state_sym::Symbol, combinatoric_ratelaw::Bool=true)
    substitutions = getsubstitutions(ih, rs, state_sym=state_sym)

    return map(Catalyst.reactions(rs)) do reac
        jrl = jumpratelaw(reac; combinatoric_ratelaw)
        jrl_s = substitute(jrl, substitutions)
        toexpr(jrl_s)
    end
end

function create_ratefuncs(rs::ReactionSystem, ih::AbstractIndexHandler; combinatoric_ratelaw::Bool=true)
    paramsyms = Symbol.(Catalyst.parameters(rs))

    return tuple(map(ex -> compile_ratefunc(ex, paramsyms),
                     build_ratefuncs(rs, ih; state_sym=:idx_in, combinatoric_ratelaw))...)
end

function compile_ratefunc(ex_rf, params)
    # Make this nicer in the future
    ex = :((idx_in, t, $(params...)) -> $(ex_rf)) |> MacroTools.flatten
    @RuntimeGeneratedFunction(ex)
end
