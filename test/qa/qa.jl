using SciMLTesting, FiniteStateProjection, Test
using JET

run_qa(
    FiniteStateProjection;
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
    # Pre-existing Aqua findings tracked in SciML/FiniteStateProjection.jl#60.
    aqua_broken = (
        :ambiguities,        # pairedindices DefaultIndexHandler{0} overlap (indexhandlers.jl 89/97/116)
        :unbound_args,       # pairedindices unbound type params (indexhandlers.jl 97/116)
        :undefined_exports,  # @reexport using Catalyst re-exports names absent from loaded deps
    ),
    # JET reports 3 issues (NaiveIndexHandler @deprecate kwcall; NullParameters used
    # but not imported in build_rhs.jl/build_rhs_ss.jl). Tracked in #60.
    jet_broken = true,
    ei_kwargs = (;
        # Names re-exported by a non-owner dependency (resolve to the owner as base
        # libraries adopt public/owner declarations).
        all_qualified_accesses_via_owners = (;
            ignore = (
                :get_systems,  # owner ModelingToolkit, accessed via Catalyst
                :scalarize,    # owner Symbolics, accessed via ModelingToolkit
                :value,        # owner Symbolics, accessed via ModelingToolkit
            ),
        ),
        # `Catalyst = "15"` pins the SciMLBase 2.x / ModelingToolkit 9.x / Symbolics 6.x
        # ecosystem, so the public-API releases (SciMLBase 3.27 / MTK 11 / Symbolics 7)
        # are NOT resolvable here. Verified on Julia 1.12 against the resolved stack
        # (Catalyst 15.0.11, ModelingToolkit 9.84.0, Symbolics 6.58.0, SciMLBase 2.153.1):
        # each of these names is still non-public in its owner at the resolvable version,
        # and MacroTools/Catalyst have not declared theirs public. Drop each as its owner
        # ships a public declaration that FiniteStateProjection can actually resolve.
        all_qualified_accesses_are_public = (;
            ignore = (
                :_symbol_to_var,   # Catalyst (non-public)
                :get_systems,      # Catalyst (owner ModelingToolkit; still non-public)
                :alias_gensyms,    # MacroTools (non-public)
                :flatten,          # MacroTools (non-public)
                :prewalk,          # MacroTools (non-public)
                :resyntax,         # MacroTools (non-public)
                :striplines,       # MacroTools (non-public)
                :scalarize,        # ModelingToolkit (owner Symbolics; still non-public)
                :value,            # ModelingToolkit (owner Symbolics; still non-public)
                :varmap_to_vars,   # ModelingToolkit (non-public)
            ),
        ),
    ),
    # Heavy `@reexport using Catalyst` plus the symbolic stack make ~23 names
    # implicit; making them all explicit is a risky mass refactor. Tracked in #60.
    ei_broken = (:no_implicit_imports,),
)
