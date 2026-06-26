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
        # Non-public names qualified from other packages. The SciMLBase/ModelingToolkit/
        # Symbolics public-API releases are held back here by `Catalyst = "15"` (which
        # pins the SciMLBase 2.x / ModelingToolkit 9.x / Symbolics 6.x ecosystem), so the
        # symbolic-stack names are still non-public in the resolved versions; MacroTools
        # and Catalyst have not made these public. Drop each as its owner releases a
        # public declaration that FiniteStateProjection can actually resolve.
        all_qualified_accesses_are_public = (;
            ignore = (
                :_symbol_to_var,   # Catalyst
                :get_systems,      # Catalyst
                :alias_gensyms,    # MacroTools
                :flatten,          # MacroTools
                :prewalk,          # MacroTools
                :resyntax,         # MacroTools
                :striplines,       # MacroTools
                :scalarize,        # ModelingToolkit (owner Symbolics)
                :value,            # ModelingToolkit (owner Symbolics)
                :varmap_to_vars,   # ModelingToolkit
            ),
        ),
    ),
    # Heavy `@reexport using Catalyst` plus the symbolic stack make ~23 names
    # implicit; making them all explicit is a risky mass refactor. Tracked in #60.
    ei_broken = (:no_implicit_imports,),
)
