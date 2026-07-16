using SciMLTesting, FiniteStateProjection, Test
using JET

const CATALYST_REEXPORTS = Tuple(intersect(names(FiniteStateProjection), names(FiniteStateProjection.Catalyst)))
const API_DOCS_IGNORE = Tuple(
    unique(
        (
            :SymbolicUtils,
            :NaiveIndexHandler,
            CATALYST_REEXPORTS...,
        )
    )
)

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
        # These names are still non-public in their resolvable owners, or are
        # re-exported from a non-owner. Drop each ignore as its owner ships a public
        # declaration that FiniteStateProjection can actually resolve.
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
                :NullParameters,   # SciMLBase (public in 3.30+, but Catalyst 15 pins 2.153.1 where it is not)
            ),
        ),
    ),
    api_docs_kwargs = (;
        rendered = true,
        # `@reexport using Catalyst` intentionally exposes Catalyst's symbolic
        # stack. Require docstrings for FiniteStateProjection-owned public API
        # here, and leave dependency API docs to their owner packages.
        ignore = API_DOCS_IGNORE,
        rendered_ignore = API_DOCS_IGNORE,
    ),
    # Heavy `@reexport using Catalyst` plus the symbolic stack make ~23 names
    # implicit; making them all explicit is a risky mass refactor. Tracked in #60.
    ei_broken = (:no_implicit_imports,),
)
