var documenterSearchIndex = {"docs":
[{"location":"examples.html#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples.html#Birth-Death-Model","page":"Examples","title":"Birth-Death Model","text":"","category":"section"},{"location":"examples.html","page":"Examples","title":"Examples","text":"This example models a linear birth-death process. The reaction network is easily defined using Catalyst.jl. Our truncated state space has length 50, which is enough for this simple system.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"This system has no conserved quantities, so we use a NaiveIndexHandler to map from a one-dimensional array with offset 1 to the state of the system. See Index Handlers for more details.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"using FiniteStateProjection, DifferentialEquations\n\n@parameters r1, r2\nrs = @reaction_network begin\n    r1, 0 --> A\n    r2, A --> 0\nend r1 r2\n\nsys = FSPSystem(rs)\n\n# Parameters for our system\nps = [ 10.0, 1.0 ]\n\n# Initial values\nu0 = zeros(50)\nu0[1] = 1.0\n\nprob = convert(ODEProblem, NaiveIndexHandler(sys, 1), sys, u0, 10.0, ps)\nsol = solve(prob, Vern7(), atol=1e-6)","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"(Image: Visualisation)","category":"page"},{"location":"examples.html#Telegraph-Model","page":"Examples","title":"Telegraph Model","text":"","category":"section"},{"location":"examples.html","page":"Examples","title":"Examples","text":"Here we showcase the telegraph model, a simplistic description of mRNA transcription in biological cells. We have one gene that transitions stochastically between an on and an off state and produces mRNA molecules while it is in the on state.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"This system technically consists of three different species, namely the two states of the gene and mRNA. It is clear, however, that these are not independent as D_on(t) + D_off(t) = 1. In order to solve the Chemical Master Equation we can therefore recover D_off(t) from the other variables and the entire state of the system is described by only two variables: D_on(t) and M(t), as well as the total number of genes, which is a constant equal to 1. The default index handler class DefaultIndexHandler does this for us automatically and maps the state of the system to a two-dimensional array. This showcases that we can often reduce the number of species in the system to make it easier to solve numerically.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"using FiniteStateProjection, DifferentialEquations\n\n@parameters r1 r2 r3 r4\nrs = @reaction_network begin\n    r1, G_on --> G_on + M\n    (r2, r3), G_on <--> G_off\n    r4, M --> 0\nend r1 r2 r3 r4\n\nsys = FSPSystem(rs)\n\n# There is one conserved quantity: G_on + G_off\ncons = conservedquantities([1,0,0], sys)\n\n# Parameters for our system\nps = [ 15.0, 0.25, 0.15, 1.0 ]\n\n# Since G_on + G_off = const. we do not have to model the latter separately\nu0 = zeros(2, 50)\nu0[1,1] = 1.0\n\nprob = convert(ODEProblem, DefaultIndexHandler(sys, 1), sys, u0, 10.0, (ps, cons))\nsol = solve(prob, Vern7(), atol=1e-6)","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"(Image: Visualisation)","category":"page"},{"location":"internal.html#Internal-API","page":"Internal API","title":"Internal API","text":"","category":"section"},{"location":"internal.html#index_handlers_internal","page":"Internal API","title":"Index Handlers","text":"","category":"section"},{"location":"internal.html#Index-Handler-Interface","page":"Internal API","title":"Index Handler Interface","text":"","category":"section"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"User-defined index handlers should inherit from AbstractIndexHandler and implement the following methods:","category":"page"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"getsubstitutions\nbuild_rhs_header\nsingleindices\npairedindices","category":"page"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"FiniteStateProjection.singleindices(::AbstractIndexHandler, arr)\nFiniteStateProjection.pairedindices\nFiniteStateProjection.getsubstitutions\nFiniteStateProjection.build_rhs_header(::AbstractIndexHandler, ::FSPSystem)","category":"page"},{"location":"internal.html#FiniteStateProjection.singleindices-Tuple{AbstractIndexHandler, Any}","page":"Internal API","title":"FiniteStateProjection.singleindices","text":"singleindices(idxhandler::AbstractIndexHandler, arr)\n\nReturns all indices I in arr. Defaults to CartesianIndices, but can be overloaded for arbitrary index handlers. \n\n\n\n\n\n","category":"method"},{"location":"internal.html#FiniteStateProjection.pairedindices","page":"Internal API","title":"FiniteStateProjection.pairedindices","text":"pairedindices(idxhandler::AbstractIndexHandler, arr, shift::CartesianIndex)\n\nReturns all pairs of indices (I .- shift, I) in arr.\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.getsubstitutions","page":"Internal API","title":"FiniteStateProjection.getsubstitutions","text":"getsubstitutions(idxhandler::AbstractIndexHandler, sys::FSPSystem; state_sym::Symbol)::Dict\n\nConstruct the map speciesname => expr that gives the species abundances in terms of the state variable state_sym. See NaiveIndexHandler for the default implementation.\n\nSee also: build_ratefuncs, build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_rhs_header-Tuple{AbstractIndexHandler, FSPSystem}","page":"Internal API","title":"FiniteStateProjection.build_rhs_header","text":"build_rhs_header(idxhandler::AbstractIndexHandler, sys::FSPSystem)::Expr\n\nReturn initialisation code for the RHS function, unpacking the parameters p supplied by DifferentialEquations. The default implementation just unpacks parameters from p.\n\nSee also: unpackparams, build_rhs\n\n\n\n\n\n","category":"method"},{"location":"internal.html#Built-in-implementations","page":"Internal API","title":"Built-in implementations","text":"","category":"section"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"elidedspecies(::AbstractMatrix{Int})\nFiniteStateProjection.elisions\nFiniteStateProjection.getsubstitutions(::NaiveIndexHandler, ::FSPSystem)\nFiniteStateProjection.getsubstitutions(::DefaultIndexHandler, ::FSPSystem)\nFiniteStateProjection.build_rhs_header(::DefaultIndexHandler, ::FSPSystem)\nFiniteStateProjection.pairedindices(::DefaultIndexHandler, ::AbstractArray, ::CartesianIndex)","category":"page"},{"location":"internal.html#FiniteStateProjection.elidedspecies-Tuple{AbstractMatrix{Int64}}","page":"Internal API","title":"FiniteStateProjection.elidedspecies","text":"elidedspecies(cons_laws::AbstractMatrix{Int})::Vector\n\nReturns a list of species  s_1   which can be removed from the reaction system description using the provided matrix of conservation laws.\n\n\n\n\n\n","category":"method"},{"location":"internal.html#FiniteStateProjection.elisions","page":"Internal API","title":"FiniteStateProjection.elisions","text":"elisions(idxhandler::DefaultIndexHandler, sys::FSPSystem)\n\nReplaces the symbols A(t), B(t), ... of elided species by N_1(t) - X(t) - Y(t), N_2(t) - U(t) - V(t), ..., where N_i(t) are the conserved quantities of the system.\n\nSee also: getsubstitutions\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.getsubstitutions-Tuple{NaiveIndexHandler, FSPSystem}","page":"Internal API","title":"FiniteStateProjection.getsubstitutions","text":"getsubstitutions(idxhandler::NaiveIndexHandler, sys::FSPSystem; state_sym::Symbol)::Dict\n\nDefines the abundance of species S_i to be state_sym[i] - offset.\n\n\n\n\n\n","category":"method"},{"location":"internal.html#FiniteStateProjection.getsubstitutions-Tuple{DefaultIndexHandler, FSPSystem}","page":"Internal API","title":"FiniteStateProjection.getsubstitutions","text":"getsubstitutions(idxhandler::DefaultIndexHandler, sys::FSPSystem; state_sym::Symbol)::Dict\n\nSimilar to its NaiveIndexHandler variant, but computes the abundances of elided species from the conserved quantities and the reduced species.\n\n\n\n\n\n","category":"method"},{"location":"internal.html#FiniteStateProjection.build_rhs_header-Tuple{DefaultIndexHandler, FSPSystem}","page":"Internal API","title":"FiniteStateProjection.build_rhs_header","text":"build_rhs_header(idxhandler::DefaultIndexHandler, sys::FSPSystem)::Expr\n\nAssumes p is of the form (params, cons::AbstractVector{Int}) where params  are the system parameters and cons the conserved quantities.\n\n\n\n\n\n","category":"method"},{"location":"internal.html#FiniteStateProjection.pairedindices-Tuple{DefaultIndexHandler, AbstractArray, CartesianIndex}","page":"Internal API","title":"FiniteStateProjection.pairedindices","text":"pairedindices(idxhandler::DefaultIndexHandler, arr::AbstractArray, shift::CartesianIndex)\n\nSimilar to its NaiveIndexHandler variant, but converts the indices into indices into the reduced state space array.\n\n\n\n\n\n","category":"method"},{"location":"internal.html#Function-Building","page":"Internal API","title":"Function Building","text":"","category":"section"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"FiniteStateProjection.build_rhs\nFiniteStateProjection.unpackparams\nFiniteStateProjection.build_ratefuncs\nFiniteStateProjection.build_rhs_firstpass\nFiniteStateProjection.build_rhs_secondpass","category":"page"},{"location":"internal.html#FiniteStateProjection.build_rhs","page":"Internal API","title":"FiniteStateProjection.build_rhs","text":"build_rhs(idxhandler::AbstractIndexHandler, sys::FSPSystem;\n          combinatoric_ratelaw::Bool)\n\nBuilds the function f(du,u,p,t) that defines the right-hand side of the CME,  for use in the ODE solver. If expression is true, returns an expression, else compiles the function. \n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.unpackparams","page":"Internal API","title":"FiniteStateProjection.unpackparams","text":"unpackparams(sys::FSPSystem, psym::Symbol)\n\nReturns code unpacking the parameters of the system from the symbol psym in the form (p1, p2, ...) = psym. This should be called in all overloads of build_rhs_header. It is assumed that the variable psym is an AbstractVector{Float64}.\n\nSee also: build_rhs_header, build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_ratefuncs","page":"Internal API","title":"FiniteStateProjection.build_ratefuncs","text":"build_ratefuncs(idxhandler::AbstractIndexHandler, sys::FSPSystem; \n                state_sym::Symbol, combinatoric_ratelaw::Bool)::Vector\n\nReturn the rate functions converted to Julia expressions in the state variable  state_sym. Abundances of the species are computed using getsubstitutions.\n\nSee also: getsubstitutions, build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_rhs_firstpass","page":"Internal API","title":"FiniteStateProjection.build_rhs_firstpass","text":"build_rhs_firstpass(sys::FSPSystem, rfs)::Expr\n\nReturn code for the first pass of the RHS function. Goes through all reactions and computes the negative part of the CME (probability flowing out of states). This is a simple array traversal and can be done in one go for all reactions.\n\nSee also: build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_rhs_secondpass","page":"Internal API","title":"FiniteStateProjection.build_rhs_secondpass","text":"build_rhs_secondpass(sys::FSPSystem, rfs)::Expr\n\nReturn code for the second pass of the RHS function. Goes through all reactions and computes the positive part of the CME (probability flowing into states). This requires accessing du and u at different locations depending on the net stoichiometries. In order to reduce  random memory access reactions are processed one by one.\n\nSee also: build_rhs\n\n\n\n\n\n","category":"function"},{"location":"mainapi.html#Main-API","page":"Main API","title":"Main API","text":"","category":"section"},{"location":"mainapi.html#Reaction-Systems","page":"Main API","title":"Reaction Systems","text":"","category":"section"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"The following describes the FSPSystem struct, a thin wrapper around ModelingToolkit's ReactionSystem for use with this package. The main functionality of this type is determining conserved quantities for a reaction system.","category":"page"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"FSPSystem\nconservationlaws(::FSPSystem)\nconservedquantities","category":"page"},{"location":"mainapi.html#FiniteStateProjection.FSPSystem","page":"Main API","title":"FiniteStateProjection.FSPSystem","text":"struct FSPSystem\n    rs::ReactionSystem\n    cons_laws::Matrix{Int}\nend\n\nThin wrapper around ModelingToolkit.ReactionSystem for use with this package.  Automatically computes a matrix of conservation laws (see  conservationlaws). \n\nConstructor: FSPSystem(rs::ReactionSystem)\n\n\n\n\n\n","category":"type"},{"location":"mainapi.html#FiniteStateProjection.conservationlaws-Tuple{FSPSystem}","page":"Main API","title":"FiniteStateProjection.conservationlaws","text":"conservationlaws(sys::FSPSystem)::Matrix{Int}\n\nReturns matrix of conservation laws associated with the system. Each row of  the matrix contains the stoichiometric coefficients of a different conserved quantity.\n\n\n\n\n\n","category":"method"},{"location":"mainapi.html#FiniteStateProjection.conservedquantities","page":"Main API","title":"FiniteStateProjection.conservedquantities","text":"conservedquantities(state, sys::FSPSystem)\n\nCompute conserved quantities for the system at the given state.\n\n\n\n\n\n","category":"function"},{"location":"mainapi.html#Creating-ODE-systems","page":"Main API","title":"Creating ODE systems","text":"","category":"section"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"The following methods are the main way to create a system of ODEs representing the time-dependent FSP. This package provides a flexible way to represent the FSP in memory via index handlers, see [Index Handlers] for more information. ","category":"page"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"Base.convert(::Type{ODEFunction}, ::AbstractIndexHandler, ::FSPSystem)\nBase.convert(::Type{ODEProblem}, ::AbstractIndexHandler, ::FSPSystem, u0, tmax, p)","category":"page"},{"location":"mainapi.html#Base.convert-Tuple{Type{ODEFunction}, AbstractIndexHandler, FSPSystem}","page":"Main API","title":"Base.convert","text":"convert(::Type{ODEFunction}, idxhandler::AbstractIndexHandler, sys::FSPSystem)\n\nReturn an ODEFunction defining the right-hand side of the CME.\n\nCombines the RHS func and its Jacobian to define an ODEFunction for  use with DifferentialEquations. This is where most of the work in the package happens; for best performance it is suggested to build an ODEFunction once for a given reaction system and reuse it instead of directly converting a reaction system to an ODEProblem (which implicitly calls this function).\n\n\n\n\n\n","category":"method"},{"location":"mainapi.html#Base.convert-Tuple{Type{ODEProblem}, AbstractIndexHandler, FSPSystem, Any, Any, Any}","page":"Main API","title":"Base.convert","text":"convert(::Type{ODEProblem}, idxhandler::AbstractIndexHandler, sys::FSPSystem, u0, tmax, p)\n\nReturn an ODEProblem for use in DifferentialEquations. This function implicitly callsconvert(ODEFunction, indexhandler, sys). It is usually more efficient to create anODEFunctionfirst and then use that to createODEProblem`s.\n\n\n\n\n\n","category":"method"},{"location":"index.html#FiniteStateProjection.jl-Documentation","page":"Home","title":"FiniteStateProjection.jl Documentation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"CurrentModule=FiniteStateProjection","category":"page"},{"location":"index.html#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"FiniteStateProjection.jl is a package that implements Finite State Projection algorithms for chemical reaction networks based on Catalyst.jl and ModelingToolkit. FiniteStateProjection.jl converts descriptions of reaction networks into ODEProblems that can be used to compute approximate solutions of the Chemical Master Equation with packages such as DifferentialEquations.jl.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"FiniteStateProjection.jl works by converting a ReactionSystem into a function that computes the right-hand side of the Chemical Master Equation:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"fracmathrmdmathrmd t P(t) = A P(t)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"This function is generated dynamically as an ODEFunction for use with DifferentialEquations.jl and specialised for each ReactionSystem. Users can use their preferred array types and provide additional features by overloading the functions in this package.","category":"page"},{"location":"index.html#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Flexible API for user-defined array types via IndexHandlers\nAutomatic dimensionality reduction for systems with conserved quantities\nOn-the-fly generation of specialised functions for improved performance","category":"page"},{"location":"indexhandlers.html#Index-Handlers","page":"Index Handlers","title":"Index Handlers","text":"","category":"section"},{"location":"indexhandlers.html","page":"Index Handlers","title":"Index Handlers","text":"The task of an index handler is to provide a mapping between the system state and the way it is stored in memory, usually as a multidimensional array. The standard approach is to represent the states of a system with s reactions as an s-dimensional array and have the index (i_1 ldots i_s) correspond to the state (n_1 = i_1 ldots n_s = i_s). This is implemented by the class NaiveIndexHandler, which accepts an offset argument to deal with Julia's 1-based indexing (so the Julia index (1ldots1) corresponds to the state with no molecules). For systems with conservation laws the DefaultIndexHandler class generally stores the data more efficiently.","category":"page"},{"location":"indexhandlers.html","page":"Index Handlers","title":"Index Handlers","text":"See the internal API on how to define your own IndexHandler type.","category":"page"},{"location":"indexhandlers.html","page":"Index Handlers","title":"Index Handlers","text":"NaiveIndexHandler\nDefaultIndexHandler\nreducedspecies\nelidedspecies(::DefaultIndexHandler)","category":"page"},{"location":"indexhandlers.html#FiniteStateProjection.NaiveIndexHandler","page":"Index Handlers","title":"FiniteStateProjection.NaiveIndexHandler","text":"struct NaiveIndexHandler <: AbstractIndexHandler\n    offset::Int\nend\n\nBasic index handler that stores the state of a system with s species in an s-dimensional array. The offset parameter denotes the offset by which the array is indexed (defaults to 1 in Julia). Use OffsetArrays.jl to enable 0-based indexing.\n\nThis is the simplest index handler, but it will not be optimal if some states cannot be reached from the initial state, e.g. due to the presence of conservation laws. It is generally better to use DefaultIndexHandler, which will automatically elide species where possible.\n\nConstructors: NaiveIndexHandler([sys::FSPSystem, offset::Int=1])\n\nSee also: DefaultIndexHandler\n\n\n\n\n\n","category":"type"},{"location":"indexhandlers.html#FiniteStateProjection.DefaultIndexHandler","page":"Index Handlers","title":"FiniteStateProjection.DefaultIndexHandler","text":"struct DefaultIndexHandler <: AbstractIndexHandler\n\nMore efficient index handler that improves upon NaiveIndexHandler by eliminating variables whose abundances can be computed from other variables using conservation laws. Describes the system using a subset of the original species which can be obtained via reducedspecies. Reduces the  dimensionality of the FSP by the number of conservation laws in the system.\n\nConstructors: DefaultIndexHandler(sys::FSPSystem[, offset::Int=1])\n\nSee also: reducedspecies, elidedspecies, NaiveIndexHandler\n\n\n\n\n\n","category":"type"},{"location":"indexhandlers.html#FiniteStateProjection.reducedspecies","page":"Index Handlers","title":"FiniteStateProjection.reducedspecies","text":"reducedspecies(idxhandler::DefaultIndexHandler)\n\nReturn indices of reduced species.\n\nSee also: elidedspecies\n\n\n\n\n\n","category":"function"},{"location":"indexhandlers.html#FiniteStateProjection.elidedspecies-Tuple{DefaultIndexHandler}","page":"Index Handlers","title":"FiniteStateProjection.elidedspecies","text":"elidedspecies(idxhandler::DefaultIndexHandler)\n\nReturn indices of elided species.\n\nSee also: reducedspecies\n\n\n\n\n\n","category":"method"}]
}
