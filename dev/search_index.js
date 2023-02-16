var documenterSearchIndex = {"docs":
[{"location":"matrix.html#matrix_conversions","page":"Matrix Conversions","title":"Matrix Conversions","text":"","category":"section"},{"location":"matrix.html","page":"Matrix Conversions","title":"Matrix Conversions","text":"CurrentModule = FiniteStateProjection","category":"page"},{"location":"matrix.html","page":"Matrix Conversions","title":"Matrix Conversions","text":"FiniteStateProjection.jl provides functionality for building the right-hand side of the CME as a (sparse) matrix. This provides another way to solve the CME in time:","category":"page"},{"location":"matrix.html","page":"Matrix Conversions","title":"Matrix Conversions","text":"...\n\nA = convert(SparseMatrixCSC, sys, dims, p, 0)\n\nprob = ODEProblem((du,u,p,t) -> mul!(du, p, u), u0, tt, A)\n\n...","category":"page"},{"location":"matrix.html","page":"Matrix Conversions","title":"Matrix Conversions","text":"This can also be done for steady-state problems:","category":"page"},{"location":"matrix.html","page":"Matrix Conversions","title":"Matrix Conversions","text":"...\n\nA = convert(SparseMatrixCSC, sys, dims, p, SteadyState())\n\nprob = SteadyStateProblem((du,u,p,t) -> mul!(vec(du), p, vec(u)), u0, A)\n\n...","category":"page"},{"location":"matrix.html","page":"Matrix Conversions","title":"Matrix Conversions","text":"Note that the matrix A has to be rebuilt for every truncation size and every set of parameters, a restriction not shared by the ODEFunction API.","category":"page"},{"location":"matrix.html","page":"Matrix Conversions","title":"Matrix Conversions","text":"Base.convert(::Type{SparseMatrixCSC}, ::FSPSystem, ::NTuple, ps, t::Real)\nBase.convert(::Type{SparseMatrixCSC}, ::FSPSystem, ::NTuple, ps, ::SteadyState)\nvec","category":"page"},{"location":"matrix.html#Base.convert-Tuple{Type{SparseMatrixCSC}, FSPSystem, Tuple{Vararg{T, N}} where {N, T}, Any, Real}","page":"Matrix Conversions","title":"Base.convert","text":"Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, dims::NTuple, ps, t::Real)\n\nConvert the reaction system into a sparse matrix defining the right-hand side of the  Chemical Master Equation. dims is a tuple denoting the dimensions of the FSP and  ps is the tuple of parameters. The sparse matrix works on the flattened version of the state obtained using vec.\n\n\n\n\n\n","category":"method"},{"location":"matrix.html#Base.convert-Tuple{Type{SparseMatrixCSC}, FSPSystem, Tuple{Vararg{T, N}} where {N, T}, Any, SteadyState}","page":"Matrix Conversions","title":"Base.convert","text":"Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, dims::NTuple, ps, ::SteadyState)\n\nConvert the reaction system into a sparse matrix defining the right-hand side of the  Chemical Master Equation, steady-state version.\n\n\n\n\n\n","category":"method"},{"location":"matrix.html#Base.vec","page":"Matrix Conversions","title":"Base.vec","text":"vec(idxhandler::AbstractIndexHandler, arr)\n\nConverts the right-hand side defining the solution of the CME into a one-dimensional vector to which a matrix can be applied.\n\nSee also: LinearIndices\n\n\n\n\n\n","category":"function"},{"location":"examples.html#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples.html#Birth-Death-Model","page":"Examples","title":"Birth-Death Model","text":"","category":"section"},{"location":"examples.html","page":"Examples","title":"Examples","text":"This example models a linear birth-death process. The reaction network is easily defined using Catalyst.jl. Our truncated state space has length 50, which is enough for this simple system.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"using FiniteStateProjection\nusing OrdinaryDiffEq\n\nrn = @reaction_network begin\n    σ, 0 --> A\n    d, A --> 0\nend σ d\n\nsys = FSPSystem(rn)\n\n# Parameters for our system\nps = [ 10.0, 1.0 ]\n\n# Initial distribution (over 1 species)\n# Here we start with 0 copies of A\nu0 = zeros(50)\nu0[1] = 1.0 \n\nprob = ODEProblem(sys, u0, (0, 10.0), ps)\nsol = solve(prob, Vern7())","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"(Image: Visualisation)","category":"page"},{"location":"examples.html#Telegraph-Model","page":"Examples","title":"Telegraph Model","text":"","category":"section"},{"location":"examples.html","page":"Examples","title":"Examples","text":"Here we showcase the telegraph model, a simplistic description of mRNA transcription in biological cells. We have one gene that transitions stochastically between an on and an off state and produces mRNA molecules while it is in the on state.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"The most straightforward description of this system includes three species: two gene states, G_on and G_off, and mRNA M. The state space for this system is 3-dimensional. We know, however, that G_on and G_off never occur at the same time, indeed the conservation law [G_on] + [G_off] = 1 allows us to express the state of the system in terms of G_on and M only. The state space of this reduced system is 2-dimensional.If we use an mRNA cutoff of 100, the state space for the original model has size 2 times 2 times 100 = 400, while the reduced state space has size 2 times 100 = 200, a two-fold saving. Since the FSP get computationally more expensive for each species in a system, eliminating redundant species as above is recommended for improved performance.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"note: Note\nThe class ReducingIndexHandler, which performed such reduction on the fly, has been deprecated and will be moved into Catalyst.jl.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"using FiniteStateProjection\nusing OrdinaryDiffEq\n\nrn = @reaction_network begin\n    σ_on * (1 - G_on), 0 --> G_on\n    σ_off, G_on --> 0\n    ρ, G_on --> G_on + M\n    d, M --> 0\nend σ_on σ_off ρ d\n\nsys = FSPSystem(rn)\n\n# Parameters for our system\nps = [ 0.25, 0.15, 15.0, 1.0 ]\n\n# Initial distribution (over two species)\n# Here we start with 0 copies of G_on and M\nu0 = zeros(2, 50)\nu0[1,1] = 1.0\n\nprob = ODEProblem(sys, u0, (0, 10.0), ps)\nsol = solve(prob, Vern7())","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"(Image: Visualisation)","category":"page"},{"location":"troubleshoot.html#troubleshoot","page":"Troubleshooting","title":"Troubleshooting","text":"","category":"section"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"Solving the Chemical Master Equation numerically is a difficult task and errors are liable to crop up. The following section presents some ways to catch common errors.","category":"page"},{"location":"troubleshoot.html#Ensure-your-state-space-has-the-right-dimension","page":"Troubleshooting","title":"Ensure your state space has the right dimension","text":"","category":"section"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"If your are solving an SIR model with three species, S, I and R, your state space will be 3-dimensional. FiniteStateProjection.jl computes probabilities for all states simultaneously and stores the results in a 3-dimensional array. In particular, u0 must have type Float64 or similar as it represents numbers between 0 and 1.","category":"page"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"# correct\nu0 = zeros(101, 101, 101)\nu0[100,2,1] = 1.0           # start with 1 infected and 99 susceptible individuals\n\n# incorrect\nu0 = [99,1,0]               # wrong type (Int) and shape (1D)","category":"page"},{"location":"troubleshoot.html#Ensure-your-state-space-is-big-enough","page":"Troubleshooting","title":"Ensure your state space is big enough","text":"","category":"section"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"A common reason for the solver returning nonsensical solutions is a state space that is too small. Since the Finite State Projection works with a finite-dimensional approximation of the system, the number of states considered can have a large impact on accuracy. The loss of accuracy due to using a smaller state space is the truncation error.","category":"page"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"A good way to check whether your state space is large enough is to solve the CME until the required time t and to sum up the probabilities for each state - this gives the probability that the system will have remained in the truncated state space from start to end. If this quantity is noticeably less than 1, the state space is likely too small. As an informal rule of thumb, a value of less than 95% indicates that the solution will not be reliable.","category":"page"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"# always true\nsum(u0) == 1                # our initial condition lies in the truncated state space\n\n# good\nsum(ut) >= 0.99             # our truncation covers most of the relevant states\n\n# bad\nsum(ut) < 0.95              # we do not have enough coverage","category":"page"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"The above does not work directly when computing steady-state probabilities as the value will usually drop to 0 for large enough t. In this case it is a good idea to redo the computation with a larger state space - the results should agree if the truncation error is small.","category":"page"},{"location":"troubleshoot.html#Ensure-your-propensities-are-positive","page":"Troubleshooting","title":"Ensure your propensities are positive","text":"","category":"section"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"This point might seem obvious, but errors in the rate functions, or an incorrectly chosen truncation, can lead to negative reaction propensities that will typically result in numerical instabilities. As an example, consider the following version of the SI model where the population size S + I = N is constant, allowing us to rewrite the system using only one species I (with S = N - I):","category":"page"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"rn = @reaction_network begin\n   σ * (N - I), I --> 2I\n   ρ, I --> 0\nend σ ρ N\n\nsys_fsp = FSPSystem(rn)","category":"page"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"Here the propensity function for the first reaction will negative if I  N, so the following may result in numerical instabilities:","category":"page"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"u0 = zeros(30)\nu0[2] = 1\n\n# N is too small for the state space!\nprob_fsp = convert(ODEProblem, sys_fsp, u0, (0, 100.), [ 1., 1., 20 ])","category":"page"},{"location":"troubleshoot.html#Ensure-you-are-using-the-right-solver","page":"Troubleshooting","title":"Ensure you are using the right solver","text":"","category":"section"},{"location":"troubleshoot.html","page":"Troubleshooting","title":"Troubleshooting","text":"The Chemical Master Equation is generally very stiff and requires a solver that can handle this stiffness, see Tips & Tricks. If your solver fails, first check if any of the above points apply. You may be able to get a different solver to work; this requires some experimentation. Anecdotally some systems, particularly oscillatory ones such as the Schlögl model, can pose significant challenges to most solvers and take inordinate amounts of time to solve. I am not aware of any solutions for this at the moment, but please consider opening an issue on GitHub if you encounter examples of this sort.","category":"page"},{"location":"internal.html#Internal-API","page":"Internal API","title":"Internal API","text":"","category":"section"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"CurrentModule = FiniteStateProjection","category":"page"},{"location":"internal.html#index_handlers_internal","page":"Internal API","title":"Index Handlers","text":"","category":"section"},{"location":"internal.html#Index-Handler-Interface","page":"Internal API","title":"Index Handler Interface","text":"","category":"section"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"User-defined index handlers should inherit from AbstractIndexHandler and implement the following methods:","category":"page"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"getsubstitutions\nbuild_rhs_header\nsingleindices\npairedindices","category":"page"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"For matrix conversions they should additionally implement:","category":"page"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"LinearIndices\nvec","category":"page"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"singleindices\npairedindices\ngetsubstitutions\nbuild_rhs_header\nLinearIndices","category":"page"},{"location":"internal.html#FiniteStateProjection.singleindices","page":"Internal API","title":"FiniteStateProjection.singleindices","text":"singleindices(idxhandler::AbstractIndexHandler, arr)\n\nReturns all indices I in arr. Defaults to CartesianIndices, but can be overloaded for arbitrary index handlers.\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.pairedindices","page":"Internal API","title":"FiniteStateProjection.pairedindices","text":"pairedindices(idxhandler::AbstractIndexHandler, arr, shift::CartesianIndex)\n\nReturns all pairs of indices (I .- shift, I) in arr.\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.getsubstitutions","page":"Internal API","title":"FiniteStateProjection.getsubstitutions","text":"getsubstitutions(idxhandler::AbstractIndexHandler, rs::ReactionSystem; state_sym::Symbol)\n\nReturns a dict of the form S_i => f_i(state_sym), where each f_i is an expression for the abundance of species S_i in terms of the state variable state_sym.\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_rhs_header","page":"Internal API","title":"FiniteStateProjection.build_rhs_header","text":"build_rhs_header(sys::FSPSystem)\n\nReturn initialisation code for the RHS function, unpacking the parameters p supplied by DifferentialEquations. The default implementation just unpacks parameters from p.\n\nSee also: unpackparams, build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#Base.LinearIndices","page":"Internal API","title":"Base.LinearIndices","text":"LinearIndices(idxhandler::AbstractIndexHandler, arr)\n\nReturns an object lind which converts indices returned from singleindices and pairedindices to linear indices compatible with vec via lind[idx_cart] = idx_lin. The indices are related via\n\narr[idx_cart] == vec(idxhandler, arr)[idx_lin]\n\nSee also: vec\n\n\n\n\n\n","category":"type"},{"location":"internal.html#Built-in-implementations","page":"Internal API","title":"Built-in implementations","text":"","category":"section"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"getsubstitutions(::DefaultIndexHandler, ::ReactionSystem)","category":"page"},{"location":"internal.html#FiniteStateProjection.getsubstitutions-Tuple{DefaultIndexHandler, ReactionSystem}","page":"Internal API","title":"FiniteStateProjection.getsubstitutions","text":"getsubstitutions(sys::FSPSystem{DefaultIndexHandler}; state_sym::Symbol)::Dict\n\nDefines the abundance of species S_i to be state_sym[i] - offset.\n\n\n\n\n\n","category":"method"},{"location":"internal.html#Function-Building","page":"Internal API","title":"Function Building","text":"","category":"section"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"build_rhs\nunpackparams\nbuild_ratefuncs\nbuild_rhs_firstpass\nbuild_rhs_secondpass","category":"page"},{"location":"internal.html#FiniteStateProjection.build_rhs","page":"Internal API","title":"FiniteStateProjection.build_rhs","text":"build_rhs(sys::FSPSystem)\n\nBuilds the function f(du,u,p,t) that defines the right-hand side of the CME for use with DifferentialEquations.jl.\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.unpackparams","page":"Internal API","title":"FiniteStateProjection.unpackparams","text":"unpackparams(sys::FSPSystem, psym::Symbol)\n\nReturns code unpacking the parameters of the system from the symbol psym in the form (p1, p2, ...) = psym. This should be called in all overloads of build_rhs_header. It is assumed that the variable psym represents an AbstractVector.\n\nSee also: build_rhs_header, build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_ratefuncs","page":"Internal API","title":"FiniteStateProjection.build_ratefuncs","text":"build_ratefuncs(rs, ih; state_sym::Symbol, combinatoric_ratelaw::Bool)::Vector\n\nReturn the rate functions converted to Julia expressions in the state variable  state_sym. Abundances of the species are computed using getsubstitutions.\n\nSee also: getsubstitutions, build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_rhs_firstpass","page":"Internal API","title":"FiniteStateProjection.build_rhs_firstpass","text":"build_rhs_firstpass(sys::FSPSystem)\n\nReturn code for the first pass of the RHS function, for the time-dependent FSP. Goes through all reactions and computes the negative part of the CME (probability flowing out of states). This is a simple array traversal and can be done in one go for all reactions.\n\nSee also: build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_rhs_secondpass","page":"Internal API","title":"FiniteStateProjection.build_rhs_secondpass","text":"build_rhs_secondpass(sys::FSPSystem)\n\nReturn code for the second pass of the RHS function. Goes through all reactions and computes the positive part of the CME (probability flowing into states). This requires accessing du and u at different locations depending on the net stoichiometries. In order to reduce random memory access reactions are processed one by one.\n\nSee also: build_rhs\n\n\n\n\n\n","category":"function"},{"location":"internal.html#Steady-State-Functions","page":"Internal API","title":"Steady-State Functions","text":"","category":"section"},{"location":"internal.html","page":"Internal API","title":"Internal API","text":"build_rhs_ss\nbuild_rhs_singlepass_ss","category":"page"},{"location":"internal.html#FiniteStateProjection.build_rhs_ss","page":"Internal API","title":"FiniteStateProjection.build_rhs_ss","text":"build_rhs_ss(sys::FSPSystem)\n\nBuilds the function f(du,u,p,t) that defines the right-hand side of the CME for use with SteadyStateProblems.\n\n\n\n\n\n","category":"function"},{"location":"internal.html#FiniteStateProjection.build_rhs_singlepass_ss","page":"Internal API","title":"FiniteStateProjection.build_rhs_singlepass_ss","text":"build_rhs_singlepass_ss(sys::FSPSystem)\n\nReturn code for the RHS function in a single pass, for use with steady state problems. Transitions out of the truncated state space are ignored.\n\nSee also: build_rhs_ss\n\n\n\n\n\n","category":"function"},{"location":"mainapi.html#Main-API","page":"Main API","title":"Main API","text":"","category":"section"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"CurrentModule = FiniteStateProjection","category":"page"},{"location":"mainapi.html#Reaction-Systems","page":"Main API","title":"Reaction Systems","text":"","category":"section"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"This section describes the FSPSystem struct, a thin wrapper around Catalyst's ReactionSystem for use with this package.","category":"page"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"FSPSystem","category":"page"},{"location":"mainapi.html#FiniteStateProjection.FSPSystem","page":"Main API","title":"FiniteStateProjection.FSPSystem","text":"Thin wrapper around Catalyst.ReactionSystem for use with this package. \n\nConstructor: FSPSystem(rs::ReactionSystem[, ih=DefaultIndexHandler(); combinatoric_ratelaw::Bool=true])\n\n\n\n\n\n","category":"type"},{"location":"mainapi.html#Creating-ODE-systems","page":"Main API","title":"Creating ODE systems","text":"","category":"section"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"The following methods convert a reaction network into a system of ODEs representing the time-dependent FSP. This package provides a flexible way to represent the FSP in memory via index handlers, see [Index Handlers] for more information. ","category":"page"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"Base.convert(::Type{ODEFunction}, ::FSPSystem)\nBase.convert(::Type{ODEProblem}, ::FSPSystem, u0, tmax, p)","category":"page"},{"location":"mainapi.html#Base.convert-Tuple{Type{ODEFunction}, FSPSystem}","page":"Main API","title":"Base.convert","text":"Base.convert(::Type{ODEFunction}, sys::FSPSystem)\n\nReturn an ODEFunction defining the right-hand side of the CME.\n\nCreates an ODEFunction for use with DifferentialEquations. This is where most of the work in the package happens; for best performance it is suggested to build an ODEFunction once for a given reaction system and reuse it instead of directly converting a reaction system to an ODEProblem (which implicitly calls this function).\n\n\n\n\n\n","category":"method"},{"location":"mainapi.html#Base.convert-Tuple{Type{ODEProblem}, FSPSystem, Any, Any, Any}","page":"Main API","title":"Base.convert","text":"Base.convert(::Type{ODEProblem}, sys::FSPSystem, u0, tmax[, p])\n\nReturn an ODEProblem for use in DifferentialEquations. This function implicitly calls convert(ODEFunction, sys). It is usually more efficient to create an ODEFunction  first and then use that to create ODEProblems.\n\n\n\n\n\n","category":"method"},{"location":"mainapi.html#Steady-State-Problems","page":"Main API","title":"Steady-State Problems","text":"","category":"section"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"Computing steady-state distributions can be done using the SteadyStateDiffEq.jl package. At the moment FiniteStateProjection.jl adjusts the rate matrix so that reactions leaving the truncated state space have propensity 0.","category":"page"},{"location":"mainapi.html","page":"Main API","title":"Main API","text":"Base.convert(::Type{ODEFunction}, ::FSPSystem, ::SteadyState)\nBase.convert(::Type{SteadyStateProblem}, ::FSPSystem, u0, p)","category":"page"},{"location":"mainapi.html#Base.convert-Tuple{Type{ODEFunction}, FSPSystem, SteadyState}","page":"Main API","title":"Base.convert","text":"Base.convert(::Type{ODEFunction}, sys::FSPSystem, ::SteadyState)\n\nReturn an ODEFunction defining the right-hand side of the CME, for use with SteadyStateProblems.\n\n\n\n\n\n","category":"method"},{"location":"mainapi.html#Base.convert-Tuple{Type{SteadyStateProblem}, FSPSystem, Any, Any}","page":"Main API","title":"Base.convert","text":"Base.convert(::Type{SteadyStateProblem}, sys::FSPSystem, u0[, p])\n\nReturn a SteadyStateProblem for use in `DifferentialEquations.\n\n\n\n\n\n","category":"method"},{"location":"index.html#FiniteStateProjection.jl-Documentation","page":"Home","title":"FiniteStateProjection.jl Documentation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"CurrentModule=FiniteStateProjection","category":"page"},{"location":"index.html#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"FiniteStateProjection.jl is a package that implements Finite State Projection algorithms for chemical reaction networks based on Catalyst.jl and ModelingToolkit. FiniteStateProjection.jl converts descriptions of reaction networks into ODEProblems that can be used to compute approximate solutions of the Chemical Master Equation with packages such as DifferentialEquations.jl.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The Chemical Master Equation allows us to compute the probability of a reaction system to be in a given state vec n at a time t. The state of a system with s species is described by a vector vec n = (n_1 n_2 ldots n_s) of length s, where n_1 n_2 ldots are the abundances of the different species. The reaction system itself can be represented by an s-dimensional array Pn_1ldotsn_s that describes the probability of the system being in state (n_1ldotsn_s). Here P depends on time as the system evolves. Given initial conditions P(0), the time evolution of P is given by the Chemical Master Equation, which is a set of linear ODEs that can be solved numerically to predict the probability distribution of the system in the future. ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Since the number of states for a system can be infinite, the Finite State Projection approximates the Chemical Master Equation by only considering a finite number of states, namely those which are most probable. There are various ways to truncate the state space, but the most common is to provide an upper limit for each species, that is, to require n_1  M_1, n_2  M_2, etc. for fixed thresholds M_1 M_2 ldots. The number of states considered will be M_1 times M_2 times ldots times M_s, and P can be represented as an array with those dimensions. While truncating the CME allows us to solve the equations numerically, the amount of space required to store P increases quickly as more species are added. ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"FiniteStateProjection.jl works by converting a ReactionSystem into a function that computes the right-hand side of the (truncated) Chemical Master Equation:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"fracmathrmdmathrmd t P(t) = A P(t)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"This function is generated dynamically as an ODEFunction for use with DifferentialEquations.jl and specialised for each ReactionSystem. Users can use their preferred array types and provide additional features by overloading the functions in this package. Alternatively the matrix A can be constructed as a SparseMatrixCSC. FiniteStateProjection.jl supports both the time-dependent Chemical Master Equation and its steady-state version.","category":"page"},{"location":"index.html#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Built on top of Catalyst.jl\nFSP equations are generated as ODEFunction/ODEProblems and can be solved with DifferentialEquations.jl, with on-the-fly generation of targeted functions for improved performance\nThe Chemical Master Equation can be generated as a SparseMatrixCSC\nFlexible API for user-defined array types via IndexHandlers","category":"page"},{"location":"index.html#Acknowledgments","page":"Home","title":"Acknowledgments","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Special thanks to Xiamong Fu, Brian Munsky and Huy Vo for their examples and suggestions and for contributing most of the content in the Tips & Tricks page!","category":"page"},{"location":"indexhandlers.html#Index-Handlers","page":"Index Handlers","title":"Index Handlers","text":"","category":"section"},{"location":"indexhandlers.html","page":"Index Handlers","title":"Index Handlers","text":"The task of an index handler is to provide a mapping between the system state and the way it is stored in memory, usually as a multidimensional array. The standard approach is to represent the states of a system with s reactions as an s-dimensional array and have the index (i_1 ldots i_s) correspond to the state (n_1 = i_1 ldots n_s = i_s). This is implemented by the class DefaultIndexHandler, which accepts an offset argument to deal with Julia's 1-based indexing (so the Julia index (1ldots1) corresponds to the state with no molecules). ","category":"page"},{"location":"indexhandlers.html","page":"Index Handlers","title":"Index Handlers","text":"See the internal API on how to define your own IndexHandler type.","category":"page"},{"location":"indexhandlers.html","page":"Index Handlers","title":"Index Handlers","text":"DefaultIndexHandler","category":"page"},{"location":"indexhandlers.html#FiniteStateProjection.DefaultIndexHandler","page":"Index Handlers","title":"FiniteStateProjection.DefaultIndexHandler","text":"struct DefaultIndexHandler <: AbstractIndexHandler\n    offset::Int\nend\n\nBasic index handler that stores the state of a system with s species in an s-dimensional array. The offset parameter denotes the offset by which the array is indexed (defaults to 1 in Julia).\n\nThis is the simplest index handler, but it will not be optimal if some states cannot be reached from the initial state, e.g. due to the presence of conservation laws. In these cases one should use ReducingIndexHandler, which will automatically elide species where possible.\n\nConstructors: DefaultIndexHandler([sys::FSPSystem, offset::Int=1])\n\n\n\n\n\n","category":"type"},{"location":"tips.html#tips","page":"Tips & Tricks","title":"Tips and Tricks","text":"","category":"section"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"The FSP approximates an infinite-dimensional system of equations by truncating it to a finite number of variables. The accuracy of the FSP therefore depends on how many variables are retained, ie.~what portion of the state space is modelled. While simple reaction networks with 1 or 2 species are not too difficult to handle using the FSP, a naive approach will require unfeasibly large truncations for systems of even moderate complexity.","category":"page"},{"location":"tips.html#Solving-Linear-Equations","page":"Tips & Tricks","title":"Solving Linear Equations","text":"","category":"section"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"The Chemical Master Equation is a set of linear ODEs, and there is a vast literature on efficiently solving these. FiniteStateProjection.jl currently offers two main approaches:","category":"page"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"Convert a reaction system into an ODEFunction for use with DifferentialEquations.jl\nConvert a reaction system into a sparse matrix for use with any linear ODE solver","category":"page"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"While there are no official benchmarks yet the first method is likely going to be more efficient for systems with time-dependent reaction rates, while the second is recommended for time-homogeneous systems.","category":"page"},{"location":"tips.html#Picking-a-Solver","page":"Tips & Tricks","title":"Picking a Solver","text":"","category":"section"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"The FSP generally results in stiff equations, so using a stiff solver in generally is recommended. At the time of writing it seems like the best all-round solver is CVODE_BDF with the GMRES linear solver (requires Sundials.jl):","category":"page"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"using Sundials\n\n(...)\n\nsol = solve(prob, CVODE_BDF(linear_solver=:GMRES), atol=1e-6) # very fast","category":"page"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"For time-homogeneous systems where the right-hand side of the CME does not depend on t, solving the FSP is equivalent to computing the exponential-vector product exp(A .* t) * vec(u0), where A is the evolution operator. Here we use the vec function to flatten the state space into a vector that can be multiplied by the matrix A. We can compute the solution very efficiently using the expmv function provided by Expokit.jl:","category":"page"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"using Expokit\n\n(...)\n\nut = similar(u0)\nexpmv!(vec(ut), t, A, vec(u0), tol=1e-6)  # really fast","category":"page"},{"location":"tips.html#Further-Comments","page":"Tips & Tricks","title":"Further Comments","text":"","category":"section"},{"location":"tips.html","page":"Tips & Tricks","title":"Tips & Tricks","text":"Choosing the right solver for large systems of ODEs can result in time savings on the order of 10-100x, and it is recommended that you experiment with a few solvers to see which works best in your case. This section is still work in progress and there has been a lot of research on accelerating the FSP and extending it to larger reaction networks which will hopefully be reviewed here soon. Feel free to share any comments or suggestions in this direction at the Github repository!","category":"page"}]
}
