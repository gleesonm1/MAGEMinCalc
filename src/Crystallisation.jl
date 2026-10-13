using MAGEMin_C
using PythonCall
using DataFrames

function findliq(; bulk :: Vector{Float64}, P_kbar :: Float64, T_start_C :: Float64 = 1400.0,
                fo2_buffer :: Union{String, Nothing} = nothing, fo2_offset :: Union{Float64, Nothing} = 0.0, Model :: String = "ig",
                suppress :: Union{Vector{String},Nothing} = nothing)
    bulk_in = bulk;

    if fo2_offset === nothing
        fo2_offset = 0.0
    end

	new_bulk = bulk/sum(bulk);
    if Model === "ig"
    	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];
    else
    	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"];
    end

    T = T_start_C;

    if suppress !== nothing
        rm_list = remove_phases(suppress, Model)
    end
    
    data = fo2_buffer !== nothing ? Initialize_MAGEMin(Model, verbose=false, buffer=fo2_buffer) : Initialize_MAGEMin(Model, verbose=false)
    if fo2_buffer !== nothing
        if suppress !== nothing
            out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, B = fo2_offset,name_solvus=true )
        else
            out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset,name_solvus=true )
        end
    else
        if suppress !== nothing
            out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list,name_solvus=true )
        else
            out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt",name_solvus=true )
        end
    end

    Liq = ["liq", "fl"];
    PhaseList = out.ph;
    if fo2_buffer !== nothing
        filter!(x -> x != fo2_buffer, PhaseList)
    end 

    i = intersect(Liq, PhaseList);

    Step = [3, 1, 0.1];
    for k in eachindex(Step)
        if length(i) < length(PhaseList)
            while length(i) < length(PhaseList)
                T = T + Step[k]
                if fo2_buffer !== nothing
                    if suppress !== nothing
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, B = fo2_offset,name_solvus=true )
                    else
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset,name_solvus=true )
                    end
                else
                    if suppress !== nothing
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list,name_solvus=true )
                    else
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt",name_solvus=true )
                    end
                end

                PhaseList = out.ph
                if fo2_buffer !== nothing
                    filter!(x -> x != fo2_buffer, PhaseList)
                end 
                i = intersect(Liq, PhaseList)
            end
        end

        if length(i) === length(PhaseList)
            while length(i) === length(PhaseList)
                T = T - Step[k]
                if fo2_buffer !== nothing
                    if suppress !== nothing
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, B = fo2_offset,name_solvus=true )
                    else
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset,name_solvus=true )
                    end
                else
                    if suppress !== nothing
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list,name_solvus=true )
                    else
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt",name_solvus=true )
                    end
                end

                PhaseList = out.ph
                if fo2_buffer !== nothing
                    filter!(x -> x != fo2_buffer, PhaseList)
                end 
                i = intersect(Liq, PhaseList)
            end
        end

        T = T + Step[k]/2
    end

    Finalize_MAGEMin(data)
    T_Liq = T
    return T_Liq
end

"""
    findLiq_julia(; bulk, P_kbar, T_start_C, Model="ig", fo2_buffer=nothing,
                  fo2_offset=0.0, suppress=nothing, bracket_step_C=25.0,
                  max_bracket_expansions=8, max_bisection_iters=15, T_tol_C=0.25)

Vectorised liquidus-temperature search for a BATCH of points, each with its
own bulk composition, pressure, and starting temperature guess. Designed to
be called directly from Python (PetThermoTools's `findLiq_multi`) for
batches of 1000s of real compositions, AND (later) from MAGEMinCalc's own
`path`/`path_main` for a single isolated liquidus lookup -- N=1 is simply a
length-1 vector, so there is no special-casing anywhere below for the
single-point vs batch case; both are the exact same code path.

Deliberately separate from `findliq` (above: the existing, slow, single-
point, fixed-1400C step-search) and from the existing, broken, unused
`findLiq_multi(bulk, P_kbar, T_start_C)` (below: references `Model`/
`fo2_buffer` without taking them as parameters) -- neither is touched or
replaced by this function. See PetThermoTools/MAINTENANCE_NOTES.md
(2026-10-10, "MAGEMin liquidus search") for the design rationale and
verification results.

Uses `multi_point_minimization`, which natively supports a DIFFERENT bulk
composition per point (its own docstring's "Example 3") and parallelises
internally across points via Julia threads, all within ONE shared
`Initialize_MAGEMin`/`Finalize_MAGEMin` context for the whole batch --
unlike alphaMELTS, MAGEMin_C needs no per-task process isolation, so no
Python-side multiprocessing is involved in calling this, and (unlike
`multi_path_julia`/`barometry`) no `Distributed`/`pmap` either -- the
parallelism is internal to the one call, across Julia threads, not across
separate worker processes. Needs `PYTHON_JULIACALL_THREADS` (or
`JULIA_NUM_THREADS`) set before `juliacall`'s first import on the Python
side for that internal thread-parallelism to actually have more than one
thread to use -- see `_ensure_julia_ready()` in PetThermoTools/GenFuncs.py.

Algorithm, per point, fully bounded (never an open-ended loop):
1. Bracket establishment: starting from T_start_C, classify (pure liquid+
   fluid, i.e. at/above the liquidus, vs some solid present, i.e. below),
   then step away from the known side (doubling the step each time,
   starting at `bracket_step_C`) until the classification flips, giving a
   (T_lo, T_hi) bracket, T_lo < T_hi. A point that never flips within
   `max_bracket_expansions` steps is marked `converged=false` and dropped --
   never retried indefinitely.
2. Bisection: `T_mid = (T_lo+T_hi)/2` for every still-active point in ONE
   `multi_point_minimization` call per iteration; update brackets; drop a
   point once its own bracket width is below `T_tol_C`. Stops after
   `max_bisection_iters` regardless of remaining width.
3. ONE final `multi_point_minimization` call AT `T_lo` (not `T_hi`) for
   every converged point -- matching alphaMELTS's own `findLiq_MELTS`
   convention of reporting the highest temperature at which solid is
   *already* present (i.e. the first point *just below* the liquidus, not
   the last point confirmed solid-free above it) -- recovers
   `liquidus_phase` (the first solid phase present, explicitly excluding
   "liq"/"fl"), `fluid_saturated`, and the liquid ("liq") phase's own
   `Comp_wt` at that point (very slightly evolved from the bulk, since a
   small amount of solid has already formed there -- this is deliberately
   NOT assumed equal to the bulk, it is read back from MAGEMin_C itself,
   the same way `equilibrate` above extracts any phase's composition via
   `out[k].PP_vec[i].Comp_wt`/`out[k].SS_vec[j].Comp_wt`).

Returns a `Dict{String,Any}` with one entry per FIELD (not one entry per
point -- easier for a Python caller to turn directly into DataFrame
columns): "T_Liq_C" (NaN, not a fake sentinel like 0.0, for anything that
did not converge -- directly avoiding the exact bug already fixed on the
MELTS side of this package), "converged" (Bool), "liquidus_phase" (String,
"" if not converged), "fluid_saturated" (Bool), "T_bracket_width_C"
(Float64, NaN if not converged -- a useful QC signal: a converged-but-wide
bracket means `T_tol_C` was not actually reached within
`max_bisection_iters`), "liquid_oxides" (Vector{String}, the oxide names
`liquid_comp` is ordered by -- shared across all points of one call),
"liquid_comp" (Vector{Vector{Float64}}, one row per point in `liquid_oxides`
order, NaN-filled for a point that did not converge).
"""
function findLiq_julia(; bulk, P_kbar::Vector{Float64}, T_start_C::Vector{Float64},
    Model::String = "ig", fo2_buffer::Union{String,Nothing} = nothing,
    fo2_offset::Union{Vector{Float64},Float64,Nothing} = 0.0,
    suppress::Union{Vector{String},Nothing} = nothing,
    bracket_step_C::Float64 = 25.0, max_bracket_expansions::Int = 8,
    max_bisection_iters::Int = 15, T_tol_C::Float64 = 0.25)

    n = length(P_kbar)
    @assert length(T_start_C) == n "P_kbar and T_start_C must be the same length"

    # --- bulk normalisation: exactly equilibrate's convention, EXCEPT we
    # always take a fresh copy of each row. equilibrate's own "else: new_bulk
    # = bulk" aliases the caller's array; this function then mutates row[9]
    # in place for fo2_buffer, which silently corrupts the caller's bulk
    # vector for any later call (confirmed: a buffered call followed by a
    # plain call that reuses the same bulk object inherits the row[9]=10
    # placeholder and returns garbage -- every phase assemblage comes back
    # wrong, even at absurd temperatures, because the "bulk" is no longer
    # the real composition). Copying costs nothing (10-11 floats per row)
    # and makes this function safe to call repeatedly on the same bulk
    # object -- required for use case (ii), isolated per-step liquidus calls
    # from path()/crystallisation code that likely reuse one evolving bulk
    # vector across many calls. (equilibrate has the identical latent bug;
    # left untouched here, flagged separately.)
    if isa(bulk, Matrix{<:AbstractFloat})
        new_bulk = [copy(bulk[i, :]) for i in 1:size(bulk, 1)]
    else
        new_bulk = [copy(row) for row in bulk]
    end
    @assert length(new_bulk) == n "bulk must have one row/vector per point"

    if fo2_offset === nothing
        fo2_offset = 0.0
    end
    if fo2_buffer !== nothing
        # "O" (ferric-iron proxy, index 9 in the 11-oxide order below) becomes
        # irrelevant once a buffer sets the oxidation state itself, but must
        # still be nonzero -- same placeholder convention as
        # equilibrate/path_main (copied, not reinvented).
        for row in new_bulk
            row[9] = 10
        end
        if !(fo2_offset isa Vector)
            fo2_offset = fill(fo2_offset, n)
        end
    end

    # --- Model name mapping: exactly equilibrate's convention ---
    Model_mapped = Model == "Weller2024" ? "igad" : "ig"

    if Model_mapped == "ig"
        new_bulk_ox = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"]
    else
        new_bulk_ox = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3"]
    end

    rm_list = suppress !== nothing ? remove_phases(suppress, Model_mapped) : nothing

    data = fo2_buffer !== nothing ? Initialize_MAGEMin(Model_mapped, verbose=false, buffer=fo2_buffer) :
                                      Initialize_MAGEMin(Model_mapped, verbose=false)

    LiqFluid = ["liq", "fl"]

    function run_minimization(idx_active::Vector{Int}, T_active::Vector{Float64})
        P_active = P_kbar[idx_active]
        X_active = new_bulk[idx_active]
        kwargs = Dict{Symbol,Any}(:X => X_active, :Xoxides => new_bulk_ox, :sys_in => "wt", :name_solvus => true,
                                   :progressbar => false)
        if fo2_buffer !== nothing
            kwargs[:B] = fo2_offset[idx_active]
        end
        if rm_list !== nothing
            kwargs[:rm_list] = rm_list
        end
        return multi_point_minimization(P_active, T_active, data; kwargs...)
    end

    # Lightweight classification used repeatedly during bracketing/bisection:
    # just above/below, nothing else extracted (cheap -- called many times).
    function classify(idx_active::Vector{Int}, T_active::Vector{Float64})
        out = run_minimization(idx_active, T_active)
        above = falses(length(idx_active))
        for k in eachindex(idx_active)
            PhaseList = copy(out[k].ph)
            if fo2_buffer !== nothing
                filter!(x -> x != fo2_buffer, PhaseList)
            end
            above[k] = length(intersect(LiqFluid, PhaseList)) == length(PhaseList)
        end
        return above
    end

    T_lo = fill(NaN, n)
    T_hi = fill(NaN, n)
    bracketed = falses(n)

    # --- Phase 1: bracket establishment ---
    active = collect(1:n)
    T_current = copy(T_start_C)
    step = fill(bracket_step_C, n)

    above0 = classify(active, T_current[active])
    for (k, idx) in enumerate(active)
        if above0[k]
            T_hi[idx] = T_current[idx]
        else
            T_lo[idx] = T_current[idx]
        end
    end

    for expansion in 1:max_bracket_expansions
        active = [idx for idx in active if !bracketed[idx]]
        isempty(active) && break

        for idx in active
            if !isnan(T_hi[idx]) && isnan(T_lo[idx])
                T_current[idx] = T_hi[idx] - step[idx]   # know "above"; probe down
            elseif !isnan(T_lo[idx]) && isnan(T_hi[idx])
                T_current[idx] = T_lo[idx] + step[idx]    # know "below"; probe up
            end
            step[idx] *= 2
        end

        above = classify(active, T_current[active])
        for (k, idx) in enumerate(active)
            if above[k]
                T_hi[idx] = T_current[idx]
            else
                T_lo[idx] = T_current[idx]
            end
            if !isnan(T_lo[idx]) && !isnan(T_hi[idx])
                bracketed[idx] = true
            end
        end
    end

    # points that never bracketed: give up cleanly, never retried indefinitely
    for idx in 1:n
        if !bracketed[idx]
            T_lo[idx] = NaN
            T_hi[idx] = NaN
        end
    end

    # --- Phase 2: bisection narrowing ---
    active = [idx for idx in 1:n if bracketed[idx]]
    for iter in 1:max_bisection_iters
        active = [idx for idx in active if (T_hi[idx] - T_lo[idx]) >= T_tol_C]
        isempty(active) && break

        T_mid = [(T_lo[idx] + T_hi[idx]) / 2 for idx in active]
        above = classify(active, T_mid)
        for (k, idx) in enumerate(active)
            if above[k]
                T_hi[idx] = T_mid[k]
            else
                T_lo[idx] = T_mid[k]
            end
        end
    end

    # --- Final extraction, AT T_lo (just below the liquidus -- matches
    # findLiq_MELTS's own "first point with solid present" convention) ---
    converged = bracketed  # alias: every bracketed point reached phase 2 and has a usable T_lo
    liquidus_phase = fill("", n)
    fluid_saturated = falses(n)
    liquid_oxides = String[]
    liquid_comp = [fill(NaN, 0) for _ in 1:n]

    final_active = [idx for idx in 1:n if converged[idx]]
    if !isempty(final_active)
        out = run_minimization(final_active, T_lo[final_active])
        for (k, idx) in enumerate(final_active)
            PhaseList = copy(out[k].ph)
            if fo2_buffer !== nothing
                filter!(x -> x != fo2_buffer, PhaseList)
            end
            Type = out[k].ph_type
            Oxides = out[k].oxides
            if isempty(liquid_oxides)
                liquid_oxides = string.(Oxides)
            end

            fluid_saturated[idx] = "fl" in PhaseList
            solids = [p for p in PhaseList if p != "liq" && p != "fl"]
            liquidus_phase[idx] = isempty(solids) ? "" : string(first(solids))

            i, j = 0, 0
            for index in eachindex(PhaseList)
                phase_name = string(PhaseList[index])
                if Type[index] == 0
                    i += 1
                    if phase_name == "liq"
                        liquid_comp[idx] = copy(out[k].PP_vec[i].Comp_wt)
                    end
                else
                    j += 1
                    if phase_name == "liq"
                        liquid_comp[idx] = copy(out[k].SS_vec[j].Comp_wt)
                    end
                end
            end
        end
    end
    if isempty(liquid_oxides)
        liquid_oxides = new_bulk_ox  # no point converged at all -- fall back to the input oxide order
    end
    for idx in 1:n
        if isempty(liquid_comp[idx])
            liquid_comp[idx] = fill(NaN, length(liquid_oxides))
        end
    end

    Finalize_MAGEMin(data)

    T_Liq_C = [converged[idx] ? T_lo[idx] : NaN for idx in 1:n]
    T_bracket_width_C = [converged[idx] ? (T_hi[idx] - T_lo[idx]) : NaN for idx in 1:n]

    return Dict(
        "T_Liq_C" => T_Liq_C,
        "converged" => converged,
        "liquidus_phase" => liquidus_phase,
        "fluid_saturated" => fluid_saturated,
        "T_bracket_width_C" => T_bracket_width_C,
        "liquid_oxides" => liquid_oxides,
        "liquid_comp" => liquid_comp,
    )
end

function path(; comp :: Dict, T_start_C :: Union{Float64, Nothing} = nothing, T_end_C :: Union{Float64, Nothing} = nothing, dt_C :: Union{Float64, Nothing} = nothing, T_C :: Union{Float64, Nothing} = nothing, T_min_C :: Union{Float64, Nothing} = nothing,
            P_start_bar :: Union{Float64, Nothing} = nothing, P_end_bar :: Union{Float64, Nothing} = nothing, dp_bar :: Union{Float64, Nothing} = nothing, P_bar :: Union{Float64, Nothing} = nothing, T_path_C :: Union{Vector{Float64}, Nothing} = nothing,
            P_path_bar :: Union{Vector{Float64}, Nothing} = nothing, frac_xtal :: Union{Float64, Bool, Nothing} = false, 
            phases :: Union{Vector{String}, Nothing} = nothing, Model :: String = "ig",
            fo2_buffer :: Union{String, Nothing} = nothing, fo2_offset :: Union{Float64, Nothing} = 0.0,
            find_liquidus :: Union{Bool, Nothing} = false, suppress :: Union{Vector{String},Nothing} = nothing)

    Results = Dict()

    if frac_xtal === nothing
        frac_xtal = false
    end

    if find_liquidus === nothing
        find_liquidus = false
    end
    
    if Model == "Weller2024"
        model = "igad"
    else
        model = "ig"
    end
    
    if !isnothing(P_bar) && isnothing(P_path_bar)
        P_path_bar = P_bar
    end
    if !isnothing(T_C) && isnothing(T_start_C)
        T_start_C = T_C
    end

    if isnothing(P_path_bar) && isnothing(P_start_bar)
        throw(Exception("Initial P system must be defined"))
    end
    if isnothing(T_path_C) && isnothing(T_start_C) && isnothing(find_liquidus)
        throw(Exception("Starting temperature must be specified or the liquidus must be found"))
    end
    
    if isnothing(comp)
        throw(Exception("No composition specified"))
    end

    if model == "ig"
        bulk = [
            comp["SiO2_Liq"], comp["Al2O3_Liq"], comp["CaO_Liq"], comp["MgO_Liq"], comp["FeOt_Liq"], 
            comp["K2O_Liq"], comp["Na2O_Liq"], comp["TiO2_Liq"], 
            comp["Fe3Fet_Liq"] * (((159.59 / 2) / 71.844) * comp["FeOt_Liq"] - comp["FeOt_Liq"]), 
            comp["Cr2O3_Liq"], comp["H2O_Liq"]
        ]
    else
        bulk = [
            comp["SiO2_Liq"], comp["Al2O3_Liq"], comp["CaO_Liq"], comp["MgO_Liq"], comp["FeOt_Liq"], 
            comp["K2O_Liq"], comp["Na2O_Liq"], comp["TiO2_Liq"], 
            comp["Fe3Fet_Liq"] * (((159.59 / 2) / 71.844) * comp["FeOt_Liq"] - comp["FeOt_Liq"]), 
            comp["Cr2O3_Liq"]
        ]
    end

    if find_liquidus
        if !isnothing(P_path_bar)
            # try
                if typeof(P_path_bar) <: AbstractVector
                    T_Liq = findliq(bulk = bulk, P_kbar = P_path_bar[1] ./ 1000.0, T_start_C = 1400.0,
                                                fo2_buffer = fo2_buffer, fo2_offset = fo2_offset, Model = model, suppress = suppress)
                else
                    T_Liq = findliq(bulk = bulk, P_kbar = P_path_bar / 1000.0, T_start_C = 1400.0,
                                            fo2_buffer = fo2_buffer, fo2_offset = fo2_offset, Model = model, suppress = suppress)
                end
            # catch
            #     return Results
            # end
        elseif !isnothing(P_start_bar)
            # try
                T_Liq = findliq(bulk = bulk, P_kbar = P_start_bar ./ 1000.0, T_start_C = 1400.0,
                        fo2_buffer = fo2_buffer, fo2_offset = fo2_offset, Model = model, suppress = suppress)
            # catch
            #     return Results
            # end
        end
    
        T_start_C = T_Liq
        if !isnothing(T_min_C)
            T_end_C = T_Liq - T_min_C
        end
    end

    if isnothing(T_path_C)
        if isnothing(T_end_C) && isnothing(dt_C)
            T = T_start_C
        elseif !isnothing(T_end_C) && !isnothing(dt_C)
            T = range(T_start_C, stop = T_end_C, length = 1 + round(Int, (T_start_C - T_end_C) / dt_C))
        end
    else
        T = T_path_C
    end
    
    if isnothing(P_path_bar)
        if isnothing(P_end_bar) && isnothing(dp_bar)
            P = P_start_bar
        elseif !isnothing(P_end_bar) && !isnothing(dp_bar)
            P = range(P_start_bar, stop = P_end_bar, length = 1 + round(Int, (P_start_bar - P_end_bar) / dp_bar))
        elseif !isnothing(P_end_bar) && isnothing(dp_bar)
            P = range(P_start_bar, stop = P_end_bar, length = length(collect(T)))
        end
    else
        P = P_path_bar
    end

    if typeof(T) <: AbstractVector && isnothing(P_end_bar) && !isnothing(dp_bar)
        P = range(P_start_bar, stop = P_start_bar - dp_bar * (length(T) - 1), length = length(T))
    elseif typeof(P) <: AbstractVector && isnothing(T_end_C) && !isnothing(dt_C)
        T = range(T_start_C, stop = T_start_C - dt_C * (length(P) - 1), length = length(P))
    elseif typeof(P) <: AbstractVector && isnothing(T_end_C) && isnothing(dt_C) && isnothing(T_path_C)
        T = range(T_start_C, stop = T_start_C, length = length(P))
    end
    
    if typeof(T) <: AbstractVector && typeof(P) <: AbstractVector && length(T) != length(P)
        throw(Exception("Length of P and T vectors are not the same. Check input parameters"))
    end

    if find_liquidus
        if !isnothing(P_path_bar) || !isnothing(T_path_C)
            if typeof(P) <: AbstractVector && typeof(T) <: AbstractVector
                T_Liq_loc = argmin(abs.(T .- T_Liq))
                if T[T_Liq_loc] > T_Liq
                    T = T[T_Liq_loc:end]
                    P = P[T_Liq_loc:end]
                else
                    T = T[T_Liq_loc-1:end]
                    P = P[T_Liq_loc-1:end]
                end
            end
        end
    end

    if typeof(T) <: AbstractVector && !(typeof(P) <: AbstractVector)
        P = fill(P, length(T))
    elseif typeof(P) <: AbstractVector && !(typeof(T) != AbstractVector)
        T = fill(T, length(P))
    end

    T = collect(T)
    P = collect(P)

    Results = path_main(bulk = bulk, T_C = T, P_kbar = P./1000.0, frac_xtal = frac_xtal, phases = phases,
                        Model = model, fo2_buffer = fo2_buffer, fo2_offset = fo2_offset, suppress = suppress)

    return Results
end

function equilibrate(; bulk :: Any, P_kbar :: Vector{Float64}, T_C :: Vector{Float64},
    fo2_buffer :: Union{String, Nothing} = nothing, fo2_offset :: Union{Vector{Float64},Float64, Nothing} = 0.0,
    Model :: String = "ig", suppress :: Union{Vector{String},Nothing} = nothing)
        # Always copy each row rather than aliasing the caller's own bulk
        # array -- same fix as findLiq_julia (see that function's own
        # comment for the full story): without this, the row[9] = 10
        # fo2_buffer placeholder mutation below permanently corrupts the
        # caller's array, silently breaking any later call that reuses it
        # (e.g. a buffered call followed by an unbuffered one on the same
        # bulk vector). This was flagged as a known latent bug here and
        # left unfixed while only findLiq_julia needed it; fixing it now
        # while adding suppress support, since both touch this exact code.
        if isa(bulk, Matrix{<:AbstractFloat})
            new_bulk = [copy(bulk[i, :]) for i in 1:size(bulk, 1)]
        else
            new_bulk = [copy(row) for row in bulk]
        end

        if fo2_buffer !== nothing
            for row in new_bulk
                row[9] = 10
            end
        end

        if fo2_offset === nothing
            fo2_offset = 0.0
        end

        if fo2_buffer !== nothing && eltype(fo2_offset) == Float64
            fo2_offset = fill(fo2_offset, length(P_kbar))
        end

        if Model == "Weller2024"
            Model = "igad"
        else
            Model = "ig"
        end


        if Model === "ig"
        	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];
        else
        	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"];
        end

        # Same pattern as path_main's own suppress handling: compute
        # rm_list once (Model is already the mapped "ig"/"igad" name by
        # this point), then branch at the one minimization call site.
        rm_list = suppress !== nothing ? remove_phases(suppress, Model) : nothing

        if fo2_buffer !== nothing
            data = Initialize_MAGEMin(Model, verbose = false, buffer = fo2_buffer)
            if suppress !== nothing
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset, rm_list = rm_list, name_solvus=true )
            else
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset, name_solvus=true )
            end
        else
            data = Initialize_MAGEMin(Model, verbose = false)
            if suppress !== nothing
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, name_solvus=true )
            else
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", name_solvus=true )
            end
        end

    
    new_bulk_ox_saved =  new_bulk_ox #["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]
    Results = Dict()
    Results["Conditions"] = DataFrame(Symbol("T_C") => T_C, 
                                    Symbol("P_kbar") => P_kbar, 
                                    Symbol("mass_g") => zeros(length(T_C)),
                                    Symbol("mass_per_mole_g/mol") => zeros(length(T_C)),
                                    Symbol("rho_kg/m3") => zeros(length(T_C)),
                                    Symbol("eta_Pa.s") => zeros(length(T_C)),
                                    Symbol("s_kJ/K") => zeros(length(T_C)),
                                    Symbol("h_kJ/mol") => zeros(length(T_C)),
                                    Symbol("alpha_1/K") => zeros(length(T_C)),
                                    Symbol("cp_J/kg/K") => zeros(length(T_C)),
                                    Symbol("log10(fO2)") => zeros(length(T_C)),
                                    Symbol("log10(dQFM)") => zeros(length(T_C)))
    # Results["sys"] = DataFrame(zeros(length(T_C), length(new_bulk_ox_saved)), :auto)
    # rename!(Results["sys"], new_bulk_ox_saved)

    for k in eachindex(T_C)
        Phase = out[k].ph
        if fo2_buffer !== nothing
            filter!(x -> x != fo2_buffer, Phase)
        end 
        Oxides = out[k].oxides
        Type = out[k].ph_type

        Results["Conditions"][k, :] = (T_C[k], P_kbar[k], 100.0,
                                    (out[k].M_sys isa AbstractVector ? out[k].M_sys[1] : out[k].M_sys),
                                    (out[k].rho isa AbstractVector ? out[k].rho[1] : out[k].rho),
                                    (out[k].eta_M isa AbstractVector ? out[k].eta_M[1] : out[k].eta_M),
                                    (out[k].entropy isa AbstractVector ? out[k].entropy[1] : out[k].entropy),
                                    (out[k].enthalpy isa AbstractVector ? out[k].enthalpy[1] : out[k].enthalpy),
                                    (out[k].alpha isa AbstractVector ? out[k].alpha[1] : out[k].alpha),
                                    (out[k].s_cp isa AbstractVector ? out[k].s_cp[1] : out[k].s_cp),
                                    out[k].fO2[1], out[k].dQFM[1])

        if !isempty(Phase)
            phase_counts = Dict{String, Int}()
            i, j = 0, 0
            for index in eachindex(Phase)
                phase_name = string(Phase[index])
                phase_counts[phase_name] = get(phase_counts, phase_name, 0) + 1
                unique_phase_name = string(phase_name, phase_counts[phase_name])

                if !(unique_phase_name in keys(Results))
                    Results[unique_phase_name] = DataFrame([zeros(length(T_C)) for _ in new_bulk_ox], new_bulk_ox)
                    Results[string(unique_phase_name, "_prop")] = DataFrame(Symbol("mass_g")=>zeros(length(T_C)),
                                                                            Symbol("mass%")=>zeros(length(T_C)),
                                                                            Symbol("mol%")=>zeros(length(T_C)),
                                                                            Symbol("vol%")=>zeros(length(T_C)),
                                                                            Symbol("rho_kg/m3")=>zeros(length(T_C)),
                                                                            Symbol("cp_J/kg/K")=>zeros(length(T_C)),
                                                                            Symbol("alpha_1/K")=>zeros(length(T_C)),
                                                                            Symbol("s_kJ/K")=>zeros(length(T_C)),
                                                                            Symbol("h_kJ/mol")=>zeros(length(T_C)))
                end

                # Frac = out[k].ph_frac_wt[index]
                # Results[string(unique_phase_name, "_prop")][k, :mass] = Frac
                if Type[index] == 0
                    i += 1
                    Comp = out[k].PP_vec[i].Comp_wt
                    Results[unique_phase_name][k, Oxides] .= Comp
                    Results[string(unique_phase_name, "_prop")][k,:] = (100.0*out[k].ph_frac_wt[index],
                                                                        out[k].ph_frac_wt[index],
                                                                        out[k].ph_frac[index],
                                                                        out[k].ph_frac_vol[index],
                                                                        out[k].PP_vec[i].rho,
                                                                        out[k].PP_vec[i].cp,
                                                                        out[k].PP_vec[i].alpha,
                                                                        out[k].PP_vec[i].entropy,
                                                                        out[k].PP_vec[i].enthalpy)
                else
                    j += 1
                    Comp = out[k].SS_vec[j].Comp_wt
                    Results[unique_phase_name][k, Oxides] .= Comp
                    Results[string(unique_phase_name, "_prop")][k,:] = (100.0*out[k].ph_frac_wt[index],
                                                                        out[k].ph_frac_wt[index],
                                                                        out[k].ph_frac[index],
                                                                        out[k].ph_frac_vol[index],
                                                                        out[k].SS_vec[j].rho,
                                                                        out[k].SS_vec[j].cp,
                                                                        out[k].SS_vec[j].alpha,
                                                                        out[k].SS_vec[j].entropy,
                                                                        out[k].SS_vec[j].enthalpy)
                end
            end
        end
    end

    Finalize_MAGEMin(data)
    Results_df = Dict(k => pytable(v) for (k, v) in Results)
    Results_df = standardize_mineral_labels(Results_df)
    # println(keys(Results_df))
    return Results_df
end

function findLiq_multi(bulk, P_kbar, T_start_C)
    if isa(bulk, Matrix{<:AbstractFloat})
        new_bulk = [bulk[i, :] for i in 1:size(bulk, 1)]
    else
        new_bulk = bulk
    end
    
    new_bulk = new_bulk/sum(new_bulk);
    if Model === "ig"
    	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];
    else
    	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"];
    end
    
    T_C = T_start_C
    T_C_high = copy(T_C)
    T_C_low = copy(T_C)

    data = fo2_buffer !== nothing ? Initialize_MAGEMin(Model, verbose=false, buffer=fo2_buffer) : Initialize_MAGEMin(Model, verbose=false)

    Liq = ["liq", "fl"]
    Step = collect(1:10)
    
    for j in eachindex(Step)
        if j == 1
            for _ in 1:10
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", name_solvus=true )
                for k in eachindex(T_C)
                    PhaseList = out[k].ph
                    i = intersect(Liq, PhaseList)
                    if length(i) == length(PhaseList)
                        T_C_high[k] = T_C[k]
                        T_C_low[k] = T_C[k] - 50
                    else
                        T_C_low[k] = T_C[k]
                        T_C_high[k] = T_C[k] + 50
                    end            
                end
                T_C .= (T_C_high .+ T_C_low) ./ 2 
            end
        else
            out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", name_solvus=true )
            for k in eachindex(T_C)
                PhaseList = out[k].ph
                i = intersect(Liq, PhaseList)
                if length(i) == length(PhaseList)
                    T_C_high[k] = T_C[k]
                else
                    T_C_low[k] = T_C[k]
                end           
            end
            T_C .= (T_C_high .+ T_C_low) ./ 2  
        end
    end
    
    Finalize_MAGEMin(data)
    return T_C
end

function path_main(; bulk::Vector{Float64}, T_C::Vector{Float64}, P_kbar::Vector{Float64}, 
    frac_xtal::Union{Float64, Bool} = false, phases::Union{Vector{String}, Nothing} = nothing,
    Model::String = "ig", fo2_buffer::Union{String, Nothing} = nothing, 
    fo2_offset::Union{Float64, Nothing} = 0.0, suppress::Union{Vector{String},String,Nothing} = nothing)

    mass = 100.0

    if fo2_buffer !== nothing
        bulk[9] = 10
    end

    bulk_in = bulk;

    if fo2_offset === nothing
        fo2_offset = 0.0
    end

    new_bulk = 100 .* bulk ./sum(bulk);
    if Model === "igad"
        new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"]
    else
        new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]
    end
    Results = Dict()
    Results["sys"] = DataFrame([zeros(length(T_C)) for _ in new_bulk_ox], new_bulk_ox)
    Results["Conditions"] = DataFrame(Symbol("T_C") => T_C, 
                                    Symbol("P_kbar") => P_kbar, 
                                    Symbol("mass_g") => zeros(length(T_C)),
                                    Symbol("mass_per_mole_g/mol") => zeros(length(T_C)),
                                    Symbol("rho_kg/m3") => zeros(length(T_C)),
                                    Symbol("eta_Pa.s") => zeros(length(T_C)),
                                    Symbol("s_kJ/K") => zeros(length(T_C)),
                                    Symbol("h_kJ/mol") => zeros(length(T_C)),
                                    Symbol("alpha_1/K") => zeros(length(T_C)),
                                    Symbol("cp_J/kg/K") => zeros(length(T_C)),
                                    Symbol("log10(fO2)") => zeros(length(T_C)),
                                    Symbol("log10(dQFM)") => zeros(length(T_C)))
    # Results["sys"] = DataFrame([zeros(length(T_C)) for _ in new_bulk_ox], new_bulk_ox)

    if fo2_buffer !== nothing
        data = Initialize_MAGEMin(Model, verbose=false, buffer=fo2_buffer)
    else
        data = Initialize_MAGEMin(Model, verbose=false)
    end
    # data = fo2_buffer !== nothing ? Initialize_MAGEMin(Model, verbose=false, buffer=fo2_buffer) : Initialize_MAGEMin(Model, verbose=false)

    if suppress !== nothing
        rm_list = remove_phases(suppress, Model)
    end

    if fo2_buffer !== nothing
        fo2_offset_vec = fill(fo2_offset, length(T_C))
    else
        fo2_offset_vec = nothing
    end

    if !frac_xtal
        if fo2_buffer !== nothing
            if suppress !== nothing
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset_vec, rm_list = rm_list, name_solvus=true )
            else
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset_vec, name_solvus=true )
            end
        else
            if suppress !== nothing
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, name_solvus=true )
            else
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", name_solvus=true )
            end
        end
        Out_all = deepcopy(out)
    end
    for k in eachindex(T_C)
        if frac_xtal
            if fo2_buffer !== nothing
                if suppress !== nothing
                    out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset, rm_list = rm_list, name_solvus=true )
                else
                    out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset, name_solvus=true )
                end
            else
                if suppress !== nothing
                    out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, name_solvus=true )
                else
                    out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", name_solvus=true )
                end
            end
        else
            out = Out_all[k]
        end

        Phase = out.ph;
        if fo2_buffer !== nothing
            filter!(x -> x != fo2_buffer, Phase)
        end 
        Oxides = out.oxides;
        Type = out.ph_type;

        Results["Conditions"][k, :] = (T_C[k], P_kbar[k], mass, 
                                    (out.M_sys isa AbstractVector ? out.M_sys[1] : out.M_sys),
                                    (out.rho isa AbstractVector ? out.rho[1] : out.rho),
                                    (out.eta_M isa AbstractVector ? out.eta_M[1] : out.eta_M),
                                    (out.entropy isa AbstractVector ? out.entropy[1] : out.entropy),
                                    (out.enthalpy isa AbstractVector ? out.enthalpy[1] : out.enthalpy),
                                    (out.alpha isa AbstractVector ? out.alpha[1] : out.alpha),
                                    (out.s_cp isa AbstractVector ? out.s_cp[1] : out.s_cp),
                                    out.fO2[1], out.dQFM[1])
        Results["sys"][k,Oxides] .= out.bulk_wt
        # Results["sys"][k, Oxides] .= out.bulk

        phase_counts = Dict{String, Int}()
        i, j = 0, 0

        for index in eachindex(Phase)
            phase_name = string(Phase[index])
            phase_counts[phase_name] = get(phase_counts, phase_name, 0) + 1
            unique_phase_name = string(phase_name, phase_counts[phase_name])

            if !(unique_phase_name in keys(Results))
                Results[unique_phase_name] = DataFrame([zeros(length(T_C)) for _ in new_bulk_ox], new_bulk_ox)
                Results[string(unique_phase_name, "_prop")] = DataFrame(Symbol("mass_g")=>zeros(length(T_C)),
                                                                        Symbol("mass%")=>zeros(length(T_C)),
                                                                        Symbol("mol%")=>zeros(length(T_C)),
                                                                        Symbol("vol%")=>zeros(length(T_C)),
                                                                        Symbol("rho_kg/m3")=>zeros(length(T_C)),
                                                                        Symbol("cp_J/kg/K")=>zeros(length(T_C)),
                                                                        Symbol("alpha_1/K")=>zeros(length(T_C)),
                                                                        Symbol("s_kJ/K")=>zeros(length(T_C)),
                                                                        Symbol("h_kJ/mol")=>zeros(length(T_C)))
            end

            # Frac = out.ph_frac_wt[index]
            # Results[string(unique_phase_name, "_prop")][k, :mass] = Frac
            
            if Type[index] == 0
                i = i + 1
                Comp = out.PP_vec[i].Comp_wt;
                Results[unique_phase_name][k, Oxides] .= Comp;
                Results[string(unique_phase_name, "_prop")][k,:] = (mass*out.ph_frac_wt[index],
                                                                    out.ph_frac_wt[index],
                                                                    out.ph_frac[index],
                                                                    out.ph_frac_vol[index],
                                                                    out.PP_vec[i].rho,
                                                                    out.PP_vec[i].cp,
                                                                    out.PP_vec[i].alpha,
                                                                    out.PP_vec[i].entropy,
                                                                    out.PP_vec[i].enthalpy)
            else
                j = j +1
                Comp = out.SS_vec[j].Comp_wt;
                Results[unique_phase_name][k, Oxides] .= Comp;
                Results[string(unique_phase_name, "_prop")][k,:] = (mass*out.ph_frac_wt[index],
                                                                    out.ph_frac_wt[index],
                                                                    out.ph_frac[index],
                                                                    out.ph_frac_vol[index],
                                                                    out.SS_vec[j].rho,
                                                                    out.SS_vec[j].cp,
                                                                    out.SS_vec[j].alpha,
                                                                    out.SS_vec[j].entropy,
                                                                    out.SS_vec[j].enthalpy)
            end
        end

        # if !frac_xtal
        #     bulk = copy(bulk_in)
        #     new_bulk = 100 .* bulk ./ sum(bulk)
        if frac_xtal
            comp = Results["liq1"][k, :]
            if Model == "ig"
                bulk = [comp["SiO2"], comp["Al2O3"], comp["CaO"], comp["MgO"], comp["FeO"], comp["K2O"], comp["Na2O"], comp["TiO2"], comp["O"], comp["Cr2O3"], comp["H2O"]]
            else
                bulk = [comp["SiO2"], comp["Al2O3"], comp["CaO"], comp["MgO"], comp["FeO"], comp["K2O"], comp["Na2O"], comp["TiO2"], comp["O"], comp["Cr2O3"]]
            end
            new_bulk = 100*bulk/sum(bulk)
            if fo2_buffer !== nothing
                new_bulk[9] = 4.0
            end

            mass = mass*Results["liq1_prop"][k,"mass%"]
        end
        if phases !== nothing && all(str in keys(Results) for str in phases)
            break
        end
    end


    Finalize_MAGEMin(data)
    Results_df = Dict(k => pytable(v) for (k, v) in Results)
    Results_df = standardize_mineral_labels(Results_df)
    return Results_df
end