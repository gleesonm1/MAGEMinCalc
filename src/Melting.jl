# using MAGEMin_C
# using DataFrames
# using BasicInterpolators
# using PythonCall

using MAGEMin_C
using DataFrames
using BasicInterpolators
using Roots
using PythonCall


function optimize_entropy_4T(T, T_new, s, P, data, bulk, bulk_ox; 
    fo2_buffer :: Union{String, Nothing} = nothing, fo2_offset :: Union{Float64, Nothing} = 0.0, 
    Tp_search :: Bool = false, tol :: Float64 = 1e-4, max_iter :: Int = 10)
"""
optimize_entropy_4T(T, T_new, s, P, data, bulk, bulk_ox; 
                    fo2_buffer::Union{String,Nothing}=nothing, 
                    fo2_offset::Union{Float64,Nothing}=0.0, 
                    Tp_search::Bool=false, 
                    tol::Float64=1e-4, 
                    max_iter::Int=10)

Iteratively determine the temperature corresponding to a target system entropy 
under fixed pressure and bulk composition.

# Description
This routine performs an iterative temperature search (bisection-style) to find 
the equilibrium temperature (`T_ave`) that yields a system entropy (`out_mean.entropy`)
matching the target value `s`. Each iteration calls 
[`single_point_minimization`](@ref) to evaluate equilibrium phase relations at 
different trial temperatures. 

The search begins between two initial temperature bounds (`T` and `T_new`), 
which are adjusted to bracket the desired entropy value. If `Tp_search` is `true`, 
the function first expands the upper temperature bound (`T`) to ensure that 
`out_high.entropy > s`. Once an appropriate bracket is established, the midpoint 
temperature is iteratively refined until convergence within the specified tolerance 
(`tol`) or until the maximum number of iterations (`max_iter`) is reached.

The function optionally supports buffering the system redox state with an 
oxygen fugacity (`fO₂`) buffer, applied as a relative offset (`fo2_offset`) to 
the chosen buffer (`fo2_buffer`).

# Arguments
- `T::Float64`: Initial upper temperature bound (°C).
- `T_new::Float64`: Initial lower temperature bound (°C).
- `s::Union{Float64,Vector{Float64}}`: Target system entropy to match.
- `P::Float64`: Pressure (kbar).
- `data`: Thermodynamic dataset object passed to `single_point_minimization`.
- `bulk::Vector{Float64}`: Bulk composition (wt%) of system components.
- `bulk_ox::Vector{String}`: List of oxide components corresponding to `bulk`.

# Keyword Arguments
- `fo2_buffer::Union{String,Nothing}=nothing`: Name of the oxygen fugacity buffer 
(e.g., `"QFM"`, `"NNO"`, `"IW"`). If `nothing`, no redox buffer is applied.
- `fo2_offset::Union{Float64,Nothing}=0.0`: Log-unit offset from the specified 
`fo2_buffer` (e.g., +1 = one log unit above the buffer).
- `Tp_search::Bool=false`: If `true`, expand the upper temperature bound (`T`) 
until `out_high.entropy > s`.
- `tol::Float64=1e-4`: Convergence tolerance in fractional entropy difference.
- `max_iter::Int=10`: Maximum number of bisection iterations allowed.

# Returns
A tuple:
1. `T_ave::Float64`: Equilibrium temperature (°C) at which modeled entropy matches `s`
within tolerance.
2. `out_mean`: Output structure from `single_point_minimization` corresponding to 
`T_ave`.

"""
    if Tp_search
        if fo2_buffer !== nothing
            out_high = single_point_minimization(P, T, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", B = fo2_offset)
            out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", B = fo2_offset)
        else
            out_high = single_point_minimization(P, T, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
            out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
        end

        if out_high.entropy < s
            while out_high.entropy < s
                T = T + 10
                if fo2_buffer !== nothing
                    out_high = single_point_minimization(P, T, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", B = fo2_offset)
                else
                    out_high = single_point_minimization(P, T, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
                end
            end
        end
    else
        if fo2_buffer !== nothing
            out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", B = fo2_offset)
        else
            out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
        end
    end

    if out_low.entropy > s
        while out_low.entropy > s
            T_new = T_new - 2
            if fo2_buffer !== nothing
                out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", B = fo2_offset)
            else
                out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
            end
        end
    end

    T_ave = (T_new + T)/2
    if fo2_buffer !== nothing
        out_mean = single_point_minimization(P, T_ave, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", B = fo2_offset)
    else
        out_mean = single_point_minimization(P, T_ave, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
    end

    for i in 1:max_iter
        ds = abs((out_mean.entropy isa AbstractVector ? out_mean.entropy[1] : out_mean.entropy) -
        (s isa AbstractVector ? s[1] : s))/(s isa AbstractVector ? s[1] : s)
        if ds < tol
            return T_ave, out_mean
        end

        if out_mean.entropy > s
            T = T_ave
        else
            T_new = T_ave
        end

        T_ave = (T_new + T)/2

        if fo2_buffer !== nothing
            out_mean = single_point_minimization(P, T_ave, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", B = fo2_offset)
        else
            out_mean = single_point_minimization(P, T_ave, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
        end
    end

    return T_ave, out_mean
end


function AdiabaticDecompressionMelting(; comp :: Dict, T_start_C :: Union{Float64, Nothing} = nothing, 
                                    P_start_kbar :: Float64 = 30.0, P_end_kbar :: Float64 = 2.0, 
                                    dp_kbar :: Float64 = 0.2, Model :: String = "ig",
                                    fo2_buffer :: Union{String, Nothing} = nothing, 
                                    fo2_offset :: Union{Float64, Nothing} = 0.0,
                                    Tp_C :: Union{Float64, Nothing} = 1350.0)
    """
        AdiabaticDecompressionMelting(; comp::Dict, 
                                        T_start_C::Union{Float64,Nothing}=nothing,
                                        P_start_kbar::Float64=30.0, 
                                        P_end_kbar::Float64=2.0, 
                                        dp_kbar::Float64=0.2, 
                                        Model::String="ig",
                                        fo2_buffer::Union{String,Nothing}=nothing,
                                        fo2_offset::Union{Float64,Nothing}=0.0,
                                        Tp_C::Union{Float64,Nothing}=1350.0)
        
        Perform an **adiabatic decompression melting simulation** for a specified bulk composition 
        using  MAGEMin. The function tracks changes in temperature, pressure, 
        phase proportions, and compositions during decompression under isentropic conditions.
        
        # Description
        
        This function models the evolution of a bulk silicate composition during adiabatic (isentropic)
        decompression between an initial and final pressure. It uses the MAGEMin equilibrium solver to 
        determine stable phase assemblages and their compositions at each pressure step, assuming either:
        
        - A specified **potential temperature** (`Tp_C`) from which the adiabatic path is optimized, or  
        - A specified **starting temperature** (`T_start_C`) at the initial pressure.
        
        The function can include or exclude an oxygen fugacity buffer, and supports both standard 
        (`"ig"`) and `"Weller2024"` (“igad”) thermodynamic models.
        
        At each pressure increment, MAGEMin is called via `optimize_entropy_4T` (or `single_point_minimization`)
        to compute equilibrium states. The results are accumulated into `DataFrame`s containing phase compositions, 
        mass fractions, and system conditions.
        
        # Arguments
        
        | Argument | Type | Default | Description |
        |-----------|------|----------|-------------|
        | `comp` | `Dict` | — | Bulk composition of the liquid, as oxide weight percentages (e.g., `"SiO2_Liq"`, `"FeOt_Liq"`, `"H2O_Liq"`, etc.). |
        | `T_start_C` | `Union{Float64,Nothing}` | `nothing` | Starting temperature in °C. Required if `Tp_C` is not provided. |
        | `P_start_kbar` | `Float64` | `30.0` | Starting pressure (kbar). |
        | `P_end_kbar` | `Float64` | `2.0` | Ending pressure (kbar). |
        | `dp_kbar` | `Float64` | `0.2` | Pressure decrement between steps (kbar). |
        | `Model` | `String` | `"ig"` | Thermodynamic model. `"Weller2024"` switches to the `"igad"` model. |
        | `fo2_buffer` | `Union{String,Nothing}` | `nothing` | Optional oxygen fugacity buffer (e.g., `"FMQ"`, `"NNO"`). |
        | `fo2_offset` | `Union{Float64,Nothing}` | `0.0` | Offset (in log units) from the specified buffer. Ignored if `fo2_buffer == nothing`. |
        | `Tp_C` | `Union{Float64,Nothing}` | `1350.0` | Mantle potential temperature (°C). If set, the initial entropy is determined by isentropic optimization from this value. |
        
        # Returns
        
        `Dict{String,PyObject}` — A dictionary where each key corresponds to a phase or system component, 
        and each value is a `pandas.DataFrame` (via `pytable`) containing the calculated evolution along the decompression path:
        
        - `"Conditions"` → Columns: `T_C`, `P_kbar`, `s` (entropy)
        - For each phase (e.g., `"ol1"`, `"sp1"`, `"liq1"`):  
            - `<phase>` → Oxide compositions through pressure
            - `<phase>_prop` → Mass fraction through pressure
        
        # Method Overview
        
        1. **Initialization**  
            - Converts the provided `comp` dictionary into a bulk composition vector (`bulk`).
            - Defines oxide names (`new_bulk_ox`) depending on the selected model.
            - Initializes the MAGEMin model (`Initialize_MAGEMin`).
        
        2. **Entropy Optimization / Starting State**  
            - If `Tp_C` is provided, computes the corresponding entropy (`s`) by minimizing the Gibbs energy 
                at low pressure and finding the temperature that reproduces the potential temperature.  
            - Otherwise, starts directly from `T_start_C`.
        
        3. **Adiabatic Decompression Loop**  
            - Iteratively decreases pressure by `dp_kbar` and computes equilibrium phase assemblages using 
                `optimize_entropy_4T`, maintaining constant entropy (`s`).
            - Stores phase compositions, proportions, and system properties.
        
        4. **Output and Cleanup**  
            - Results are collated into a dictionary of `DataFrame`s for easy access.
            - MAGEMin is finalized (`Finalize_MAGEMin`).
        
        # Example
        
        ```julia
        composition = Dict(
            "SiO2_Liq" => 49.0, "Al2O3_Liq" => 15.0, "CaO_Liq" => 10.0,
            "MgO_Liq" => 10.0, "FeOt_Liq" => 8.0, "K2O_Liq" => 1.0,
            "Na2O_Liq" => 2.0, "TiO2_Liq" => 1.0, "Fe3Fet_Liq" => 0.05,
            "Cr2O3_Liq" => 0.2, "H2O_Liq" => 1.0
        )
        
        results = AdiabaticDecompressionMelting(; 
            comp = composition,
            P_start_kbar = 30.0,
            P_end_kbar = 2.0,
            dp_kbar = 0.5,
            Tp_C = 1350.0,
            fo2_buffer = "FMQ",
            fo2_offset = -0.5
        )
        
        results["Conditions"]  # view P–T–S path
        results["liq1_prop"]   # melt mass fraction vs P
    """
    P = collect(range(P_start_kbar, P_end_kbar, 1+round(Int,(P_start_kbar - P_end_kbar)/dp_kbar)));

    if Model == "Weller2024"
        Model = "igad"
    else
        Model = "ig"
    end

    if fo2_offset === nothing
        fo2_offset = 0.0
    end

    if Model == "ig"
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

    if fo2_buffer !== nothing
        bulk[9] = 10
    end

    bulk_in = bulk;

	new_bulk = 100 .* bulk ./sum(bulk);
    if Model === "igad"
        new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"]
    else
        new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]
    end

    Results = Dict()
    
    data = fo2_buffer !== nothing ? Initialize_MAGEMin(Model, verbose=false, buffer=fo2_buffer) : Initialize_MAGEMin(Model, verbose=false)
    if Tp_C !== nothing
        rm_list = remove_phases(["liq"], Model)
        if fo2_buffer !== nothing
            out = single_point_minimization(0.001, Tp_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset, rm_list = rm_list)
        else
            out = single_point_minimization(0.001, Tp_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list)
        end
        s = out.entropy
    end
    
    if Tp_C !== nothing
        if fo2_buffer !== nothing
            T, out = optimize_entropy_4T(Tp_C + 0.8*(P[1]*(100000)/(3300*9.81)), Tp_C, 
                                    s, P[1], data, new_bulk, new_bulk_ox, 
                                    fo2_buffer=fo2_buffer, fo2_offset = fo2_offset, Tp_search=true)
        else
            T, out = optimize_entropy_4T(Tp_C + 0.8*(P[1]*(100000)/(3300*9.81)), Tp_C,
                                    s, P[1], data, new_bulk, new_bulk_ox, Tp_search=true)
        end
    else
        T = T_start_C
        if fo2_buffer !== nothing
            out = single_point_minimization(P[1], T_start_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
        else
            out = single_point_minimization(P[1], T_start_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
        end
        s = out.entropy
    end

    Results["Conditions"] = DataFrame(T_C = zeros(length(P)), P_kbar = zeros(length(P)),s = zeros(length(P)));
    # Results["sys"] = DataFrame([zeros(length(P)) for _ in new_bulk_ox], new_bulk_ox)

    # helper to safely extract first element or return scalar
    getfirst(x) = x isa AbstractVector ? x[1] : x

    for k in eachindex(P)
        bulk = bulk_in;
        new_bulk = bulk/sum(bulk);

        if k > 1
            if fo2_buffer !== nothing
                T, out = optimize_entropy_4T(T, T-3, s, P[k], data, new_bulk, new_bulk_ox, fo2_buffer=fo2_buffer, fo2_offset = fo2_offset)
            else
                T, out = optimize_entropy_4T(T, T-3, s, P[k], data, new_bulk, new_bulk_ox)
            end
        end   

        Phase = out.ph;
        if fo2_buffer !== nothing
            filter!(x -> x != fo2_buffer, Phase)
        end 
        Oxides = out.oxides;
        Type = out.ph_type;

        Results["Conditions"][k, :] = (T, P[k], getfirst(out.entropy))
        # Results["sys"][k, Oxides] .= out.bulk

        phase_counts = Dict{String, Int}()
        i, j = 0, 0

        for index in eachindex(Phase)
            phase_name = string(Phase[index])
            phase_counts[phase_name] = get(phase_counts, phase_name, 0) + 1
            unique_phase_name = string(phase_name, phase_counts[phase_name])

            if !(unique_phase_name in keys(Results))
                Results[unique_phase_name] = DataFrame([zeros(length(P)) for _ in new_bulk_ox], new_bulk_ox)
                Results[string(unique_phase_name, "_prop")] = DataFrame(mass=zeros(length(P)))
            end

            Frac = out.ph_frac_wt[index]
            Results[string(unique_phase_name, "_prop")][k, :mass] = Frac
            
            if Type[index] == 0
                i = i + 1
                Comp = out.PP_vec[i].Comp_wt;
                Results[unique_phase_name][k, Oxides] .= Comp;
            else
                j = j +1
                Comp = out.SS_vec[j].Comp_wt;
                Results[unique_phase_name][k, Oxides] .= Comp;
            end
        end
    end
    Finalize_MAGEMin(data)
    Results_df = Dict(k => pytable(v) for (k, v) in Results)
	return Results_df
end