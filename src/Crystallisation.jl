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
            out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, B = fo2_offset)
        else
            out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
        end
    else
        if suppress !== nothing
            out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list)
        else
            out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
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
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, B = fo2_offset)
                    else
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
                    end
                else
                    if suppress !== nothing
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list)
                    else
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
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
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list, B = fo2_offset)
                    else
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
                    end
                else
                    if suppress !== nothing
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list)
                    else
                        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
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
    fo2_buffer :: Union{String, Nothing} = nothing, fo2_offset :: Union{Float64, Nothing} = 0.0, 
    Model :: String = "ig")
        if isa(bulk, Matrix{<:AbstractFloat})
            new_bulk = [bulk[i, :] for i in 1:size(bulk, 1)]
        else
            new_bulk = bulk
        end
    
        if fo2_buffer !== nothing
            if isa(new_bulk, Matrix{<:AbstractFloat})
                new_bulk[:, 9] .= 10
            elseif isa(new_bulk, Vector{<:AbstractVector})
                for row in new_bulk
                    row[9] = 10
                end
            else
                error("Unsupported format for new_bulk")
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
    
        if fo2_buffer !== nothing
            data = Initialize_MAGEMin(Model, verbose = false, buffer = fo2_buffer)
            out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
        else
            data = Initialize_MAGEMin(Model, verbose = false)
            out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
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
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
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
            out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
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

    data = fo2_buffer !== nothing ? Initialize_MAGEMin(Model, verbose=false, buffer=fo2_buffer) : Initialize_MAGEMin(Model, verbose=false)

    if suppress !== nothing
        rm_list = remove_phases(suppress, Model)
    end

    if !frac_xtal
        if fo2_buffer !== nothing
            if suppress !== nothing
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset, rm_list = rm_list)
            else
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
            end
        else
            if suppress !== nothing
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list)
            else
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
            end
        end
        Out_all = deepcopy(out)
    end
    for k in eachindex(T_C)
        if frac_xtal
            if fo2_buffer !== nothing
                if suppress !== nothing
                    out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset, rm_list = rm_list)
                else
                    out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
                end
            else
                if suppress !== nothing
                    out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", rm_list = rm_list)
                else
                    out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
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
    return Results_df
end