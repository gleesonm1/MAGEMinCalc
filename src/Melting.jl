using MAGEMin_C
using DataFrames
using BasicInterpolators
using PythonCall

# function polyfit_manual(x::Vector{Float64}, y::Vector{Float64}, deg::Int)
#     X = [x .^ i for i in 0:deg]'  # Vandermonde matrix
#     return X \ y  # Least squares fit
# end

# function polyval(c::Vector{Float64}, x::Float64)
#     sum(c[i] * x^(i-1) for i in 1:length(c))
# end

function optimize_entropy(T, s, P, data, bulk, bulk_ox, n)
    T_save = collect(range(T, step=-0.75, length=n))
    P_test = fill(P, n)
    out_save = multi_point_minimization(P_test, T_save, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", progressbar=false)

    # Filter and fit polynomial to estimate T
    s_save = [out.entropy for out in out_save]
    mask = .!ismissing.(s_save)
    coeffs = fit(s_save[mask], T_save[mask], 3)
    return coeffs(s)
end

function optimize_entropy_4T(T, T_new, s, P, data, bulk, bulk_ox; 
    fo2_buffer :: Union{String, Nothing} = nothing, fo2_offset :: Union{Float64, Nothing} = 0.0)
    out_high = single_point_minimization(P, T, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
    out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")

    if out_low.entropy > s
        while out_low.entropy > s
            T_new = T_new - 2
            out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
        end
    end

    T_ave = (T_new + T)/2
    out_mean = single_point_minimization(P, T_ave, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")

    if out_mean.entropy > s
        T = T_ave
    else
        T_new = T_ave
    end

    threshold = 0.0001
    if abs(out_mean.entropy - s)/s > threshold
        count = 0
        while abs(out_mean.entropy - s)/s > threshold
            out_high = single_point_minimization(P, T, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
            out_low = single_point_minimization(P, T_new, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
            
            T_ave = (T_new + T)/2
            out_mean = single_point_minimization(P, T_ave, data, X=bulk, Xoxides=bulk_ox, sys_in="wt")
        
            if out_mean.entropy > s
                T = T_ave
            else
                T_new = T_ave
            end
            
            count = count + 1
            if count > 10
                break
            end
        end
    end

    # Filter and fit polynomial to estimate T
    return T_ave
end


function AdiabaticDecompressionMelting(; comp :: Dict, T_start_C :: Float64, P_start_kbar :: Float64, 
                                    P_end_kbar :: Float64, dp_kbar :: Float64, Model :: String = "ig",
                                    fo2_buffer :: Union{String, Nothing} = nothing, fo2_offset :: Union{Float64, Nothing} = 0.0)

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
    # Results["Conditions"] = DataFrame(T_C = zeros(length(P)), P_kbar = P)
    # Results["sys"] = DataFrame([zeros(length(P)) for _ in ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"])

    data = fo2_buffer !== nothing ? Initialize_MAGEMin(Model, verbose=false, buffer=fo2_buffer) : Initialize_MAGEMin(Model, verbose=false)
    
    T = T_start_C
    # out = single_point_minimization(P[1], T_start_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
    if fo2_buffer !== nothing
        out = single_point_minimization(P[1], T_start_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
    else
        out = single_point_minimization(P[1], T_start_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
    end
    s = out.entropy
    println(s)

    # Results["Conditions"] = create_dataframe(["T_C", "P_kbar"], length(P))
    # Results["sys"] = create_dataframe(new_bulk_ox, length(P))    
    Results["Conditions"] = DataFrame(T_C = zeros(length(P)), P_kbar = zeros(length(P)),s = zeros(length(P)));
    Results["sys"] = DataFrame([zeros(length(P)) for _ in new_bulk_ox], new_bulk_ox)
    for k in eachindex(P)
        bulk = bulk_in;
        new_bulk = bulk/sum(bulk);

        if k > 1
            if fo2_buffer !== nothing
                T = optimize_entropy_4T(T, T-5, s, P[k], data, new_bulk, new_bulk_ox, fo2_buffer=fo2_buffer, fo2_offset = fo2_offset)
            else
                T = optimize_entropy_4T(T, T-5, s, P[k], data, new_bulk, new_bulk_ox)
            end
        end   
        
        # out = single_point_minimization(P[k], T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
        if fo2_buffer !== nothing
            out = single_point_minimization(P[k], T_start_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", B = fo2_offset)
        else
            out = single_point_minimization(P[k], T_start_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
        end

        Phase = out.ph;
        if fo2_buffer !== nothing
            filter!(x -> x != fo2_buffer, Phase)
        end 
        Oxides = out.oxides;
        Type = out.ph_type;

        Results["Conditions"][k, :] = (T, P[k], out.entropy)
        Results["sys"][k, Oxides] .= out.bulk

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