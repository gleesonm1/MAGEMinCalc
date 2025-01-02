using MAGEMin_C
# using DataFrames
using Polynomials
using Pandas

function equilibrate(bulk, P_kbar, T_C)
    if isa(bulk, Matrix{<:AbstractFloat})
        new_bulk = [bulk[i, :] for i in 1:size(bulk, 1)]
    else
        new_bulk = bulk
    end

    println(new_bulk)
    
    println(typeof(new_bulk))
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    data = Initialize_MAGEMin("ig", verbose = false)
    out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")

    Results = Dict()
    Results["Conditions"] = DataFrame(columns = ["T_C", "P_kbar"], data = zeros(length(T_C), 2));
    Results["sys"] = DataFrame(columns = new_bulk_ox, data = zeros(length(T_C), length(new_bulk_ox)));

    for k in eachindex(T_C)
        Phase = out[k].ph;
        Oxides = out[k].oxides;
        Type = out[k].ph_type;
        
        iloc(Results["Conditions"])[k] = Dict("T_C" => T_C[k], "P_kbar" => P_kbar[k]);
        iloc(Results["sys"])[k] = Dict(zip(Oxides, out[k].bulk));
        
        if length(Phase) > 0	
            phase_counts = Dict{String, Int}()  # Counter for phases
            i = 0
            j = 0
            for index in eachindex(Phase)
                phase_name = string(Phase[index])
                phase_counts[phase_name] = get(phase_counts, phase_name, 0) + 1
                unique_phase_name = string(phase_name, phase_counts[phase_name])

                if !(unique_phase_name in keys(Results))
                    Results[unique_phase_name] = DataFrame(columns = new_bulk_ox, data = zeros(length(T_C), length(new_bulk_ox)));
                    Results[string(unique_phase_name,"_prop")] = DataFrame(columns = ["mass"], data = zeros(length(T_C), 1));
                end

                Frac = out[k].ph_frac_wt[index];
                iloc(Results[string(unique_phase_name,"_prop")])[k] = Dict("mass" => Frac);
                if Type[index] == 0
                    i = i + 1
                    Comp = out[k].PP_vec[i].Comp_wt;
                    iloc(Results[unique_phase_name])[k] = Dict(zip(Oxides,Comp));
                else
                    j = j +1
                    Comp = out[k].SS_vec[j].Comp_wt;
                    iloc(Results[unique_phase_name])[k] = Dict(zip(Oxides,Comp));
                end
            end
        end
    end

    Finalize_MAGEMin(data)
	return Results

end

function findLiq_multi(bulk, P_kbar, T_start_C)
    if isa(bulk, Matrix{<:AbstractFloat})
        new_bulk = [bulk[i, :] for i in 1:size(bulk, 1)]
    else
        new_bulk = bulk
    end

    println(new_bulk)
    
    println(typeof(new_bulk))
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    T_C = T_start_C
    T_C_high = copy(T_C)
    T_C_low = copy(T_C)

    data = Initialize_MAGEMin("ig", verbose = false)

    Liq = ["liq", "fl"];
    Step = collect(1:10);
    for j in eachindex(Step)
        if j === 1
            for s in eachindex(collect(1:10))
                out = multi_point_minimization(P_kbar, T_C, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
                for k in eachindex(T_C)
                    PhaseList = out[k].ph
                    i = intersect(Liq,PhaseList)
                    if length(i) === length(PhaseList)
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
                i = intersect(Liq,PhaseList)
                if length(i) === length(PhaseList)
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

function findliq(bulk, P_kbar, T_start_C)
    bulk_in = bulk;

	new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];

    T = T_start_C;

    data = Initialize_MAGEMin("ig", verbose = false);
    out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt");

    Liq = ["liq", "fl"];
    PhaseList = out.ph;

    i = intersect(Liq, PhaseList);

    Step = [3, 1, 0.1];
    for k in eachindex(Step)
        if length(i) < length(PhaseList)
            while length(i) < length(PhaseList)
                T = T + Step[k]
                out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")

                PhaseList = out.ph
                i = intersect(Liq, PhaseList)
            end
        end

        if length(i) === length(PhaseList)
            while length(i) === length(PhaseList)
                T = T - Step[k]
                out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")

                PhaseList = out.ph
                i = intersect(Liq, PhaseList)
            end
        end

        T = T + 2*Step[k]+2
        out = single_point_minimization(P_kbar, T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")

        PhaseList = out.ph
        i = intersect(Liq, PhaseList)

        T = T - Step[k]-2
    end

    Finalize_MAGEMin(data)
    T_Liq = T
    return T_Liq
end

function path(bulk, T_C, P_kbar, Frac, phases)
    Choice = Frac;
    bulk_in = bulk;

	new_bulk = 100*bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    Results = Dict()
    Results["Conditions"] = DataFrame(columns = ["T_C", "P_kbar"], data = zeros(length(T_C), 2));
    Results["sys"] = DataFrame(columns = new_bulk_ox, data = zeros(length(T_C), length(new_bulk_ox)));
    
    data = Initialize_MAGEMin("ig", verbose = false)

    for k in eachindex(T_C)
        out = single_point_minimization(P_kbar[k], T_C[k], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
    
        Phase = out.ph;
        Oxides = out.oxides;
        Type = out.ph_type;
        
        iloc(Results["Conditions"])[k] = Dict("T_C" => T_C[k], "P_kbar" => P_kbar[k]);
        iloc(Results["sys"])[k] = Dict(zip(Oxides, out.bulk));
        
        if length(Phase) > 0
            phase_counts = Dict{String, Int}()  # Counter for phases

            i = 0
            j = 0
            for index in eachindex(Phase)
                phase_name = string(Phase[index])
                phase_counts[phase_name] = get(phase_counts, phase_name, 0) + 1
                unique_phase_name = string(phase_name, phase_counts[phase_name])

                if !(unique_phase_name in keys(Results))
                    Results[unique_phase_name] = DataFrame(columns = new_bulk_ox, data = zeros(length(T_C), length(new_bulk_ox)));
                    Results[string(unique_phase_name,"_prop")] = DataFrame(columns = ["mass"], data = zeros(length(T_C), 1));
                end

                Frac = out.ph_frac_wt[index];
                iloc(Results[string(unique_phase_name,"_prop")])[k] = Dict("mass" => Frac);
                if Type[index] == 0
                    i = i + 1
                    Comp = out.PP_vec[i].Comp_wt;
                    iloc(Results[unique_phase_name])[k] = Dict(zip(Oxides,Comp));
                else
                    j = j +1
                    Comp = out.SS_vec[j].Comp_wt;
                    iloc(Results[unique_phase_name])[k] = Dict(zip(Oxides,Comp));
                end
            end
        end

        if Choice == 0
            bulk = bulk_in;
            new_bulk = 100*bulk/sum(bulk);
        end

        if Choice == 1
            comp = iloc(Results["liq1"])[k]
            bulk = [comp["SiO2"], comp["Al2O3"], comp["CaO"], comp["MgO"], comp["FeO"], comp["K2O"], comp["Na2O"], comp["TiO2"], comp["O"], comp["Cr2O3"], comp["H2O"]]
            new_bulk = 100*bulk/sum(bulk)
            # new_bulk = round.(new_bulk, digits = 3)
            # if new_bulk[10] < 0.01
            #     new_bulk[10] = 0.01
            # end
            # print(new_bulk)
        end

        if phases !== 0
            found = 0
            for str in keys(Results)
                if str in phases
                    found = found + 1
                end
            end
            
            if found === length(phases)
                break
            end
        end
    end

    Finalize_MAGEMin(data)
	return Results
end