using MAGEMin_C
using Roots
using DataFrames
using Polynomials

function create_dataframe(columns, n)
    df = DataFrame()
    for col in columns
        df[!, col] = zeros(n)
    end
    return df
end

function AdiabaticDecompressionMelting(bulk, T_start_C, P_start_kbar, P_end_kbar, dp_kbar)
    P = collect(range(P_start_kbar, P_end_kbar, round(Int,(P_start_kbar - P_end_kbar)/dp_kbar)))

    gv, z_b, DB, splx_data = init_MAGEMin("ig");
	sys_in = "wt";

	gv.verbose = -1;

	new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");	

    T = T_start_C
    Results = Dict()
    out = point_wise_minimization(P[1], T, gv, z_b, DB, splx_data, sys_in);
    s = out.entropy
    k = 0

    Results["Conditions"] = create_dataframe(["T_C", "P_kbar"], length(P))
    Results["sys"] = create_dataframe(new_bulk_ox, length(P))    
    for k in eachindex(P)
        if k > 1
            T_save = zeros(4)
            s_save = zeros(4)
            for i in eachindex(T_save)
                T_save[i] = T - i*2
                out = point_wise_minimization(P[k], T_save[i], gv, z_b, DB, splx_data, sys_in);
                s_save[i] = out.entropy
            end

            coeffs = fit(s_save, T_save, 2)
            T = coeffs(s)

            out = point_wise_minimization(P[k], T, gv, z_b, DB, splx_data, sys_in);
        end        
    
        Phase = out.ph;
        Oxides = out.oxides;
        Type = out.ph_type;
        
        Results["Conditions"][k,:] = Dict("T_C" => T, "P_kbar" => P[k])
        Results["sys"][k,:] = Dict(zip(Oxides, out.bulk))
        
        if length(Phase) > 0	
            i = 0
            j = 0
            for index in eachindex(Phase)
                if !(Phase[index] in keys(Results))
                    Results[Phase[index]] = create_dataframe(new_bulk_ox, length(P))  
                    Results[Phase[index]*"_prop"] = create_dataframe(["Mass"], length(P)) 
                end

                # Frac = out.ph_frac_wt[index];
                # Results[Phase[index]*"_prop"][k,:] = Dict("Mass" => Frac) 
                # if Type[index] == 0
                #     i = i + 1
                #     Comp = out.PP_vec[i].Comp_wt;
                #     Results[Phase[index]][k,:] = Dict("Comp" => Dict(zip(Oxides,Comp)));
                # else
                #     j = j +1
                #     Comp = out.SS_vec[j].Comp_wt;
                #     Results[Phase[index]][k,:] = Dict("Comp" => Dict(zip(Oxides,Comp)));
                # end
            end
        end
    end

	finalize_MAGEMin(gv, DB);
	return Results
end