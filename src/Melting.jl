using MAGEMin_C
using Roots

function AdiabaticDecompressionMelting(bulk, T_start_C, P_start_kbar, P_end_kbar, dp_kbar)
    P = P_start_kbar

    gv, z_b, DB, splx_data = init_MAGEMin("ig");
	sys_in = "wt";

	gv.verbose = -1;

	new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");	

    T = T_start_C
    Results = Dict()
    k = 0
    while P > P_end_kbar
        k = k + 1
        println(k)
        if k == 1
            out = point_wise_minimization(P, T, gv, z_b, DB, splx_data, sys_in);
            s = out.entropy
            print(entropy)
        else
            P_kbar = P
            temperature = fzero(T_C -> entropy_equation(T_C, P_kbar, bulk, s), T)
            T = temperature
            out = point_wise_minimization(P, T, gv, z_b, DB, splx_data, sys_in);
        end        
    
        Phase = out.ph;
        Oxides = out.oxides;
        Type = out.ph_type;
        
        Ret = Dict();
        if length(Phase) > 0		
            Ret["sys"] = Dict("Temperature" => T, "Pressure" => P, "Phase" => Phase, "Oxides" => Oxides, "Comp" => out.bulk, "Entropy" => out.entropy);
            
            i = 0
            j = 0
            for index in 1:length(Phase)
                Frac = out.ph_frac_wt[index];
                if Type[index] == 0
                    i = i + 1
                    Comp = out.PP_vec[i].Comp_wt;
                    Ret[Phase[index]] = Dict("Frac" => Frac, "Comp" => Dict(zip(Oxides,Comp)));
                else
                    j = j +1
                    Comp = out.SS_vec[j].Comp_wt;
                    Ret[Phase[index]] = Dict("Frac" => Frac, "Comp" => Dict(zip(Oxides,Comp)));
                end
            end
        end
        Results[k] = Ret
        P = P - dp_kbar
    end

	finalize_MAGEMin(gv, DB);
	return Results
end