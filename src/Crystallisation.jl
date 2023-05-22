using MAGEMin_C
using Roots
# using DataFrames
using Polynomials
using Pandas

function path(bulk, T_C, P_kbar, Frac)
    bulk_in = bulk
    gv, z_b, DB, splx_data = init_MAGEMin("ig");
	sys_in = "wt";

	gv.verbose = -1;

	new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    Results = Dict()
    Results["Conditions"] = DataFrame(columns = ["T_C", "P_kbar"], data = zeros(length(P), 2));
    Results["sys"] = DataFrame(columns = new_bulk_ox, data = zeros(length(P), length(new_bulk_ox)));
    
    for k in eachindex(T_C)
        gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");	
        out = point_wise_minimization(P_kbar[k], T_C[k], gv, z_b, DB, splx_data, sys_in);

        Phase = out.ph;
        Oxides = out.oxides;
        Type = out.ph_type;
        
        iloc(Results["Conditions"])[k] = Dict("T_C" => T_C[k], "P_kbar" => P_kbar[k]);
        iloc(Results["sys"])[k] = Dict(zip(Oxides, out.bulk));
        
        if length(Phase) > 0	
            i = 0
            j = 0
            for index in eachindex(Phase)
                if !(Phase[index] in keys(Results))
                    Results[string(Phase[index])] = DataFrame(columns = new_bulk_ox, data = zeros(length(P), length(new_bulk_ox)));
                    Results[string(Phase[index],"_prop")] = DataFrame(columns = ["Mass"], data = zeros(length(P), 1));
                end

                Frac = out.ph_frac_wt[index];
                iloc(Results[string(Phase[index],"_prop")])[k] = Dict("Mass" => Frac);
                if Type[index] == 0
                    i = i + 1
                    Comp = out.PP_vec[i].Comp_wt;
                    iloc(Results[Phase[index]])[k] = Dict(zip(Oxides,Comp));
                else
                    j = j +1
                    Comp = out.SS_vec[j].Comp_wt;
                    iloc(Results[Phase[index]])[k] = Dict(zip(Oxides,Comp));
                end
            end
        end

        if Frac == 0
            bulk = bulk_in
            new_bulk = bulk/sum(bulk)
        end

        if Frac == 1
            bulk = out.SS_vec["liq"].Comp_wt;
            new_bulk = bulk/sum(bulk)
        end
    end

	finalize_MAGEMin(gv, DB);
	return Results
end