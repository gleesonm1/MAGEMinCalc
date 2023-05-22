using MAGEMin_C
using Roots
# using DataFrames
using Polynomials
using Pandas

function findliq(bulk, P_kbar, T_start_C)
    bulk_in = bulk;
    gv, z_b, DB, splx_data = init_MAGEMin("ig");
	sys_in = "wt";

	gv.verbose = -1;

	new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    T = T_start_C

    gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");	
    out = point_wise_minimization(P_kbar, T, gv, z_b, DB, splx_data, sys_in);

    Liq = ["liq", "fl"]
    PhaseList = out.ph

    i = intersect(Liq, PhaseList)

    Step = [3, 1, 0.1]
    for k in eachindex(Step)
        if length(i) == length(PhaseList)
            while length(i) == length(PhaseList)
                bulk = bulk_in
                new_bulk = bulk/sum(bulk)
                T = T - Step[k]
                gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");	
                out = point_wise_minimization(P_kbar, T, gv, z_b, DB, splx_data, sys_in);

                PhaseList = out.ph
                i = intersect(Liq, PhaseList)
            end
        end
        if length(i) < length(PhaseList)
            while length(i) < length(PhaseList)
                bulk = bulk_in
                new_bulk = bulk/sum(bulk)
                T = T + Step[k]
                gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");	
                out = point_wise_minimization(P_kbar, T, gv, z_b, DB, splx_data, sys_in);

                PhaseList = out.ph
                i = intersect(Liq, PhaseList)
            end

        end
    end

    finalize_MAGEMin(gv, DB);
    T_Liq = T
    return T_Liq
end

function path(bulk, T_C, P_kbar, Frac)
    Choice = Frac
    bulk_in = bulk
    gv, z_b, DB, splx_data = init_MAGEMin("ig");
	sys_in = "wt";

	gv.verbose = -1;

	new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    Results = Dict()
    Results["Conditions"] = DataFrame(columns = ["T_C", "P_kbar"], data = zeros(length(T_C), 2));
    Results["sys"] = DataFrame(columns = new_bulk_ox, data = zeros(length(T_C), length(new_bulk_ox)));
    
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
                    Results[string(Phase[index])] = DataFrame(columns = new_bulk_ox, data = zeros(length(T_C), length(new_bulk_ox)));
                    Results[string(Phase[index],"_prop")] = DataFrame(columns = ["Mass"], data = zeros(length(T_C), 1));
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

        if Choice == 0
            bulk = bulk_in
            new_bulk = bulk/sum(bulk)
        end

        if Choice == 1
            comp = iloc(Results["liq"])[k]
            bulk = 100*[comp["SiO2"], comp["Al2O3"], comp["CaO"], comp["MgO"], comp["FeO"], comp["K2O"], comp["Na2O"], comp["TiO2"], comp["O"], comp["Cr2O3"], comp["H2O"]]
            new_bulk = bulk/sum(bulk)
        end
    end

	finalize_MAGEMin(gv, DB);
	return Results
end