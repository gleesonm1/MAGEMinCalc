using MAGEMin_C
using Roots
# using DataFrames
using Polynomials
using Pandas

function create_dataframe(columns, n)
    df = DataFrame()
    for col in columns
        df[!, col] = zeros(n)
    end
    return df
end

function point_minimisation(P, T, gv, z_b, DB, splx_data, sys_in)
    out = point_wise_minimization(P, T, gv, z_b, DB, splx_data, sys_in);
    return out
end

function run_minimisation_with_timeout(P, T, gv, z_b, DB, splx_data, sys_in)
    timeout = 30
    out = nothing
    Test = 0
    elapsed_time = 0.0

    elapsed_time = @elapsed begin
        task = @async begin
            out = point_minimisation(P, T, gv, z_b, DB, splx_data, sys_in)
        end

        while !istaskdone(task) && elapsed_time < timeout
            sleep(0.1)
            elapsed_time = @elapsed begin
                if istaskdone(task)
                    out = fetch(task)
                end
            end
        end

        if !istaskdone(task)
            # Timeout occurred
            Test = 1
        end
    end

    return out, Test
end

function AdiabaticDecompressionMelting(bulk, T_start_C, P_start_kbar, P_end_kbar, dp_kbar)
    P = collect(range(P_start_kbar, P_end_kbar, round(Int,(P_start_kbar - P_end_kbar)/dp_kbar)));
    bulk_in = bulk
    # gv, z_b, DB, splx_data = init_MAGEMin("ig");
	# sys_in = "wt";

	# gv.verbose = -1;

	# new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    # gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");	

    T = T_start_C
    Results = Dict()
    out = PT_minimisation(P[1], T, bulk); #point_wise_minimization(P[1], T, gv, z_b, DB, splx_data, sys_in);
    s = out.entropy

    # Results["Conditions"] = create_dataframe(["T_C", "P_kbar"], length(P))
    # Results["sys"] = create_dataframe(new_bulk_ox, length(P))    
    Results["Conditions"] = DataFrame(columns = ["T_C", "P_kbar"], data = zeros(length(P), 2));
    Results["sys"] = DataFrame(columns = new_bulk_ox, data = zeros(length(P), length(new_bulk_ox)));
    for k in eachindex(P)
        bulk = bulk_in
        if k > 1
            T_save = zeros(3)
            s_save = zeros(3)
            for i in eachindex(T_save)
                T_save[i] = T - (i-1)*0.75
                out = PT_minimisation(P[k], T_save[i], bulk); #point_wise_minimization(P[k], T_save[i], gv, z_b, DB, splx_data, sys_in);
                s_save[i] = out.entropy;
            end
            print(T_save)

            coeffs = fit(s_save, T_save, 2);
            T = coeffs(s);
            
            print(T)
            out = PT_minimisation(P[k], T, bulk); #point_wise_minimization(P[k], T, gv, z_b, DB, splx_data, sys_in);
        end  
        println(P[k])      
    
        Phase = out.ph;
        Oxides = out.oxides;
        Type = out.ph_type;
        
        iloc(Results["Conditions"])[k] = Dict("T_C" => T, "P_kbar" => P[k]);
        iloc(Results["sys"])[k] = Dict(zip(Oxides, out.bulk));
        
        if length(Phase) > 0	
            i = 0
            j = 0
            for index in eachindex(Phase)
                if !(Phase[index] in keys(Results))
                    # Results[string(Phase[index])] = create_dataframe(new_bulk_ox, length(P))  
                    # Results[string(Phase[index],"_prop")] = create_dataframe(["Mass"], length(P)) 
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
    end

	finalize_MAGEMin(gv, DB);
	return Results
end