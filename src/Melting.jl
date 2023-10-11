using MAGEMin_C
# using DataFrames
using Polynomials
using Pandas
using BasicInterpolators

# function create_dataframe(columns, n)
#     df = DataFrame()
#     for col in columns
#         df[!, col] = zeros(n)
#     end
#     return df
# end

# function point_minimisation(P, T, gv, z_b, DB, splx_data, sys_in)
#     out = point_wise_minimization(P, T, gv, z_b, DB, splx_data, sys_in);
#     return out
# end

# function run_minimisation_with_timeout(P, T, gv, z_b, DB, splx_data, sys_in)
#     timeout = 30
#     out = nothing
#     Test = 0
#     elapsed_time = 0.0

#     elapsed_time = @elapsed begin
#         task = @async begin
#             out = point_minimisation(P, T, gv, z_b, DB, splx_data, sys_in)
#         end

#         while !istaskdone(task) && elapsed_time < timeout
#             sleep(0.1)
#             elapsed_time = @elapsed begin
#                 if istaskdone(task)
#                     out = fetch(task)
#                 end
#             end
#         end

#         if !istaskdone(task)
#             # Timeout occurred
#             Test = 1
#         end
#     end

#     return out, Test
# end

function AdiabaticDecompressionMelting(bulk, T_start_C, P_start_kbar, P_end_kbar, dp_kbar, Frac)
    P = collect(range(P_start_kbar, P_end_kbar, 1+round(Int,(P_start_kbar - P_end_kbar)/dp_kbar)));
    bulk_in = bulk
    # gv, z_b, DB, splx_data = init_MAGEMin("ig");
	# sys_in = "wt";

	# gv.verbose = -1;

	new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]

    # gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");	

    data = Initialize_MAGEMin("ig", verbose = false)

    T = T_start_C
    Results = Dict()
    # out = point_wise_minimization(P[1], T, gv, z_b, DB, splx_data, sys_in);#PT_minimisation(P[1], T, bulk); #point_wise_minimization(P[1], T, gv, z_b, DB, splx_data, sys_in);
    out = single_point_minimization(P[1], T, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
    s = out.entropy
    println(s)

    # Results["Conditions"] = create_dataframe(["T_C", "P_kbar"], length(P))
    # Results["sys"] = create_dataframe(new_bulk_ox, length(P))    
    Results["Conditions"] = DataFrame(columns = ["T_C", "P_kbar", "s"], data = zeros(length(P), 3));
    Results["sys"] = DataFrame(columns = new_bulk_ox, data = zeros(length(P), length(new_bulk_ox)));
    for k in eachindex(P)
        bulk = bulk_in;
        new_bulk = bulk/sum(bulk);
        
        print(P[k])
        if k > 1
            # s_check = 0
            # T_next = T-1
            # while abs(s_check - s)/s > 0.0001
            #     gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");
            #     out = point_wise_minimization(P[k], T_next, gv, z_b, DB, splx_data, sys_in); #PT_minimisation(P[k], T_save[i], bulk); #point_wise_minimization(P[k], T_save[i], gv, z_b, DB, splx_data, sys_in);
            #     s_check = out.entropy;
            #     println(s_check)
            #     if s_check < s
            #         T_next = T-abs(T-T_next)/2
            #     else
            #         diff = abs(T-T_next)
            #         T = T_next
            #         T_next = T_next - diff
            #     end
            # end

            n = 4
            T_save = collect(range(T, T - 3.0, n));
            P_test = collect(range(P[k], P[k], n));
            out_save = multi_point_minimization(P_test, T_save, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", progressbar = false);
            
            s_save = zeros(n);
            for i in eachindex(out_save)
                s_save[i] = out_save[i].entropy;
            end

            # Create a Boolean mask to identify non-missing values in s_save
            mask = .!ismissing.(s_save)
            T_save = T_save[mask]
            s_save = s_save[mask]

            # T_save = zeros(3)
            # s_save = zeros(3)
            # for i in eachindex(T_save)
            #     T_save[i] = T - (i-1)*1
            #     # gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");
            #     # out = point_wise_minimization(P[k], T_save[i], gv, z_b, DB, splx_data, sys_in); #PT_minimisation(P[k], T_save[i], bulk); #point_wise_minimization(P[k], T_save[i], gv, z_b, DB, splx_data, sys_in);
            #     out = single_point_minimization(P[k], T_save[i], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
            #     s_save[i] = out.entropy;
            # end
            # print(T_save)

            coeffs = fit(s_save, T_save, 3);
            T_next = coeffs(s);

            # p = CubicSplineInterpolator(reverse(T_save),reverse(s_save),NoBoundaries());
            # T_check = collect(range(maximum(T_save), minimum(T_save), 101));
            # s_check = zeros(101)
            # for i in eachindex(s_check)
            #     s_check[i] = p(T_check[i]);
            # end
            # s_diff = abs.(s_check .- s);
            # T_next = T_check[argmin(s_diff)];
            # T_next = p(s);
            
            
            # print(T_next)
            # gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");
            # out = point_wise_minimization(P[k], T_next, gv, z_b, DB, splx_data, sys_in); #PT_minimisation(P[k], T, bulk); #point_wise_minimization(P[k], T, gv, z_b, DB, splx_data, sys_in);
            out = single_point_minimization(P[k], T_next, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")

            s_check = out.entropy;
            if abs(s_check - s)/s > 0.0001
                while abs(s_check - s)/s > 0.0001
                    n = n*2
                    T_save = collect(range(T, T - 3.0, n));
                    P_test = collect(range(P[k], P[k], n));
                    out_save = multi_point_minimization(P_test, T_save, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt", progressbar = false);
            
                    s_save = zeros(n);
                    for i in eachindex(out_save)
                        s_save[i] = out_save[i].entropy;
                    end
                    
                    # Create a Boolean mask to identify non-missing values in s_save
                    mask = .!ismissing.(s_save)
                    T_save = T_save[mask]
                    s_save = s_save[mask]

                    # T_new = zeros(n)
                    # s_new = zeros(n)
                    # for i in eachindex(T_new)
                    #     T_new[i] = T - (i-1)*(2.5/n)
                    #     # gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");
                    #     # out = point_wise_minimization(P[k], T_new[i], gv, z_b, DB, splx_data, sys_in); #PT_minimisation(P[k], T_save[i], bulk); #point_wise_minimization(P[k], T_save[i], gv, z_b, DB, splx_data, sys_in);
                    #     out = single_point_minimization(P[k], T_new[i], data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
                    #     s_new[i] = out.entropy;
                    # end
                    coeffs = fit(s_save, T_save, 3);
                    T_next = coeffs(s);

                    # p = CubicSplineInterpolator(reverse(T_new),reverse(s_new), NoBoundaries())
                    # T_check = collect(range(maximum(T_new), minimum(T_new), 101));
                    # s_check = zeros(101)
                    # for i in eachindex(s_check)
                    #     s_check[i] = p(T_check[i]);
                    # end
                    # s_diff = abs.(s_check .- s);
                    # T_next = T_check[argmin(s_diff)];
                    
                    # print(T)
                    # gv = define_bulk_rock(gv, new_bulk, new_bulk_ox, sys_in, "ig");
                    # out = point_wise_minimization(P[k], T_next, gv, z_b, DB, splx_data, sys_in); #PT_minimisation(P[k], T, bulk); #point_wise_minimization(P[k], T, gv, z_b, DB, splx_data, sys_in);
                    out = single_point_minimization(P[k], T_next, data, X = new_bulk, Xoxides = new_bulk_ox, sys_in = "wt")
                    s_check = out.entropy;
                    if n > 20
                        T = T_next
                        break
                    end
                end
            end
            T = T_next
            
        end  
        # println(P[k])      
    
        Phase = out.ph;
        Oxides = out.oxides;
        Type = out.ph_type;
        
        iloc(Results["Conditions"])[k] = Dict("T_C" => T, "P_kbar" => P[k], "s" => out.entropy);
        iloc(Results["sys"])[k] = Dict(zip(Oxides, out.bulk));
        
        if length(Phase) > 0	
            i = 0
            j = 0
            for index in eachindex(Phase)
                if !(Phase[index] in keys(Results))
                    # Results[string(Phase[index])] = create_dataframe(new_bulk_ox, length(P))  
                    # Results[string(Phase[index],"_prop")] = create_dataframe(["Mass"], length(P)) 
                    Results[string(Phase[index])] = DataFrame(columns = new_bulk_ox, data = zeros(length(P), length(new_bulk_ox)));
                    Results[string(Phase[index],"_prop")] = DataFrame(columns = ["mass"], data = zeros(length(P), 1));
                end

                Fraction = out.ph_frac_wt[index];
                iloc(Results[string(Phase[index],"_prop")])[k] = Dict("mass" => Fraction);
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

	# finalize_MAGEMin(gv, DB);
    Finalize_MAGEMin(data)
	return Results
end