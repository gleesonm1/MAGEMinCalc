using MAGEMin_C
using DataFrames
# using Polynomials
# using Pandas
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

function optimize_entropy(T, s, P, data, bulk, bulk_ox, n)
    T_save = range(T, step=-0.75, length=n)
    P_test = fill(P, n)
    out_save = multi_point_minimization(P_test, T_save, data, X=bulk, Xoxides=bulk_ox, sys_in="wt", progressbar=false)

    # Filter and fit polynomial to estimate T
    s_save = [out.entropy for out in out_save]
    mask = .!ismissing.(s_save)
    coeffs = fit(s_save[mask], T_save[mask], 3)
    return coeffs(s)
end

function AdiabaticDecompressionMelting_new(bulk::Vector{Float64}, T_start_C::Float64, P_start_kbar::Float64, 
                                        P_end_kbar::Float64, dp_kbar::Float64, Frac::Float64)
    # Precompute pressure range
    P = range(P_start_kbar, step=-dp_kbar, stop=P_end_kbar)

    # Normalize bulk composition
    new_bulk = bulk / sum(bulk)
    new_bulk_ox = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"]

    # Initialize system
    data = Initialize_MAGEMin("ig", verbose=false)
    Results = Dict("Conditions" => DataFrame((T_C=zeros(length(P)), P_kbar=zeros(length(P)), s=zeros(length(P)))),
        "sys" => DataFrame(zeros(length(P), length(new_bulk_ox)), columns=new_bulk_ox))

    T = T_start_C
    out = single_point_minimization(P[1], T, data, X=new_bulk, Xoxides=new_bulk_ox, sys_in="wt")
    s = out.entropy

    for (k, pressure) in enumerate(P)
        # Update bulk
        new_bulk = bulk / sum(bulk)

        if k > 1
            # Converge entropy via minimization
            n = 4
            T_next = optimize_entropy(T, s, pressure, data, new_bulk, new_bulk_ox, n)

            # Update temperature for next step
            out = single_point_minimization(pressure, T_next, data, X=new_bulk, Xoxides=new_bulk_ox, sys_in="wt")

            if abs(out.entropy - s)/s > 0.0001
                while abs(out.entropy - s)/s > 0.0001
                    n = n*2
                    T_next, s_check = optimize_entropy(T, s, pressure, data, new_bulk, new_bulk_ox, n)

                    # Update temperature for next step
                    out = single_point_minimization(pressure, T_next, data, X=new_bulk, Xoxides=new_bulk_ox, sys_in="wt")
                    if n > 30
                        break
                    end
                end
            end
            T = T_next
        end

        # Record conditions and system results
        Results["Conditions"][k, :] .= [T, pressure, out.entropy]
        Results["sys"][k, :] .= out.bulk

        phase_counts = Dict{String, Int}()

        # Handle phases
        for (index, phase) in enumerate(out.ph)  # Counter for PhaseList
            phase_name = string(phase)
            phase_counts[phase_name] = get(phase_counts, phase_name, 0) + 1
            unique_phase_name = string(phase_name, phase_counts[phase_name])

            if !(phase in keys(Results))
                Results[unique_phase_name] = DataFrame(zeros(length(P), length(new_bulk_ox)), columns=new_bulk_ox)
                Results["$(unique_phase_name)_prop"] = DataFrame(mass=zeros(length(P)))
            end
            Results["$(unique_phase_name)_prop"][k, :] .= out.ph_frac_wt[index]
            Results[unique_phase_name][k, :] .= out.ph_type[index] == 0 ? out.PP_vec[index].Comp_wt : out.SS_vec[index].Comp_wt
        end
    end

    Finalize_MAGEMin(data)
    return Results
end


function AdiabaticDecompressionMelting(bulk, T_start_C, P_start_kbar, P_end_kbar, dp_kbar, Frac)
    P = collect(range(P_start_kbar, P_end_kbar, 1+round(Int,(P_start_kbar - P_end_kbar)/dp_kbar)));
    bulk_in = bulk

	new_bulk = bulk/sum(bulk);
	new_bulk_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]	

    data = Initialize_MAGEMin("ig", verbose = false)

    T = T_start_C
    Results = Dict()
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
    
        Phase = out.ph;
        Oxides = out.oxides;
        Type = out.ph_type;
        
        iloc(Results["Conditions"])[k] = Dict("T_C" => T, "P_kbar" => P[k], "s" => out.entropy);
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
                    # Results[string(Phase[index])] = create_dataframe(new_bulk_ox, length(P))  
                    # Results[string(Phase[index],"_prop")] = create_dataframe(["Mass"], length(P)) 
                    Results[unique_phase_name] = DataFrame(columns = new_bulk_ox, data = zeros(length(P), length(new_bulk_ox)));
                    Results[string(unique_phase_name,"_prop")] = DataFrame(columns = ["mass"], data = zeros(length(P), 1));
                end

                Fraction = out.ph_frac_wt[index];
                iloc(Results[string(unique_phase_name,"_prop")])[k] = Dict("mass" => Fraction);
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
    end

    Finalize_MAGEMin(data)
	return Results
end