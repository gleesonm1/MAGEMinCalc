using .Threads
using Distributed

function barometry(args::Dict)
    args_sym = Dict(Symbol(k) => v for (k,v) in args)

    Model              = get(args_sym, :Model, nothing)
    bulk               = get(args_sym, :comp, nothing)
    cores              = get(args_sym, :cores, 1)

    T_initial_C        = get(args_sym, :T_initial_C, 1400.0)
    T_maxdrop_C        = get(args_sym, :T_maxdrop_C, 25.0)
    dt_C               = get(args_sym, :dt_C, 2.0)

    P_bar              = get(args_sym, :P_bar, [1000.0,2000.0])

    fO2_buffer         = get(args_sym, :fO2_buffer, nothing)
    fO2_offset         = get(args_sym, :fO2_offset, nothing)

    results_list = pmap(1:length(P_bar)) do i
        path(comp = bulk, T_start_C = T_initial_C, dt_C = dt_C,  T_min_C = T_maxdrop_C,
            P_bar = P_bar[i], Model = Model,
            fo2_buffer = fO2_buffer, fo2_offset = fO2_offset,
            find_liquidus = true)
    end

    Results = Dict()
    for i in 1:length(P_bar)
        Results["Run $(i-1)"] = results_list[i]
    end

    return Results
end

function multi_path_julia(args::Dict)
    args_sym = Dict(Symbol(k) => v for (k,v) in args)

    Model              = get(args_sym, :Model, nothing)
    bulk               = get(args_sym, :bulk, nothing)
    comp               = get(args_sym, :comp, nothing)

    Frac_solid         = get(args_sym, :Frac_solid, nothing)
    Frac_fluid         = get(args_sym, :Frac_fluid, nothing)

    T_C                = get(args_sym, :T_C, nothing)
    T_path_C           = get(args_sym, :T_path_C, nothing)
    T_start_C          = get(args_sym, :T_start_C, nothing)
    T_end_C            = get(args_sym, :T_end_C, nothing)
    dt_C               = get(args_sym, :dt_C, nothing)
    find_liquidus      = get(args_sym, :find_liquidus,false)

    P_bar              = get(args_sym, :P_bar, nothing)
    P_path_bar         = get(args_sym, :P_path_bar, nothing)
    P_start_bar        = get(args_sym, :P_start_bar, nothing)
    P_end_bar          = get(args_sym, :P_end_bar, nothing)
    dp_bar             = get(args_sym, :dp_bar, nothing)

    fO2_buffer         = get(args_sym, :fO2_buffer, nothing)
    fO2_offset         = get(args_sym, :fO2_offset, nothing)

    isenthalpic        = get(args_sym, :isenthalpic, false)
    isentropic         = get(args_sym, :isentropic, false)
    isochoric          = get(args_sym, :isochoric, false)

    fluid_sat          = get(args_sym, :fluid_sat, false)
    Crystallinity_limit= get(args_sym, :Crystallinity_limit, nothing)
    len                = get(args_sym, :length, 1)

    timeout            = get(args_sym, :timeout, 300)
    suppress           = get(args_sym, :suppress, nothing)

    if len > 1
        # Use pmap (Parallel Map) to run 'path' on different workers
        # This is cleaner for returning a list of results
        results_list = pmap(1:len) do i
            path(
                comp = comp, T_start_C = T_start_C[i], T_end_C = T_end_C[i], dt_C = dt_C[i], T_C = T_C[i],
                P_start_bar = P_start_bar[i], P_end_bar = P_end_bar[i], dp_bar = dp_bar[i], 
                P_bar = P_bar[i], 
                frac_xtal = Frac_solid, Model = Model,
                fo2_buffer = fO2_buffer, fo2_offset = fO2_offset[i],
                find_liquidus = find_liquidus, suppress = suppress
            )
        end

        # Map the list back to your dictionary
        Results = Dict()
        for i in 1:len
            Results["Run $(i-1)"] = results_list[i]
        end
    else
        Results = path(comp = comp, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, T_C = T_C,
                P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, P_bar = P_bar, 
                frac_xtal = Frac_solid, Model =Model,
                fo2_buffer = fO2_buffer, fo2_offset = fO2_offset,
                find_liquidus = find_liquidus, suppress = suppress)
    end
    return Results
end