using .Threads
using Distributed

function split_vector(v::AbstractVector, n::Int)
    len = length(v)
    size = div(len, n)
    rem = mod(len, n)
    
    indices = Int[]
    curr = 1
    for i in 1:n
        push!(indices, curr)
        curr += size + (i <= rem ? 1 : 0)
    end
    push!(indices, len + 1)
    
    return [v[indices[i]:indices[i+1]-1] for i in 1:n]
end

function phaseDiagram_calc(args::Dict)
    args_sym = Dict(Symbol(k) => v for (k,v) in args)

    Model              = get(args_sym, :Model, nothing)
    comp               = get(args_sym, :comp, nothing)
    cores              = get(args_sym, :cores, 1)

    T_C                = get(args_sym, :T_C, nothing)
    P_bar              = get(args_sym, :P_bar, nothing)

    fO2_buffer         = get(args_sym, :fO2_buffer, nothing)
    fO2_offset         = get(args_sym, :fO2_offset, nothing)

    suppress           = get(args_sym, :suppress, nothing)

    T_new = split_vector(T_C, cores)
    P_new = split_vector(P_bar, cores)

    results_list = pmap(1:cores) do i
        path(
            comp = comp, T_path_C = T_new[i], P_path_bar = P_new[i], 
            Model = Model,
            fo2_buffer = fO2_buffer, fo2_offset = fO2_offset,
            suppress = suppress
        )
    end

    # Map the list back to your dictionary
    Results = Dict()
    for i in 1:cores
        Results["Run $(i-1)"] = results_list[i]
    end

    return Results
end