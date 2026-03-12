

"""
    standardize_mineral_labels(output_dict, dict_ss)

Renames keys in `output_dict` based on grouping rules defined in `dict_ss`.
Groups similar phases (e.g., cpx and dio) into a single category with sequential numbering.
"""
function standardize_mineral_labels(output_dict::Dict)

    ## phase list taken from MAGEMin_C
    dict_ss = Dict{String, Any}()
    dict_ss["sp"] = "spinel",[ "sp" "spinel";
                                "mt" "magnetite" ] 

    dict_ss["spl"] = "spinel",[ "spl" "spinel";
                                "cm" "chromite";
                                "usp" "uvospinel";
                                "mgt" "magnetite" ] 

    # dict_ss["fsp"] = "feldspar",[   "afs" "alkali-feldspar";
    #                                 "pl" "plagioclase" ] 

    dict_ss["mu"] = "muscovite",[   "pat" "paragonite";
                                    "mu" "muscovite" ] 

    dict_ss["amp"] = "amphibole",[  "gl" "glaucophane";
                                    "act" "actinolite";
                                    "cumm" "cummingtonite";
                                    "tr" "tremolite";
                                    "amp" "amphibole" ] 

    dict_ss["ilm"] = "ilmenite",[   "hem" "hematite";
                                    "ilm" "ilmenite" ] 

    dict_ss["ilmm"] = "ilmenite",[   "hemm" "hematite";
                                    "ilmm" "ilmenite" ] 

    dict_ss["nph"] = "nepheline",[  "K-nph" "K-rich nepheline";
                                    "nph" "nepheline" ] 

    dict_ss["cpx"] = "clinopyroxene",[  "pig" "pigeonite";
                                        "Na-cpx" "Na-rich clinopyroxene";
                                        "cpx" "clinopyroxene" ] 

    dict_ss["dio"] = "diopside",[   "dio" "diopside";
                                    "omph" "omphacite";
                                    "jd" "jadeite" ] 

    dict_ss["occm"] = "carbonate",[ "sid" "siderite";
                                    "ank" "ankerite";
                                    "mag" "magnesite";
                                    "cc" "calcite"] 

    dict_ss["oamp"] = "orthoamhibole",[ "anth" "anthophyllite";
                                        "ged" "gedrite" ]

    # 1. Define the Group Merges as requested
    # Maps specific dict_ss keys to a unified general label
    group_merges = Dict(
        "sp"   => "spl",
        "spl"  => "spl",
        "cpx"  => "cpx",
        "dio"  => "cpx",
        "ilm"  => "ilm",
        "ilmm" => "ilm",
        "amp"  => "amp",
        "oamp" => "amp"
    )

    # 2. Build a Reverse Lookup Map
    # Maps specific sub-labels (e.g., "mgt", "usp", "pig") to the unified general label
    reverse_lookup = Dict{String, String}()
    for (group_key, data) in dict_ss
        target_label = get(group_merges, group_key, group_key)
        
        # Access the matrix of abbreviations (second element of the tuple in dict_ss)
        abbr_matrix = data[2]
        for i in 1:size(abbr_matrix, 1)
            abbr = abbr_matrix[i, 1]
            reverse_lookup[abbr] = target_label
        end
    end

    # 3. Process the output_dict keys
    new_dict = Dict{String, Any}()
    
    # We need to process base names (e.g., "mgt1") and their "_prop" counterparts together
    # Identify unique base keys by stripping numbers and "_prop"
    # Mapping: "original_base" => "new_base"
    rename_map = Dict{String, String}()
    
    # Keep track of counts for each unified group (e.g., "spl" => 1, "spl" => 2)
    group_counters = Dict{String, Int}()

    # Get all unique phase-instance keys (e.g., "pl1", "mgt1", "usp1")
    # We ignore "_prop", "sys", and "Conditions" for the counting phase
    all_keys = collect(keys(output_dict))
    base_instances = unique([replace(k, "_prop" => "") for k in all_keys 
                             if !occursin("_prop", k) && k != "sys" && k != "Conditions"])

    for orig_instance in base_instances
        # Extract the prefix (letters) and suffix (numbers)
        m = match(r"([a-zA-Z]+)([0-9]*)", orig_instance)
        prefix = m.captures[1]
        
        if haskey(reverse_lookup, prefix)
            target_group = reverse_lookup[prefix]
            
            # Increment counter for this group
            group_counters[target_group] = get(group_counters, target_group, 0) + 1
            new_instance = "$(target_group)$(group_counters[target_group])"
            
            rename_map[orig_instance] = new_instance
        else
            # If not in our mapping, keep it as is
            rename_map[orig_instance] = orig_instance
        end
    end

    # 4. Construct the final dictionary
    for (old_key, value) in output_dict
        if old_key == "sys" || old_key == "Conditions"
            new_dict[old_key] = value
            continue
        end

        is_prop = occursin("_prop", old_key)
        base = replace(old_key, "_prop" => "")
        
        if haskey(rename_map, base)
            new_key = is_prop ? "$(rename_map[base])_prop" : rename_map[base]
            new_dict[new_key] = value
        else
            new_dict[old_key] = value
        end
    end

    return new_dict
end