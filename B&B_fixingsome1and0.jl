include("solve_SDP.jl")
include("solve_minlp.jl")

struct node
    fixing::Array{Int64}
    father_lb::Float64
end

struct open_node
    fixing::Array{Int64}
    node_lb::Float64
end

function extract_deepfirst(node_list)
    max_fixed_var = 0
    selected_node = first(node_list)
    for node in node_list
        nb_fixed_var = 0
        for elem in node.fixing
            if elem == 1 || elem == 0
                nb_fixed_var +=1
            end
        end
        if nb_fixed_var > max_fixed_var
            max_fixed_var = nb_fixed_var
            selected_node = node
        elseif nb_fixed_var == max_fixed_var
            nb_1_node = length([value for value in node.fixing if value==1])
            nb_1_selected_node = length([value for value in selected_node.fixing if value==1])
            if nb_1_node > nb_1_selected_node
                selected_node = node
            end
        end
    end
    return selected_node
end

function extract_bestlb(node_list)
    best_lb = Inf
    selected_node = first(node_list)
    for node in node_list
        lb = node.father_lb
        if lb < best_lb
            selected_node = node
            best_lb = lb
        end
    end
    return selected_node
end


function select_05(value_bins)
    min_distance_to_05 = 1
    var = first(keys(value_bins))
    for (varname, value) in value_bins
        distance_to_05 = abs(value-0.5)
        if distance_to_05 <= min_distance_to_05
            min_distance_to_05 = distance_to_05
            var = varname
        end
    end
    return var
end


function select_1(value_bins)
    u_values = [(value, varname) for (varname,value) in value_bins if abs(1-value) > 10.0^(-6)]
    sort!(u_values, rev=true)
    tuple = u_values[1]
    var = tuple[2]
    return var
end
