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


function BandB_fixingsome1and0(ROPF, BB_parameters, max_time, nb_max_moves, typeofconstraint::T) where T<: Type{MAXkmovesConstraint}
    output_instance_path = ROPF.output_instance_path
    data_path = ROPF.output_decomposition_path
    formulation = ROPF.decomposition
    instance = ROPF.instance_name
    search_strategy = BB_parameters.search_strategy
    branch_strategy = BB_parameters.branch_strategy
    seuil_u = BB_parameters.seuil_u
    seuil_l = BB_parameters.seuil_l
    #initilization
    start_time = time()
    elapsed_time = 0.0
    #SDP
    # dat_file_path = "D:\\repo\\ComplexPowerModels\\test\\ROPFu\\$instance.dat"
    instance_dat_file_path = joinpath(output_instance_path, "$instance.dat")
    SDP_var_list, Bin_var_list, dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO = read_dat_file(instance_dat_file_path)
    blocks_dict, CLIQUE_TREE = read_blocks(data_path, formulation, instance)
    index_var = Dict{String, Int64}()
    nb_bin = length(Bin_var_list)
    for i=1:nb_bin
        var = Bin_var_list[i]
        index_var[var] = i
    end
    output_file = "maxtime$(max_time)_BandBmax$(nb_max_moves)moves_fixingsome1_$(seuil_u)_and0_$(seuil_l)_$(search_strategy)_$(branch_strategy)_$instance.txt"
    isdir("BandB_runs") || mkpath("BandB_runs")
    f = open(joinpath("BandB_runs", output_file), "w")
    write(f, "nb_max_moves = $nb_max_moves \n")
    #B&B
    mosek_pt = readdlm("Mosek_solutions_max$(nb_max_moves)moves\\$instance.dat")
    fixing_from_SDP = -ones(nb_bin)
    shunt_values = [(mosek_pt[i,2], mosek_pt[i,1]) for i in 1:size(mosek_pt,1) if mosek_pt[i,1][1:5]=="Shunt"]
    sort!(shunt_values, rev=true)
    nb_bin_unfixed = nb_bin
    for i in 1:nb_bin
        tuple = shunt_values[i]
        varname = tuple[2]
        if tuple[1] > seuil_u
            fixing_from_SDP[index_var[varname]] = 1
            nb_bin_unfixed -= 1
        elseif tuple[1] < seuil_l
                fixing_from_SDP[index_var[varname]] = 0
                nb_bin_unfixed -= 1
        end
    end
    write(f, "NB binary variables : $nb_bin_unfixed \n")
    println(nb_bin_unfixed)
    best_ub = solve_minlp(ROPF, fixing_from_SDP, index_var, nb_max_moves, typeofconstraint)
    println(best_ub)
    node0 = node(fixing_from_SDP, -Inf)
    node_list = Set([node0])
    open_nodes = Set([])
    mv("Knitro_solution.csv", joinpath("BandB_runs", "BEST_solution_$instance.csv"), force=true)
    write(f, "UB : $best_ub \n")
    write(f, "var=$Bin_var_list \n")
    close(f)
    log = open(joinpath("BandB_runs",output_file[1:end-3]*"csv"), "w")
    write(log, "Nb noeuds explores ; Nb noeuds a explorer ; Best UB ; Best LB ; Gap (%) ; Temps ; Nb noeuds gap optimalite \n")
    close(log)
    nb_explored_nodes = 0
    #exploring nodes
    while length(node_list) > 0 && elapsed_time < max_time
        if search_strategy == "bestfirst"
            next_node = extract_bestlb(node_list)
        elseif search_strategy == "deepfirst"
            next_node = extract_deepfirst(node_list)
        end
        nb_explored_nodes +=1
        println("Node $nb_explored_nodes \n")
        lb_father = next_node.father_lb
        fixing = copy(next_node.fixing)
        delete!(node_list, next_node)
        dict_variables_to_fix = Dict(var => fixing[index_var[var]] for var in Bin_var_list)
        value_bins, value_SDP_var, opt_value, primal_status = construct_SDP(typeofconstraint, instance, blocks_dict, CLIQUE_TREE, SDP_var_list, Bin_var_list,
           dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO, dict_variables_to_fix, nb_max_moves, "")
        if opt_value < (1-10^(-4))*best_ub && (primal_status == MOI.FEASIBLE_POINT || primal_status == MOI.NEARLY_FEASIBLE_POINT)
            if maximum(min(1-b, b) for (name,b) in value_bins) <= 10^(-6)
                #STOP: no nodes after
                #NOTE : vérifier car la relaxation SDP ne donne pas forcément une solution réalisable pour le QCQP non convexe
                ub =  solve_minlp(ROPF, fixing, index_var, nb_max_moves, typeofconstraint)
                f = open(joinpath("BandB_runs", output_file), "a")
                if (ub-opt_value)/ub > 10^(-4) #gap
                    write(f, "OPEN NODE ")
                    push!(open_nodes, open_node(fixing, opt_value))
                end
                write(f, "Node $nb_explored_nodes : fixing=$fixing ; UB=$ub ; LB=$opt_value \n")
                close(f)
                if ub < best_ub
                    best_ub = ub
                    for op_node in node_list
                        if op_node.father_lb >= (1-10^(-4))*best_ub
                            write(f, "CUT : node fixing $(op_node.fixing) \n")
                            delete!(node_list, op_node)
                        end
                    end
                    for op_node in open_nodes
                        if op_node.node_lb >= (1-10^(-4))*best_ub
                            delete!(open_nodes, op_node)
                        end
                    end
                    mv("Knitro_solution.csv", joinpath("BandB_runs", "BEST_solution_$instance.csv"), force=true)
                end
            else
                f = open(joinpath("BandB_runs", output_file), "a")
                write(f, "Node $nb_explored_nodes : fixing=$fixing ;  LB=$opt_value \n")
                close(f)
                if branch_strategy == "0.5"
                    var_to_branch = select_05(value_bins)
                elseif branch_strategy == "1"
                    var_to_branch = select_1(value_bins)
                end
                next_fixing0 = copy(fixing)
                next_fixing0[index_var[var_to_branch]] = 0
                node_left = node(next_fixing0, opt_value)
                push!(node_list, node_left)
                next_fixing1 = copy(fixing)
                next_fixing1[index_var[var_to_branch]] = 1
                node_right = node(next_fixing1, opt_value)
                push!(node_list, node_right)
            end
        else
            f = open(joinpath("BandB_runs", output_file), "a")
            write(f, "Node $nb_explored_nodes : NON REALISABLE ou borne inf>best_ub fixing=$fixing ;  LB=$opt_value \n")
            close(f)
        end

        best_lb = -Inf
        if length(node_list) > 0
            best_lb = (minimum(n.father_lb for n in node_list))
        else
            best_lb = opt_value
        end
        gap = ((best_ub-best_lb)/best_ub*100)
        elapsed_time = time() - start_time
        log = open(joinpath("BandB_runs",output_file[1:end-3]*"csv"), "a")
        write(log, "$nb_explored_nodes ; $(length(node_list)) ; $best_ub ; $best_lb ; $gap ; $elapsed_time ; $(length(open_nodes)) \n")
        close(log)
    end
    f = open(joinpath("BandB_runs", output_file), "a")
    if length(node_list) == 0
        write(f, "B&B DONE \n")
    else
        write(f, "B&B stopped because of time limit, nb nodes not explored : $(length(node_list)) \n")
    end
    write(f, "BEST UB=$(best_ub) \nNb explored nodes : $nb_explored_nodes \n")
    write(f, "Time : $elapsed_time \n")
    write(f, "Nb open nodes : $(length(open_nodes)) \n")
    close(f)
    return best_ub, nb_explored_nodes, open_nodes
end

function BandB_fixingsome1and0(ROPF, BB_parameters, max_time, nb_max_moves, typeofconstraint::T) where T<: Type{NoConstraint}
    output_instance_path = ROPF.output_instance_path
    data_path = ROPF.output_decomposition_path
    formulation = ROPF.decomposition
    instance = ROPF.instance_name
    search_strategy = BB_parameters.search_strategy
    branch_strategy = BB_parameters.branch_strategy
    seuil_u = BB_parameters.seuil_u
    seuil_l = BB_parameters.seuil_l
    #initilization
    start_time = time()
    elapsed_time = 0.0
    #SDP
    # dat_file_path = "D:\\repo\\ComplexPowerModels\\test\\ROPFu\\$instance.dat"
    instance_dat_file_path = joinpath(output_instance_path, "$instance.dat")
    SDP_var_list, Bin_var_list, dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO = read_dat_file(instance_dat_file_path)
    blocks_dict, CLIQUE_TREE = read_blocks(data_path, formulation, instance)
    index_var = Dict{String, Int64}()
    nb_bin = length(Bin_var_list)
    for i=1:nb_bin
        var = Bin_var_list[i]
        index_var[var] = i
    end
    output_file = "maxtime$(max_time)_BandB_fixingsome1_$(seuil_u)_and0_$(seuil_l)_$(search_strategy)_$(branch_strategy)_$instance.txt"
    isdir("BandB_runs") || mkpath("BandB_runs")
    f = open(joinpath("BandB_runs", output_file), "w")
    #B&B
    mosek_pt = readdlm("Mosek_solutions\\$instance.dat")
    fixing_from_SDP = -ones(nb_bin)
    shunt_values = [(mosek_pt[i,2], mosek_pt[i,1]) for i in 1:size(mosek_pt,1) if mosek_pt[i,1][1:5]=="Shunt"]
    sort!(shunt_values, rev=true)
    nb_bin_unfixed = nb_bin
    for i in 1:nb_bin
        tuple = shunt_values[i]
        varname = tuple[2]
        if tuple[1] > seuil_u
            fixing_from_SDP[index_var[varname]] = 1
            nb_bin_unfixed -= 1
        elseif tuple[1] < seuil_l
                fixing_from_SDP[index_var[varname]] = 0
                nb_bin_unfixed -= 1
        end
    end
    write(f, "NB binary variables : $nb_bin_unfixed \n")
    println(nb_bin_unfixed)
    best_ub = solve_minlp(ROPF, fixing_from_SDP, index_var, nb_max_moves, typeofconstraint)
    println(best_ub)
    node0 = node(fixing_from_SDP, -Inf)
    node_list = Set([node0])
    open_nodes = Set([])
    mv("Knitro_solution.csv", joinpath("BandB_runs", "BEST_solution_$instance.csv"), force=true)
    write(f, "UB : $best_ub \n")
    write(f, "var=$Bin_var_list \n")
    close(f)
    log = open(joinpath("BandB_runs",output_file[1:end-3]*"csv"), "w")
    write(log, "Nb noeuds explores ; Nb noeuds a explorer ; Best UB ; Best LB ; Gap (%) ; Temps ; Nb noeuds gap optimalite \n")
    close(log)
    nb_explored_nodes = 0
    #exploring nodes
    while length(node_list) > 0 && elapsed_time < max_time
        if search_strategy == "bestfirst"
            next_node = extract_bestlb(node_list)
        elseif search_strategy == "deepfirst"
            next_node = extract_deepfirst(node_list)
        end
        nb_explored_nodes +=1
        println("Node $nb_explored_nodes \n")
        lb_father = next_node.father_lb
        fixing = copy(next_node.fixing)
        delete!(node_list, next_node)
        dict_variables_to_fix = Dict(var => fixing[index_var[var]] for var in Bin_var_list)
        value_bins, value_SDP_var, opt_value, primal_status = construct_SDP(typeofconstraint, instance, blocks_dict, CLIQUE_TREE, SDP_var_list, Bin_var_list,
           dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO, dict_variables_to_fix, nb_max_moves, "")
        if opt_value < (1-10^(-4))*best_ub && (primal_status == MOI.FEASIBLE_POINT || primal_status == MOI.NEARLY_FEASIBLE_POINT)
            if maximum(min(1-b, b) for (name,b) in value_bins) <= 10^(-6)
                #STOP: no nodes after
                #NOTE : vérifier car la relaxation SDP ne donne pas forcément une solution réalisable pour le QCQP non convexe
                ub =  solve_minlp(ROPF, fixing, index_var, nb_max_moves, typeofconstraint)
                f = open(joinpath("BandB_runs", output_file), "a")
                if (ub-opt_value)/ub > 10^(-4) #gap
                    write(f, "OPEN NODE ")
                    push!(open_nodes, open_node(fixing, opt_value))
                end
                write(f, "Node $nb_explored_nodes : fixing=$fixing ; UB=$ub ; LB=$opt_value \n")
                close(f)
                if ub < best_ub
                    best_ub = ub
                    for op_node in node_list
                        if op_node.father_lb >= (1-10^(-4))*best_ub
                            write(f, "CUT : node fixing $(op_node.fixing) \n")
                            delete!(node_list, op_node)
                        end
                    end
                    for op_node in open_nodes
                        if op_node.node_lb >= (1-10^(-4))*best_ub
                            delete!(open_nodes, op_node)
                        end
                    end
                    mv("Knitro_solution.csv", joinpath("BandB_runs", "BEST_solution_$instance.csv"), force=true)
                end
            else
                f = open(joinpath("BandB_runs", output_file), "a")
                write(f, "Node $nb_explored_nodes : fixing=$fixing ;  LB=$opt_value \n")
                close(f)
                if branch_strategy == "0.5"
                    var_to_branch = select_05(value_bins)
                elseif branch_strategy == "1"
                    var_to_branch = select_1(value_bins)
                end
                next_fixing0 = copy(fixing)
                next_fixing0[index_var[var_to_branch]] = 0
                node_left = node(next_fixing0, opt_value)
                push!(node_list, node_left)
                next_fixing1 = copy(fixing)
                next_fixing1[index_var[var_to_branch]] = 1
                node_right = node(next_fixing1, opt_value)
                push!(node_list, node_right)
            end
        else
            f = open(joinpath("BandB_runs", output_file), "a")
            write(f, "Node $nb_explored_nodes : NON REALISABLE ou borne inf>best_ub fixing=$fixing ;  LB=$opt_value \n")
            close(f)
        end

        best_lb = -Inf
        if length(node_list) > 0
            best_lb = (minimum(n.father_lb for n in node_list))
        else
            best_lb = opt_value
        end
        gap = ((best_ub-best_lb)/best_ub*100)
        elapsed_time = time() - start_time
        log = open(joinpath("BandB_runs",output_file[1:end-3]*"csv"), "a")
        write(log, "$nb_explored_nodes ; $(length(node_list)) ; $best_ub ; $best_lb ; $gap ; $elapsed_time ; $(length(open_nodes)) \n")
        close(log)
    end
    f = open(joinpath("BandB_runs", output_file), "a")
    if length(node_list) == 0
        write(f, "B&B DONE \n")
    else
        write(f, "B&B stopped because of time limit, nb nodes not explored : $(length(node_list)) \n")
    end
    write(f, "BEST UB=$(best_ub) \nNb explored nodes : $nb_explored_nodes \n")
    write(f, "Time : $elapsed_time \n")
    write(f, "Nb open nodes : $(length(open_nodes)) \n")
    close(f)
    return best_ub, nb_explored_nodes, open_nodes
end



function BandB_fixingsome1and0(ROPF, BB_parameters, max_time, nb_max_shunts, typeofconstraint::T) where T<: Type{MAXkshuntsConstraint}
    output_instance_path = ROPF.output_instance_path
    data_path = ROPF.output_decomposition_path
    formulation = ROPF.decomposition
    instance = ROPF.instance_name
    search_strategy = BB_parameters.search_strategy
    branch_strategy = BB_parameters.branch_strategy
    seuil_u = BB_parameters.seuil_u
    seuil_l = BB_parameters.seuil_l
    #initilization
    start_time = time()
    elapsed_time = 0.0
    #SDP
    # dat_file_path = "D:\\repo\\ComplexPowerModels\\test\\ROPFu\\$instance.dat"
    instance_dat_file_path = joinpath(output_instance_path, "$instance.dat")
    SDP_var_list, Bin_var_list, dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO = read_dat_file(instance_dat_file_path)
    blocks_dict, CLIQUE_TREE = read_blocks(data_path, formulation, instance)
    index_var = Dict{String, Int64}()
    nb_bin = length(Bin_var_list)
    for i=1:nb_bin
        var = Bin_var_list[i]
        index_var[var] = i
    end
    output_file = "maxtime$(max_time)_BandBmax$(nb_max_shunts)_fixingsome1_$(seuil_u)_and0_$(seuil_l)_$(search_strategy)_$(branch_strategy)_$instance.txt"
    isdir("BandB_runs") || mkpath("BandB_runs")
    f = open(joinpath("BandB_runs", output_file), "w")
    write(f, "nb_max_shunts = $nb_max_shunts \n")
    #B&B
    if nb_max_shunts == 4
        mosek_pt = readdlm("Mosek_solutions_maxkshunts\\$instance.dat")
    else
        mosek_pt = readdlm("Mosek_solutions_maxkshunts\\max$(nb_max_shunts)_$instance.dat")
    end
    fixing_from_SDP = -ones(nb_bin)
    shunt_values = [(mosek_pt[i,2], mosek_pt[i,1]) for i in 1:size(mosek_pt,1) if mosek_pt[i,1][1:5]=="Shunt"]
    sort!(shunt_values, rev=true)
    nb_bin_unfixed = nb_bin
    for i in 1:nb_bin
        tuple = shunt_values[i]
        varname = tuple[2]
        if tuple[1] > seuil_u
            fixing_from_SDP[index_var[varname]] = 1
            nb_bin_unfixed -= 1
        elseif tuple[1] < seuil_l
                fixing_from_SDP[index_var[varname]] = 0
                nb_bin_unfixed -= 1
        end
    end
    write(f, "NB binary variables : $nb_bin_unfixed \n")
    println(nb_bin_unfixed)
    best_ub = solve_minlp(ROPF, fixing_from_SDP, index_var, nb_max_shunts, typeofconstraint)
    println(best_ub)
    node0 = node(fixing_from_SDP, -Inf)
    node_list = Set([node0])
    open_nodes = Set([])
    mv("Knitro_solution.csv", joinpath("BandB_runs", "BEST_solution_$instance.csv"), force=true)
    write(f, "UB : $best_ub \n")
    write(f, "var=$Bin_var_list \n")
    close(f)
    log = open(joinpath("BandB_runs",output_file[1:end-3]*"csv"), "w")
    write(log, "Nb noeuds explores ; Nb noeuds a explorer ; Best UB ; Best LB ; Gap (%) ; Temps ; Nb noeuds gap optimalite \n")
    close(log)
    nb_explored_nodes = 0
    #exploring nodes
    while length(node_list) > 0 && elapsed_time < max_time
        if search_strategy == "bestfirst"
            next_node = extract_bestlb(node_list)
        elseif search_strategy == "deepfirst"
            next_node = extract_deepfirst(node_list)
        end
        nb_explored_nodes +=1
        println("Node $nb_explored_nodes \n")
        lb_father = next_node.father_lb
        fixing = copy(next_node.fixing)
        delete!(node_list, next_node)
        dict_variables_to_fix = Dict(var => fixing[index_var[var]] for var in Bin_var_list)
        value_bins, value_SDP_var, opt_value, primal_status = construct_SDP(typeofconstraint, instance, blocks_dict, CLIQUE_TREE, SDP_var_list, Bin_var_list,
           dict_quad_ctr, dict_bounds_ctr, dict_constants_ctr, dict_Bin_ctr, dict_MONO, dict_variables_to_fix, nb_max_shunts, "")
        if opt_value < (1-10^(-4))*best_ub && (primal_status == MOI.FEASIBLE_POINT || primal_status == MOI.NEARLY_FEASIBLE_POINT)
            if maximum(min(1-b, b) for (name,b) in value_bins) <= 10^(-6)
                #STOP: no nodes after
                #NOTE : vérifier car la relaxation SDP ne donne pas forcément une solution réalisable pour le QCQP non convexe
                ub =  solve_minlp(ROPF, fixing, index_var, nb_max_shunts, typeofconstraint)
                f = open(joinpath("BandB_runs", output_file), "a")
                if (ub-opt_value)/ub > 10^(-4) #gap
                    write(f, "OPEN NODE ")
                    push!(open_nodes, open_node(fixing, opt_value))
                end
                write(f, "Node $nb_explored_nodes : fixing=$fixing ; UB=$ub ; LB=$opt_value \n")
                close(f)
                if ub < best_ub
                    best_ub = ub
                    for op_node in node_list
                        if op_node.father_lb >= (1-10^(-4))*best_ub
                            write(f, "CUT : node fixing $(op_node.fixing) \n")
                            delete!(node_list, op_node)
                        end
                    end
                    for op_node in open_nodes
                        if op_node.node_lb >= (1-10^(-4))*best_ub
                            delete!(open_nodes, op_node)
                        end
                    end
                    mv("Knitro_solution.csv", joinpath("BandB_runs", "BEST_solution_$instance.csv"), force=true)
                end
            else
                f = open(joinpath("BandB_runs", output_file), "a")
                write(f, "Node $nb_explored_nodes : fixing=$fixing ;  LB=$opt_value \n")
                close(f)
                if branch_strategy == "0.5"
                    var_to_branch = select_05(value_bins)
                elseif branch_strategy == "1"
                    var_to_branch = select_1(value_bins)
                end
                next_fixing0 = copy(fixing)
                next_fixing0[index_var[var_to_branch]] = 0
                node_left = node(next_fixing0, opt_value)
                push!(node_list, node_left)
                next_fixing1 = copy(fixing)
                next_fixing1[index_var[var_to_branch]] = 1
                if length([value for value in next_fixing1 if value==1]) == nb_max_shunts
                    for i in 1:nb_bin
                        if next_fixing1[i] != 1
                            next_fixing1[i] = 0
                        end
                    end
                end
                node_right = node(next_fixing1, opt_value)
                push!(node_list, node_right)
            end
        else
            f = open(joinpath("BandB_runs", output_file), "a")
            write(f, "Node $nb_explored_nodes : NON REALISABLE ou borne inf>best_ub fixing=$fixing ;  LB=$opt_value \n")
            close(f)
        end

        best_lb = -Inf
        if length(node_list) > 0
            best_lb = (minimum(n.father_lb for n in node_list))
        else
            best_lb = opt_value
        end
        gap = ((best_ub-best_lb)/best_ub*100)
        elapsed_time = time() - start_time
        log = open(joinpath("BandB_runs",output_file[1:end-3]*"csv"), "a")
        write(log, "$nb_explored_nodes ; $(length(node_list)) ; $best_ub ; $best_lb ; $gap ; $elapsed_time ; $(length(open_nodes)) \n")
        close(log)
    end
    f = open(joinpath("BandB_runs", output_file), "a")
    if length(node_list) == 0
        write(f, "B&B DONE \n")
    else
        write(f, "B&B stopped because of time limit, nb nodes not explored : $(length(node_list)) \n")
    end
    write(f, "BEST UB=$(best_ub) \nNb explored nodes : $nb_explored_nodes \n")
    write(f, "Time : $elapsed_time \n")
    write(f, "Nb open nodes : $(length(open_nodes)) \n")
    close(f)
    return best_ub, nb_explored_nodes, open_nodes
end
