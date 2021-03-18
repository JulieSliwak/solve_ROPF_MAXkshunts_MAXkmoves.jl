function solve_minlp(ROPF, fixing, index_var, nb_max_shunts, typeofconstraint::T) where T<: Type{NoConstraint}
    instance = ROPF.instance_name
    root = pwd()
    # date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
    isdir("knitro_runs") || mkpath("knitro_runs")
    if length(fixing) > 0
        outlog = joinpath(pwd(), "knitro_runs", "BB_$(instance).log")
    else
        outlog = joinpath(pwd(), "knitro_runs", "$(instance).log")
    end
    pb_path = ROPF.output_instance_path
    src_ampl_path = joinpath(pwd(), "src_ampl")

    cd(pb_path)

    instance_dat_file = "$(instance).dat"
    cp(instance_dat_file, "minlp_instance.dat", force=true)

    f = open("fixing.dat", "w")
    all_fixing = true
    for (var, index) in index_var
        value = fixing[index]
        if value != -1
            write(f, "$var    $value \n")
        else
            all_fixing = false
        end
    end
    close(f)
    if length(fixing) == 0
        all_fixing = false
    end

    if all_fixing
        open("phase3.run", "w") do f
            println(f, "include $(joinpath(src_ampl_path, "phase3.run"));")
        end
    else
        open("minlp.run", "w") do f
            println(f, "include $(joinpath(src_ampl_path, "minlp.run"));")
        end
    end

    open("minlp.mod", "w") do f
      println(f, "include $(joinpath(src_ampl_path, "minlp.mod"));")
    end

    try
        if all_fixing
            run(`cmd /c ampl phase3.run '>' $(outlog)`)
        else
            run(`cmd /c ampl minlp.run '>' $(outlog)`)
        end
        # run(`cmd /c ampl real_minlp.run `)
    catch
        @warn("AMPL/Knitro failed, returning.")
    end

    mv("knitro_solution.csv", joinpath(root, "knitro_solution.csv"), force=true)
    # mv("knitro_solution.csv", "D:\\repo\\ROPF.jl\\knitro_optimal_solutions\\$(instance).csv", force=true)
    rm("minlp_instance.dat")
    rm("fixing.dat")
    cd(root)
    return read_obj_outlog(outlog)
end

function read_obj_outlog(outlog)
    objective_value = 0.0
    status = ""
    time = 0.0
    lines = readlines(outlog)
    for line in lines
        splitted_line = split(line, ":")
        if splitted_line[1] == "EXIT"
            status = splitted_line[2]
        end
        splitted_line = split(line, "=")
        if splitted_line[1] == "Final objective value               "
            objective_value = parse(Float64,splitted_line[end])
        elseif splitted_line[1] == "Total program time (secs)           "
            t = split(splitted_line[2])[1]
            time += parse(Float64,t)
        end
    end

    if status == " Locally optimal solution found." || status == " Optimal solution found."
        return objective_value
    else
        return +Inf
    end
end


function solve_minlp(ROPF, fixing, index_var, nb_max_shunts, typeofconstraint::T) where T<: Type{MAXkshuntsConstraint}
    instance = ROPF.instance_name
    root = pwd()
    # date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
    isdir("knitro_runs") || mkpath("knitro_runs")
    if length(fixing) > 0
        outlog = joinpath(pwd(), "knitro_runs", "MAX$(nb_max_shunts)shunts_BB_$(instance).log")
    else
        outlog = joinpath(pwd(), "knitro_runs", "MAX$(nb_max_shunts)shunts_$(instance).log")
    end
    pb_path = ROPF.output_instance_path
    src_ampl_path = joinpath(pwd(), "src_ampl")

    cd(pb_path)

    instance_dat_file = "$(instance).dat"
    cp(instance_dat_file, "minlp_instance.dat", force=true)

    f = open("k.dat", "w")
    write(f, "$nb_max_shunts")
    close(f)

    f = open("fixing.dat", "w")
    all_fixing = true
    for (var, index) in index_var
        value = fixing[index]
        if value != -1
            write(f, "$var    $value \n")
        else
            all_fixing = false
        end
    end
    close(f)
    if length(fixing) == 0
        all_fixing = false
    end

    if all_fixing
        open("phase3_maxk.run", "w") do f
            println(f, "include $(joinpath(src_ampl_path, "phase3_maxk.run"));")
        end
    else
        open("minlp_maxk.run", "w") do f
            println(f, "include $(joinpath(src_ampl_path, "minlp_maxk.run"));")
        end
    end

    open("minlp_maxk.mod", "w") do f
      println(f, "include $(joinpath(src_ampl_path, "minlp_maxk.mod"));")
    end

    try
        if all_fixing
            run(`cmd /c ampl phase3_maxk.run '>' $(outlog)`)
        else
            run(`cmd /c ampl minlp_maxk.run '>' $(outlog)`)
        end
        # run(`cmd /c ampl real_minlp.run `)
    catch
        @warn("AMPL/Knitro failed, returning.")
    end

    mv("knitro_solution.csv", joinpath(root, "knitro_solution.csv"), force=true)
    # mv("knitro_solution.csv", "D:\\repo\\ROPF.jl\\knitro_optimal_solutions\\$(instance).csv", force=true)
    rm("minlp_instance.dat")
    rm("fixing.dat")
    cd(root)
    return read_obj_outlog(outlog)
end

function solve_minlp(ROPF, fixing, index_var, nb_max_shunts, typeofconstraint::T) where T<: Type{MAXkmovesConstraint}
    instance = ROPF.instance_name
    root = pwd()
    # date = Dates.format(now(), "yy_u_dd_HH_MM_SS")
    isdir("knitro_runs") || mkpath("knitro_runs")
    if length(fixing) > 0
        outlog = joinpath(pwd(), "knitro_runs", "MAX$(nb_max_shunts)moves_BB_$(instance).log")
    else
        outlog = joinpath(pwd(), "knitro_runs", "MAX$(nb_max_shunts)moves_$(instance).log")
    end
    pb_path = ROPF.output_instance_path
    src_ampl_path = joinpath(pwd(), "src_ampl")

    cd(pb_path)

    instance_dat_file = "$(instance).dat"
    cp(instance_dat_file, "minlp_instance.dat", force=true)

    f = open("k.dat", "w")
    write(f, "$nb_max_shunts")
    close(f)

    cp(joinpath(root,"initial_state\\$(instance).dat"), "init_pt.dat", force=true)


    f = open("fixing.dat", "w")
    all_fixing = true
    for (var, index) in index_var
        value = fixing[index]
        if value != -1
            write(f, "$var    $value \n")
        else
            all_fixing = false
        end
    end
    close(f)
    if length(fixing) == 0
        all_fixing = false
    end

    if all_fixing
        open("phase3_maxkmoves.run", "w") do f
            println(f, "include $(joinpath(src_ampl_path, "phase3_maxkmoves.run"));")
        end
    else
        open("minlp_maxkmoves.run", "w") do f
            println(f, "include $(joinpath(src_ampl_path, "minlp_maxkmoves.run"));")
        end
    end

    open("minlp_maxkmoves.mod", "w") do f
      println(f, "include $(joinpath(src_ampl_path, "minlp_maxkmoves.mod"));")
    end

    try
        if all_fixing
            run(`cmd /c ampl phase3_maxkmoves.run '>' $(outlog)`)
        else
            run(`cmd /c ampl minlp_maxkmoves.run '>' $(outlog)`)
        end
        # run(`cmd /c ampl real_minlp.run `)
    catch
        @warn("AMPL/Knitro failed, returning.")
    end

    mv("knitro_solution.csv", joinpath(root, "knitro_solution.csv"), force=true)
    # mv("knitro_solution.csv", "D:\\repo\\ROPF.jl\\knitro_optimal_solutions\\$(instance).csv", force=true)
    rm("minlp_instance.dat")
    rm("fixing.dat")
    cd(root)
    return read_obj_outlog(outlog)
end
