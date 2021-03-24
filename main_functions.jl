struct ROPF_infos
    instance_name::String
    matpower_instance_path::String
    output_instance_path::String
    decomposition::String
    output_decomposition_path::String
end

struct BB_infos
    search_strategy::String
    branch_strategy::String
    seuil_u::Float64
    seuil_l::Float64
end

abstract type AbstractROPFConstraint end

abstract type NoConstraint <: AbstractROPFConstraint end
abstract type MAXkshuntsConstraint <:AbstractROPFConstraint end
abstract type MAXkmovesConstraint <: AbstractROPFConstraint end

include("solve_SDP.jl")
include("solve_minlp.jl")
include("B&B_fixingsome1and0.jl")


function optvalue_bounds(ROPF::ROPF_infos, nb_max_shunts, typeofconstraint)
    LB, status = solve_SDP(ROPF, typeofconstraint, nb_max_shunts)
    UB = solve_minlp(ROPF, [], Dict{String,Float64}(), nb_max_shunts, typeofconstraint)
    return UB, LB
end


function solve_BandB(ROPF::ROPF_infos, max_time::Int64, BB_parameters::BB_infos, nb_max_shunts, typeofconstraint)
    LB, status = solve_SDP(ROPF, typeofconstraint, nb_max_shunts)
    (UB, nb_nodes, open_nodes) = BandB_fixingsome1and0(ROPF, BB_parameters, max_time, nb_max_shunts, typeofconstraint)
    return UB, LB
end
