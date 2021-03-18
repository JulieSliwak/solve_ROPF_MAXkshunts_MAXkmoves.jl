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
    seui_l::Float64
end

abstract type AbstractConstraint end

abstract type NoConstraint <: AbstractConstraint end
abstract type MAXkshuntsConstraint <:AbstractConstraint end
abstract type MAXkmovesConstraint <: AbstractConstraint end

include("solve_SDP.jl")
include("solve_minlp.jl")
include("B&B_fixingsome1and0.jl")


function optvalue_bounds(ROPF::ROPF_infos, typeofconstraint::AbstractConstraint)
    LB, status = solve_SDP(ROPF, typeofconstraint)
    UB = solve_minlp(ROPF, [], Dict{String,Float64}(), typeofconstraint)
    return UB, LB
end


function solve2(ROPF::ROPF_infos, max_time::Float64, BB_parameters::BB_infos, typeofconstraint::AbstractConstraint)
    LB, status = solve_SDP(ROPF, typeofconstraint)
    (UB, nb_nodes, open_nodes) = BandB_maxk_fixingsome1and0(ROPF, BB_parameters, max_time, typeofconstraint)
    return UB, LB
end
