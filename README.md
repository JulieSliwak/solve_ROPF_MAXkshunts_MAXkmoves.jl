# solve_ROPF_MAXkshunts_MAXkmoves.jl

## ROPF instance
Based on a .dat file generated from package ComplexOPF.jl.
Information summarized in a structure :

`struct ROPF_infos
        instance_name::String
        matpower_instance_path::String
        output_instance_path::String
        decomposition::String
        output_decomposition_path::String
end`

Example :
`instance_name = "case1888rte"
matpower_instance_path = "D:\\repo\\data\\data_Matpower\\matpower\\case14.m"
output_instance_path = "D:\\repo\\data\\data_ROPF\\ROPFu"
output_decomposition_path = "D:\\repo\\data\\data_sdp"`

The .dat file can founded at `joinpath(output_instance_path, $(instance_name).dat)`.


## Three optional constraints
* No optional constraint (`typeofconstraint = NoConstraint`)
* Constraint MAXkshunts (`typeofconstraint = MAXkshuntsConstraint`) : at most k shunts activated
* Constraint MAXkmoves (`typeofconstraint = MAXkmovesConstraint`) : at most k moves allowed from an initial state for shunts

## Two types of resolution

### A bounding of the optimal value
With function `UB, LB = optvalue_bounds(ROPF, nb_max_moves_or_shunts, typeofconstraint)`

UB computed with AMPL and KNITRO (3-steps procedure, use of MPEC constraints)

LB computed thanks to an SDP relaxation solved with MOSEK

* Argument `ROPF` has a `ROPF_infos` structure.
* Argument `nb_max_moves_or_shunts` is an integer less or equal than the number of shunts. Useless if `typeofconstraint = NoConstraint`.
* Argument `typeofconstraint` is equal to `NoConstraint`, `MAXkshuntsConstraint` or `MAXkmovesConstraint.`

### A Branch-and-Bound strategy
With function `UB, LB = solve_BandB(ROPF, max_time, BB_parameters, nb_max_moves_or_shunts, typeofconstraint)`

* Argument `ROPF` has a `ROPF_infos` structure.
* Argument `max_time` is a float determining the time limit (in seconds) for the B&B algorithm.
* Argument `BB_parameters` has a `BB_infos` structure with 4 attributes:

`struct BB_infos
    search_strategy::String
    branch_strategy::String
    seuil_u::Float64
    seuil_l::Float64
end`
Example : `BB_parameters = BB_infos("deepfirst", "1", 0.75, 0)`
* Argument `nb_max_moves_or_shunts` is an integer less or equal than the number of shunts. Useless if `typeofconstraint = NoConstraint`.
* Argument `typeofconstraint` is equal to `NoConstraint`, `MAXkshuntsConstraint` or `MAXkmovesConstraint.`
