include("main_functions.jl")

#example
instance_name = "case1888rte"
matpower_instance_path = "D:\\repo\\data\\data_Matpower\\matpower\\case14.m"
output_instance_path = "D:\\repo\\data\\data_ROPF\\ROPFu"
output_decomposition_path = "D:\\repo\\data\\data_sdp"
max_time = 300 #1 hour
ROPF = ROPF_infos(instance_name,
matpower_instance_path,
output_instance_path,
"cholesky",
output_decomposition_path)


# typeofconstraint = NoConstraint
typeofconstraint = MAXkshuntsConstraint
# typeofconstraint = MAXkmovesConstraint
nb_max_moves_or_shunts = 4


# UB, LB = optvalue_bounds(ROPF, nb_max_moves_or_shunts, typeofconstraint)

BB_parameters = BB_infos("deepfirst", "1", 0.75, 0)

UB, LB = solve_BandB(ROPF, max_time, BB_parameters, nb_max_moves_or_shunts, typeofconstraint)
