include("main_functions.jl")

#example
instance_name = "case14"
matpower_instance_path = "D:\\repo\\data\\data_Matpower\\matpower\\case14.m"
output_instance_path = "D:\\repo\\data\\data_ROPF\\ROPFu"
output_decomposition_path = "D:\\repo\\data\\data_sdp"
max_time = 3600 #1 hour
ROPF = ROPF_infos(instance_name,
matpower_instance_path,
output_instance_path,
"cholesky",
output_decomposition_path)


typeofconstraint = MAXkmovesConstraint
nb_max_moves_or_shunts = 4

solve_SDP(ROPF, typeofconstraint, nb_max_moves_or_shunts)
