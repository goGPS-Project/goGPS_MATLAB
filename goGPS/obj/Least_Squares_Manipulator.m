classdef Least_Squares_Manipulator < handle
    properties
        A % Values of the normal matrices [n_obs x n_param_per_epoch]
        A_idx % index of the paramter [n_obs x n_param_per_epoch]
        y % observations  [ n_obs x 1]
        epoch % epoch of the obseravtions and of the A lines [ n_obs x 1]
        param_class % [n_param x 1] each paramter can be part of a class [ 1 : x , 2 : y , 3 : z, 4:, clock, 5 : tropo]
    en
end