function [h, t, y] = ODE45( func, tstart, tend, initial_value )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
abstol = repmat([1e-6], length(initial_value), 1)
reltol = repmat([1e-6], length(initial_value), 1)
h_min = eps
h_max = .1 * (tend - tstart)
max_steps = 1e+4  
rk_matrix = [
        [0., 0., 0., 0., 0., 0.];
        [1./4, 0., 0., 0., 0., 0.];
        [3./32, 9./32, 0., 0., 0., 0.];
        [1932./2197, -7200./2197, 7296./2197, 0., 0., 0.];
        [439./216, -8., 3680./513, -845./4104, 0., 0.];
        [-8./27, 2., -3544./2565, 1859./4104, -11./40, 0.]
];
rk_weights = [25./216, 0., 1408./2565, 2197./4104, -1./5, 0.];
rk_weights_tilde = [
            16./135, 
            0., 
            6656./12825, 
            28561./56430, 
            -9./50, 
            2./55
];

rk_nodes = sum(rk_matrix, 2);

%Time measurement initialisation
tic

%Storage
t = tstart;
while t < tend
    if t + h >= tend
        h = tend - t;
    end

end

toc
end
