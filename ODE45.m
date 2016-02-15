function [h, t, y] = ODE45( func, tstart, tend, initial_value )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
abstol = repmat([1e-6], length(initial_value), 1);
reltol = repmat([1e-6], length(initial_value), 1);
h_min = eps;
h_max = .1 * (tend - tstart);
max_steps = 1e+4;
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

%%%% disgusting ugly code all around WARNING
time_increments = [];
time = [];
facmax = 2;
facmin = 0.5;
power = 5;

rk_nodes = sum(rk_matrix, 2);

%Time measurement initialisation
tic

%Storage
steps = 0;
t = tstart;
while t < tend
    if t + h >= tend
        h = tend - t;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif (h > h_max || h < h_min)
        display('Solving might not be successful.');
    elseif length(time_increments) > max_steps'
        display('Solving has not been successful.');
        break;
    end 
    
    increments = [];%compute_increments(y0, t, h, *self.settings['params']) % todo yannis =)
    
    y1 = y0 + h * dot(increments, rk_weights);
    y_hat1 = y0 + h * dot(increments, rk_weights_tilde);
                        
    scaling = abstol; % unecessary....
                        
    for i = 1:n
        factor = max(abs(y0(i)) , abs(y1(i)));
        scaling(i) = scaling(i) + (factor * reltol(i));

        error = sqrt( 1/n * sum(((y1 - y_hat1)./scaling).^2) );
        
        %r = 0;
        if error >= realmin %%%%
            r = min(facmax, max(0.1, fac * (1/error)^(1/power)) );
        else
            r = facmax;
        end

        h_new = h * r;                   
                        
        if error <= 1
        	time_increments = [time_increments,h];
                
            %Local extrapolation
            y0 = y_hat1;
                
            t = t+h;
            h = h_new;
            facmax = 5; % hardcoded reset
            y = [y,y0];

            time = [time,t];   
        else
            rejected = rejected+1; % todo
            h = h_new;
            facmax = 1;
        end
        
    end
  
    steps = steps + 1;
    if steps >= max_steps
        break;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

end
