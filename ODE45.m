% ODE_45:
%       Solves a differential using the Embbedded Runge-Kuta Method
%
% parameters: ( odefun, timespan, bound)
%   odefun -> f(t,y), the right-hand-side of the first order explicit non-stiff 
%             differential equation y'(t) = f(t, y(t)). If the function is vector-
%             valued, f must return a column-vector.
%   timespan -> [start_time, end_time]
%   y0       -> starting value


function [time, y] = ODE45( odefun, timespan, y0)

% Default settings
abstol = repmat([1e-6], length(y0), 1);
reltol = repmat([1e-6], length(y0), 1);
h_min = eps;
h_max = .1 * (timespan(2) - timespan(1));
max_steps = 1e+4;

rk_weights = [
            25./216, 
            0., 
            1408./2565, 
            2197./4104, 
            -1./5, 
            0.
];
rk_weights_tilde = [
            16./135, 
            0., 
            6656./12825, 
            28561./56430, 
            -9./50, 
            2./55
];

% Reshape initial value to columnvector
y0 = reshape(y0, [length(y0), 1]);

t = timespan(1);

% Storage
y = [y0];
time_increments = [];
time = [t];

% Constants
initial_facmax = 2;
successfull_facmax = 5;
unsuccessfull_facmax = 1;
facmin = 0.5;
power = 5;
fac = (.25)^(1./power);
n = length(y0);

% Statistics
successfull_steps = 0;
rejected_steps = 0;

% Time measurement initialisation
tic;

% Calculation of initial stepsize
evaluation1 = odefun(t, y0);
         
scaling = zeros(n, 1);
            
for i = 1:n
    scaling(i) = abstol(i) + abs(y0(i)) * reltol(i);
end
    
d0 = sqrt(1./n * sum((y0./scaling).^2));
d1 = sqrt(1./n * sum((evaluation1./scaling).^2));

if d0 < 1e-5 || d1 < 1e-5
    h0 = 1e-6;
else
    h0 = 1e-2 * (d0/d1);
end
    
y1 = y0 + h0 * evaluation1;
evaluation2 = odefun(y0 + h0, y1);
d2 = sqrt(1./n * sum(((evaluation2 - evaluation1)./scaling).^2))/h0;
            
if max(d1,d2) <= 1e-15
   h1 = max(1e-6, h0 * 1e-3);
else
   h1 = (1e-2/max(d1,d2))^(1./power);
end

h =  min(1e+2 * h0, h1);

% Solving
while t < timespan(2)
    if t + h >= timespan(2)
        h = timespan(2) - t;
    elseif (h > h_max || h < h_min)
        display('Solving might not be successful.');
    elseif length(time_increments) > max_steps
        display('Solving has not been successful.');
        break;
    end 
    
    increments = compute_increments(odefun, y0, t, h);
    
    y1 = y0 + h * (increments * rk_weights);
    y_hat1 = y0 + h * (increments * rk_weights_tilde);
                        
    scaling = zeros(n, 1);
                        
    for i = 1:n
        factor = max(abs(y0(i)) , abs(y1(i)));
        scaling(i) = abstol(i) + (factor * reltol(i));
    end
    
    error = sqrt( 1/n * sum(((y1 - y_hat1)./scaling).^2) );

    if error >= realmin
        r = min(initial_facmax, max(0.1, fac * (1/error)^(1/power)) );
    else
        r = initial_facmax;
    end

    h_new = h * r;                   
                        
    if error <= 1
        % Local extrapolation
        y0 = y_hat1;
                
        t = t+h;
        h = h_new;
        initial_facmax = successfull_facmax;
        y = [y, y0];
        time = [time, t];
    else
        rejected_steps = rejected_steps + 1;
        h = h_new;
        initial_facmax = unsuccessfull_facmax;
    end
  
    successfull_steps = successfull_steps + 1;
    if successfull_steps >= max_steps
        break;
    end
    
end

% Terminating and statistics
disp('ODE45 statistics:')

toc;

disp(sprintf('The number of successfull steps is %d', successfull_steps));
disp(sprintf('The number of rejected steps is %d', rejected_steps));

end

% Increment Computation
function [increments] = compute_increments(func, y, t, h)
    rk_matrix = [
        [0., 0., 0., 0., 0., 0.];
        [1./4, 0., 0., 0., 0., 0.];
        [3./32, 9./32, 0., 0., 0., 0.];
        [1932./2197, -7200./2197, 7296./2197, 0., 0., 0.];
        [439./216, -8., 3680./513, -845./4104, 0., 0.];
        [-8./27, 2., -3544./2565, 1859./4104, -11./40, 0.]
    ];
    
    % Matrix shape
    [m,n] = size(rk_matrix);

    rk_nodes = sum(rk_matrix, 2);
    
    increments = zeros(length(y), m);
    
    for i = 1:m
        increment = increments(:,1:i) * rk_matrix(i,1:i)';
        increments(:,i) = func(t + h * rk_nodes(i), y + h * increment);
    end
end
