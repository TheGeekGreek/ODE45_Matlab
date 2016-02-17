% direction_field:  
%       Draws a normed direction field using an ODE 
%
% parameters: ( function, timespan, y0, x)
%   function -> f(t,y)
%   timespan -> [start_time, end_time]
%   y0       -> starting value
%   x        -> lower_bound:step_size:upper_bound
%
% The vector field will be on a square grid on the interval of x, y0
% defines the starting value of the differential, timespan is the interval
% on which we will solve, and f is a function f(t,y)


function [] = direction_field( func, timespan, y0, x )    
    % create the grid
    [X, Y] = meshgrid(x, x);
    V = func(X, Y);
    U = ones(size(X));
    
    % norm the vectors
    normed = 1./sqrt(V.^2 + U.^2);
    U = normed .* U;
    V = normed .* V;
    
    % graph it
    quiver(X, Y, U, V, 'color', [0, 0, 0], ...
                'ShowArrowHead', 'off', 'Marker', 'o');
    hold on;
    [t, y] = ODE45(func, timespan, y0);   
    plot(t, y, 'color', 'r', 'Marker', 'o');
    xlabel('$t$','Interpreter','LaTex', 'FontSize', 20);
    ylabel('$y(t)$','Interpreter','LaTex', 'FontSize', 20);
    xlim([min(x), max(x)]);
    ylim([min(x), max(x)]);
end
