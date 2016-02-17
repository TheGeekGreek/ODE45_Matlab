% explicit_euler:
%       Solves a differential using the Explicit Euler
%
% parameters: ( function , timespan, y0, bound)
%   function -> f(t,y)
%   timespan -> [start_time, end_time]
%   y0       -> starting value
%   steps    -> the number of steps

function [ t, y ] = explicit_euler( f, timespan, y0, steps )
    h = (timespan(2) - timespan(1))/steps;
    t = timespan(1):h:timespan(2);
    y = [y0];
    
    for k = 1:steps
        y(k + 1) = y(k) + h * f(t, y(k));
    end    
end
