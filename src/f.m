% Sample Functions

function [ out ] = f( t, y )
    out = -.5 * y;
    %Function for ODE45 only, cannot be used for direction field
    %out = zeros(length(y),1);
    %out(1) = y(2);
    %out(2) = -sin(y(1));
end
