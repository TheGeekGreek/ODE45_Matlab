function [ t, y ] = explicit_euler( f, y0, tstart, tend, steps )
    h = (tend - tstart)/steps;
    t = tstart:h:tend;
    y = [y0];
    for k = 1:steps
        y(k + 1) = y(k) + h * f(t, y(k));
    end    
end
