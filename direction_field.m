function [] = direction_field( func, y0, tstart, tend, bound )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    delta = .1 * bound;
    x = -bound:delta:bound;
    [X, Y] = meshgrid(x, x);
    V = func(X, Y);
    U = ones(size(X));
    normed = 1./sqrt(V.^2 + U.^2);
    quiver(X, Y, normed .* U, normed .* V, 'color', [0, 0, 0], 'ShowArrowHead', 'off', 'Marker', 'o');
    hold on;
    [t, y] = ODE45(func, tstart, tend, y0);   
    plot(t, y, 'color', 'r', 'Marker', 'o');
    xlim([-bound, bound]);
    ylim([-bound, bound]);
    ax = gca;
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.Box = 'off';
end
