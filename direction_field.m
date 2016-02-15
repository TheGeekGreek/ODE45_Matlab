x = 0:1:10;
y = x;
[X, Y] = meshgrid(x, y)
V = f(t, Y);
U = ones(size(X));
quiver(X, Y, U, V);
hold on;
tstart = 0;
tend = 10;
N = 10;
y0 = 1;
[t, y] = explicit_euler(@f, y0, tstart, tend, N);
plot(t, y, 'r');
