x = -1:.1:1;
y = x;
[X, Y] = meshgrid(x, y);
lambda = .5;
V = -lambda * Y;
U = ones(size(X));
quiver(X, Y, U, V);
