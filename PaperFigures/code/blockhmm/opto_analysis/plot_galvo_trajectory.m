points = 10:-1:1;
nreps = 10;

points = repmat(points, [1 nreps]);

figure;
plot(points, 'k');
hold on
plot(points, 'k.', 'MarkerSize', 30);
mymakeaxis()