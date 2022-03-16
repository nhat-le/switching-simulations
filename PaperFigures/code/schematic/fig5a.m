%%
gammalst = linspace(0, 1.5, 25);
epslst = linspace(0, 0.5, 20);
figure;
[xx,yy] = meshgrid(epslst, gammalst);
plot(xx, yy, 'b.', 'MarkerSize', 10)
hold on
plot(epslst(5), gammalst(9), 'k.', 'MarkerSize', 60);
mymakeaxis('x_label', '\epsilon', 'y_label', '\gamma',...
    'yticks', 0:0.5:1.5, 'font_size', 30)

%%
psw = linspace(0.01, 0.45, 15);
prew = linspace(0.55, 0.99, 10);
figure;
[xx,yy] = meshgrid(prew, psw);
plot(xx, yy, 'b.', 'MarkerSize', 10)
hold on
plot(prew(10), psw(15), 'k.', 'MarkerSize', 60);
ylim([0, 0.45])
mymakeaxis('x_label', 'P_{rew}', 'y_label', 'P_{switch}',...
    'yticks', 0:0.1:0.4, 'font_size', 30)


