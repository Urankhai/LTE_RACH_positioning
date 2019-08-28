close all
% ax + by + c = 0; a^2 + b^2 = 1

a = 1; b = 1; c = 1;
norm_ab = sqrt(a^2 + b^2);
a = a/norm_ab;
b = b/norm_ab;

x = -100:100;
y = -(c + a*x)/b;
asd = [x(1);y(1)]-[x(10);y(10)]; % directional vector

xx = 100*(rand(1)-0.5);
yy = 100*(rand(1)-0.5);
M = [a b; b -a];
intersection = M*[-c; b*xx - a*yy];

v = intersection - [xx;yy];
rfl = v + intersection;
zxc = [xx;yy] - rfl;


figure
hold on
grid on
axis square
plot(x,y)
plot(xx,yy,'ro')
plot(intersection(1),intersection(2),'r*')
plot(rfl(1),rfl(2),'go')

