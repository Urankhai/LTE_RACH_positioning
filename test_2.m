clear all
close all

dist = 0.00001;
t = 100+(0:dist:838*dist)';
% t = 2*pi* 1./lambda;
y = 2*exp(-1j*4*t) + sqrt(2)*exp(-1j*7*t);
x_in = [2, 1j*4, sqrt(2), 1j*7];
% axis([0 2 -0.5 6])
% hold on
figure
hold on
plot(t,real(y),'c','linewidth',3)
plot(t,imag(y),'c','linewidth',3)
title('Data points')
% hold off

% F = @(x, xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata);
% x0 = [2 0.1+4i 1.4 7.2i];

x02 = [3.7i 7.2i];
F2 = @(x,t) fitvector(x,t,y);


% options = optimoptions('lsqcurvefit', 'TolX', 1e-15, 'MaxIter', 1000);
options.TolFun = 10^(-10);
options.MaxFunEvals = 10000;
options.MaxIter = 4000;
% [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,t,y, -inf*ones(size(x0)), inf*ones(size(x0)), options);
[x,resnorm,~,exitflag,output] = lsqcurvefit(F2,x02,t,y, -inf*ones(size(x02)), inf*ones(size(x02)), options);

hold on
plot(t,real(F2(x,t)),'k')
plot(t,imag(F2(x,t)),'k')
hold off
disp(['original = [', num2str(x_in)])
disp(['solution = [', num2str(x)])%,'];', '  diff initial - est = ',num2str(norm(x_in - x))])