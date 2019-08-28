clear all
close all

L = 3;
c = 3e8;
t = [20, 23, 35]/c;
a = [1, .7, .3]';
rot = [0, pi/6, pi/3]';

att_rot_coeff = a(1:L);%.*exp(-1j*phi(1:L));

k = (1:839)';
f = 1250;

A_org = [att_rot_coeff, 2*pi * f * t(1:L)'/(2*pi*f/c)];

NUFT = exp(-1j*k * (2*pi * f * t(1:L)));
h = NUFT * att_rot_coeff + 0*(randn(839,1) + 1j*randn(839,1));

figure
hold on
plot(k,abs(h),'g','linewidth', 3)


F = @(x, xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata) + x(5)*exp(-x(6)*xdata);
% x0 = [0.9, 0.0005, 0.5, 0.0006, 0.5, 0.0009];
x0 = [1, 1j*21/c*1250*2*pi, 0.5, 1j*23.2/c*1250*2*pi, 0.5, 1j*34/c*1250*2*pi];
A0 = [real(x0(1)), imag(x0(2))/(2*pi*f/c); real(x0(3)), imag(x0(4))/(2*pi*f/c); real(x0(5)), imag(x0(6))/(2*pi*f/c)];

options.TolFun = 10^(-19);
options.MaxFunEvals = 10000;
options.MaxIter = 4000;

[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0, k , h,  -inf*ones(size(x0)), inf*ones(size(x0)), options);
output

A_est = [real(x(1)), imag(x(2))/(2*pi*f/c); real(x(3)), imag(x(4))/(2*pi*f/c); real(x(5)), imag(x(6))/(2*pi*f/c)];
plot(k,abs(F(x,k)),'k')

[A_org, A0, A_est]


