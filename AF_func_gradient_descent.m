function [P1, x, y, flag] = AF_func_gradient_descent(alpha, initial, impulse_signal, Ant_array, udt, c, w, iter)

P0 = initial(1);        % power at the current point
phi0 = initial(2);      % the value of cuttent phi
Ant0 = AF_func_turn(Ant_array, phi0/180*pi);

theta0 = initial(3);    % the value of cuttent theta

diff_step = 1/(iter);

phi1 = phi0 + diff_step;
Ant1 = AF_func_turn(Ant_array, phi1/180*pi);
theta1 = theta0 + diff_step;


Shifts01 = (Ant0(:,2) - min(Ant0(:,2)))*sin(theta1/180*pi);
Shifts10 = (Ant1(:,2) - min(Ant1(:,2)))*sin(theta0/180*pi);

% P(phi0, theta1) = P01, P(phi1, theta0) = P10
% Calculation the power on the pozition P01
P01 = -AF_func_powerprofile(Shifts01, impulse_signal, udt, c, w);
plot3(theta0 + diff_step, phi0, P01, 's', 'MarkerFaceColor','r', 'MarkerSize',5)
% P(phi0, theta1) = P01, P(phi1, theta0) = P10
% Calculation the power on the pozition P10
P10 = -AF_func_powerprofile(Shifts10, impulse_signal, udt, c, w);
plot3(theta0, phi0 + diff_step, P10, 's', 'MarkerFaceColor','g', 'MarkerSize',5)

diff01 = (P01 - P0)/(diff_step/180*pi);
diff10 = (P10 - P0)/(diff_step/180*pi);

% keyboard

x = (theta0/180*pi - alpha*diff01)/pi*180;
y = (phi0/180*pi - alpha*diff10)/pi*180;

Ant = AF_func_turn(Ant_array, y/180*pi);
Shifts = (Ant(:,2) - min(Ant(:,2)))*sin(x/180*pi);
P1 = -AF_func_powerprofile(Shifts, impulse_signal, udt, c, w);
plot3(x, y, P1, 's', 'MarkerFaceColor','c', 'MarkerSize',5)

epsilon = 0.0001;
if abs(theta0 - x) + abs(phi0 - y) <= 2*epsilon
    flag = 1; % stop criteria
else
    flag = 0; % continue
end





