function Ant_array = AF_func_antenna(Num_x, Num_y, dist)
% MIMO Rx Antenna configuration

Ant_Num = Num_x*Num_y;
% Ant_Dist = dist;         % meters

Ant_array = zeros(Ant_Num, 2); % The coordinates of antenna elements

k = 0;
for i = 1:Num_x
    for j = 1:Num_y
        k = i + (j - 1)*Num_x;
        Ant_array(k,:) = [i*dist, j*dist] - [0.5*(Num_x + 1)*dist, 0.5*(Num_y + 1)*dist];
    end
end

% figure(101)
% axis square
% hold on
% grid on
% axis([-0.5 0.5 -0.5 0.5])
% plot(Ant_array(:,1), Ant_array(:,2),'*')