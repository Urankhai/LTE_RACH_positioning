% plane scenario, rectangular area is sized as 500 meters in width and 300
% meters in length

colors = {'--c', '--r', '--g'};
figure(100+Ant_Num)
X_len = 500;
Y_len = 500;
axis equal
% axis([0 X_len 0 Y_len])
hold on
grid on

rx_p = [0, 0, 1];                     % Polarization vector of BS
rx_p = rx_p/norm(rx_p);
BS.pol = rx_p;

BS.loc = [X_len/2, Y_len - 200];                          % BS's position
plot(BS.loc(1), BS.loc(2),'s', 'MarkerFaceColor','g','MarkerSize',3)

OBS = [-1 0; 0 -1];

% UE_Num is defined in main function
% pUe = [X_len*rand(UE_Num,1), 0.8*Y_len*rand(UE_Num,1)];
% pUe = [245+10*rand(UE_Num,1), [300-10*rand(1); 300-10*rand(1); 300-10*rand(1)]];

% radius = [3000.05;9000.05;15000.05];%
% vect_dist = 30;
% vect_var =  100;

R_mean = 30;    %vect_dist(scenario_ind);
 R_var = 100;   %vect_var(scenario_ind);


radius = R_mean + R_var*rand(1,1);
tot_ang = pi/3;
angles = tot_ang*rand(1,1);

UE_org = [250 + radius*sin(angles - tot_ang/2),300 - radius*cos(angles - tot_ang/2)];
line_array = -10:10;
% est_p = UE_org +[0,3];
% plot(est_p(1),est_p(2),'*')

%distance from user to the reflectors
d1 = 15;   % ground rerflection
rUE1 = UE_org - sqrt(2)*[d1,d1];
line1 = UE_org - (1/sqrt(2))*[d1,d1];

a1 = -line1(1)-line1(2);
line_array1 = -(line_array + line1(1)) - a1;
plot(line_array + line1(1), line_array1)

d2 = 40;     % reflection from the nearest building
rUE2 = UE_org - sqrt(2)*[-d2,d2];
line2 = UE_org - (1/sqrt(2))*[-d2,d2];

a2 = line2(1)-line2(2);
line_array2 = (line_array + line2(1)) - a2;
plot(line_array + line2(1), line_array2)


% % calculation of reflection
% sp1 = [est_p(1),est_p(2) + a1]*[-1;1]/norm([est_p(1),est_p(2) + a1])/sqrt(2);
% dd1 = norm([est_p(1),est_p(2) + a1])*sin(acos(sp1));
% ref1 = est_p + sqrt(2)*dd1*[-1, -1];
% plot(ref1(1),ref1(2), 'd')
% 
% sp2 = [est_p(1),est_p(2) + a2]*[1;1]/norm([est_p(1),est_p(2) + a2])/sqrt(2);
% dd2 = norm([est_p(1),est_p(2) + a2])*sin(acos(sp2));
% ref2 = est_p + sqrt(2)*dd2*[1, -1];
% plot(ref2(1),ref2(2), 'd')




pUe = [UE_org; rUE1; rUE2];
% radius = [1.660528024863594*10^2;0.193374253950820*10^2];
% angles = [0.237863368557995; 0.877500220634629];

% radius = [1.075382056891205*10^2; 1.806839292570790*10^2];
% angles = [0.896223816009413; 0.951484865122977];

% pUe = [250 + radius.*sin(angles - tot_ang/2),300 - radius.*cos(angles - tot_ang/2)];

UE.loc = zeros(UE_Num, 2);
UE.pol = zeros(UE_Num, 3);

for i = 1:UE_Num
    plot(pUe(i,1),pUe(i,2),'s', 'MarkerFaceColor','r','MarkerSize',3)
%     disp(['UE', num2str(i),' location ', num2str(pUe(i,:))])
    
    % Users location 
    UE.loc(i,:) = pUe(i,:);
    LoS_vect = (BS.loc - pUe(i,:));
    LoS_dist = norm(BS.loc - pUe(i,:));
    plot([BS.loc(1), pUe(i,1)], [BS.loc(2), pUe(i,2)], colors{i})
    
    LoS_dec = -OBS*LoS_vect'/LoS_dist;
    
    theta = 180*acos(LoS_dec(1))/pi;
    
%     disp(['LoS angles: theta = ',num2str(theta), '; LoS dist = ', num2str(LoS_dist)])
    
    % Users polarization
    tx_p = [0, 0, 1];                     % Polarization vector of UE
    UE.pol(i,:) = tx_p/norm(tx_p);
end
% keyboard
% LoS angles: theta = 17.077; phi = -72.2415
% LoS angles: theta = 27.7592; phi = -145.8336
% LoS angles: theta = 52.7218; phi = -163.0913