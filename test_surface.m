

z = temp_phases';
H = [(1:Ant_Num)', ones(Ant_Num,1)];
x = (H'*H)^(-1)*H'*z;
zz = H*x;

% figure
% hold on
% plot(temp_phases,'-o')
% plot(zz,'-o')
% 
% figure
% hold on
% plot(zeros(size(zz)))
% diff_phase = temp_phases - zz';
% plot(diff_phase,'-*')
% 
% 
% 
% 
% figure
% hold on
% plot(zeros(size(zz)))
% plot(diff(temp_phases - zz'),'-*')

% new_ph = zeros(size(temp_phases));

good_set = 8:Ant_Num;
Good_Ant_Num = length(good_set);
new_ph = temp_phases(good_set);