
UE1_rach = reshape(r_rach(:,1,:,:), Ant_Num, size(planes,1) + 1, []);
UE2_rach = reshape(r_rach(:,2,:,:), Ant_Num, size(planes,1) + 1, []);
UE3_rach = reshape(r_rach(:,3,:,:), Ant_Num, size(planes,1) + 1, []);

UE1_sign = reshape(r_sign(:,1,:,:), Ant_Num, size(planes,1) + 1, []);
UE2_sign = reshape(r_sign(:,2,:,:), Ant_Num, size(planes,1) + 1, []);
UE3_sign = reshape(r_sign(:,3,:,:), Ant_Num, size(planes,1) + 1, []);

% save('rec.mat', 'UE1_rach', 'UE1_sign', 'UE2_rach', 'UE2_sign', 'UE3_rach', 'UE3_sign')
