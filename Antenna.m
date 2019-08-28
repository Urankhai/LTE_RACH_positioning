% ant_vect = 4; %[100, 81, 64, 49, 36, 25, 16, 9, 4, 1];
% Ant_Num = ant_vect(1);
Ant_Dist = lambda/2;         % meters, the value has to be less than half of carrier's wavelength
disp(['Ant Num = ', num2str(Ant_Num)])

Ant_Num_x = Ant_Num;
Ant_Num_y = Ant_Num/Ant_Num_x;
Ant_array = AF_func_antenna(Ant_Num_x, Ant_Num_y, Ant_Dist);



ant_coord = Ant_array*OBS+ones(Ant_Num,1)*BS.loc; % rotation of the antenna by Oa

figure(100+Ant_Num)
plot(ant_coord(:,1),ant_coord(:,2),'x')
normal_vector_of_antenna = [BS.loc;BS.loc+100*OBS(2,:)];
plot(normal_vector_of_antenna(:,1),normal_vector_of_antenna(:,2),'linewidth',2)

BS.ant = ant_coord;
