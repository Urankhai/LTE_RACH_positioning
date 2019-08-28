number_of_calculations = 0;
for iteration = 1:4
    if iteration == 1
        delta_degree = 15;
        phi_arr = (-180:delta_degree:179);% * pi/180;
        theta_arr = (0:delta_degree:90);% * pi/180;
    elseif iteration == 2
        delta_degree = 5;
        phi_arr = phi_max*180/pi + (-15:delta_degree:15);% * pi/180;
        theta_arr = theta_max*180/pi + (-15:delta_degree:15);% * pi/180;
    elseif iteration == 3
        delta_degree = 1;
        phi_arr = phi_max*180/pi + (-5:delta_degree:5);% * pi/180;
        theta_arr = theta_max*180/pi + (-5:delta_degree:5);% * pi/180;
    else
        delta_degree = 0.1;
        phi_arr = phi_max*180/pi + (-1:delta_degree:1);% * pi/180;
        theta_arr = theta_max*180/pi + (-1:delta_degree:1);% * pi/180;
    end
    
    [phi_max, theta_max, shift_max, power_profile, P_max] = AF_func_find_path(impulse_signal, Ant_array, udt, c, w, phi_arr, theta_arr);
    disp(['Rough Power = ', num2str(P_max),'; theta = ',num2str(theta_max*180/pi), '; phi = ',num2str(phi_max*180/pi)])
    % disp(['number of operations = ', num2str(size(phi_arr,2)*size(theta_arr,2))])
    number_of_calculations = number_of_calculations + size(phi_arr,2)*size(theta_arr,2);
    
%     x = theta_max*180/pi;
%     y = phi_max*180/pi;
%     initial = [-P_max, y, x];
%     hold on
%     plot3(x,y,-P_max, 's','MarkerFaceColor', 'm', 'MarkerSize',10)
%     alpha = 5*10^13;
%     
%     results = [];
%     flag = 0;
%     %for l = 1:1000
%     time = 0;
%     while flag ~= 1
%         time = time + 1;
%         [P1, x, y, flag] = AF_func_gradient_descent(alpha, initial, impulse_signal, Ant_array, udt, c, w, time);
%         disp(['Power = ',num2str(P1), '; theta = ', num2str(x), '; phi = ',num2str(y)])
%         results = [results; [x, y, P1]];
%         % plot3(x, y, P1, 'or')
%         
%         initial = [P1, y, x];
%     end
%     figure
%     subplot(3,1,1), plot(results(:,1))
%     subplot(3,1,2), plot(results(:,2))
%     subplot(3,1,3), plot(results(:,3))
%     keyboard
end
disp(['in total we spent ', num2str(number_of_calculations), ' unit calculations instead of ',num2str(360*90)])
disp('  ')
disp(['Final thetta = ', num2str(theta_max/pi*180),'; phi = ', num2str(phi_max/pi*180)])
% figure
% plot3(x,y,z,'o')