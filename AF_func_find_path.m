function [phi_max, thetta_max, shift_max, power_profile, P_max] = AF_func_find_path(impulse_signal, Ant_array, udt, c, w, phi_arr, theta_arr)

Num_ant = size(Ant_array,1);
% delta_degree = delta_degree;
% delta_degree = 30;


power_profile = zeros(length(phi_arr), length(theta_arr));
% Smooth_power_profile = zeros(length(phi_arr), length(theta_arr));
% keyboard

% Fn = 1/udt;                                         % Nyuist frequency, frequency of cut
% step_Fn = Fn/length(impulse_signal);                % the number of Hertz between frequency steps
% f = 0:step_Fn:(length(impulse_signal)-1)*step_Fn;



max_val = 0;
polarplot_step = 0;
for l = 1:length(phi_arr)
%     if mod(l,10)==0
%         disp(['l = ', num2str(l)])
%     end
    phi = phi_arr(l)/180*pi;
    Ant = AF_func_turn(Ant_array, phi);
    
    for k = 1:length(theta_arr)
        thetta = theta_arr(k)/180*pi;
%         keyboard
        Shifts = (Ant(:,2) - min(Ant(:,2)))*sin(thetta);
        
        
        
        if power_profile(l,k) > max_val
            max_val = power_profile(l,k);
            P_max = max_val;
            phi_max = phi;
            thetta_max = thetta;
            shift_max = Shifts;
        end
        
        Geom_shift = Shifts/c;
        discret_part = ceil(Geom_shift/udt); % dt*N, where N is integer
        
        discret_r = zeros(size(impulse_signal));
        for lk = 1:Num_ant
            discret_r(lk,:) = [zeros(1,discret_part(lk)), impulse_signal(lk,1:end-discret_part(lk))];
        end
        
%         keyboard
        
        
                
        nondcrt_part = udt*discret_part - Geom_shift;
        Nondscrt_shift = exp(1j*w*nondcrt_part);
        ant_sum = transpose(Nondscrt_shift)*discret_r;
        
        
        power_profile(l,k) = ant_sum*ant_sum';
        polarplot_step = polarplot_step + 1;
%         x(polarplot_step) = thetta*sin(phi);
%         y(polarplot_step) = thetta*cos(phi);
%         z(polarplot_step) = power_profile(l,k);
        
        % realization of shift in time domain
        
        
        if power_profile(l,k) > max_val
            max_val = power_profile(l,k);
            P_max = max_val;
            phi_max = phi;
            thetta_max = thetta;
            shift_max = Shifts;
        end
        
    end
end
% disp(['phi = ',num2str(phi_max/pi*180),'; thetta = ', num2str(thetta_max/pi*180)])
figure
% keyboard
surf(theta_arr, phi_arr, -power_profile)
% axis equal
% figure
% plot3(x,y,z)
