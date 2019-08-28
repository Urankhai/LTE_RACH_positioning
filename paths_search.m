% load('rec.mat')
% output UE1_rach UE2_rach UE3_rach UE1_sign UE2_sign UE3_sign
% search1 = UE1_rach;
% search2 = UE2_rach;
% search3 = UE3_rach;

search1 = UE.rach(1,:);
search2 = UE.rach(2,:);
search3 = UE.rach(3,:);

x=size(search1,1);
y=size(search1,2);
z=size(search1,3);
% signal from one user
UE1 = zeros(x,z);
UE2 = UE1;
UE3 = UE1;
for i = 1:y
    temp_M1 = reshape(search1(:,i,:),x,z);
    temp_M2 = reshape(search2(:,i,:),x,z);
    temp_M3 = reshape(search3(:,i,:),x,z);
    UE1 = UE1 + temp_M1;
    UE2 = UE2 + temp_M2;
    UE3 = UE3 + temp_M3;
end
M = reshape(search1(:,1,:),x,z);%+40*UE2+40*UE3;
% M = reshape(search2(:,1,:),x,z)+reshape(search3(:,1,:),x,z);

GT = 1;
%GT:10:GT+6000;
LTE_sampling = GT:10:307200;
rand_seq = LTE_sampling(randperm(size(LTE_sampling,2),300));
impulse_signal = M(:,rand_seq);

% test plate borders

lattice_step = 0.5; % How is small the checking step


edges = 5;
x_min = -edges; x_max = edges;
y_min = -edges; y_max = edges;
z_plate = 1;

x_plate = x_min:lattice_step:x_max;
x_steps = size(x_plate,2);

y_plate = y_min:lattice_step:y_max;
y_steps = size(y_plate,2);

power_profile = zeros(x_steps, y_steps);
max_val = 0;
for x_ind = 1:x_steps
    for y_ind = 1:y_steps
        VS_vect = -[x_plate(x_ind), y_plate(y_ind), z_plate];
        VS_orth = VS_vect/norm(VS_vect);
        
        theta_plate = acos(-VS_orth(3));
        if VS_orth(3) == 1
            phi_plate = 0;
        else
            phi_plate = atan2(VS_orth(2),VS_orth(1));
        end
        
        
        Ant = AF_func_turn(Ant_array, phi_plate);
        Shifts = (Ant(:,2) - min(Ant(:,2)))*sin(theta_plate);
        test_P = AF_func_powerprofile(Shifts, impulse_signal, udt, c, w);
        
        %         % collecting signals from all antenna elements
        %         Geom_shift = Shifts/c;
        %         discret_part = ceil(Geom_shift/udt); % dt*N, where N is integer
        %         discret_r = zeros(size(impulse_signal));
        %         for lk = 1:size(Ant_array,1);
        %             discret_r(lk,:) = [zeros(1,discret_part(lk)), impulse_signal(lk,1:end-discret_part(lk))];
        %         end
        %
        %         nondcrt_part = udt*discret_part - Geom_shift;
        %         Nondscrt_shift = exp(1j*w*nondcrt_part);
        %         ant_sum = transpose(Nondscrt_shift)*discret_r;
        
        
        power_profile(x_ind,y_ind) = test_P;
        %         disp(['difference = ', num2str(power_profile(x_ind,y_ind) - test_P)])
        
        if power_profile(x_ind,y_ind) > max_val
            max_val = power_profile(x_ind,y_ind);
            P_max = max_val;
            phi_max = phi_plate;
            theta_max = theta_plate;
            Max_point = [-x_plate(x_ind), -y_plate(y_ind), max_val];
            shift_max = Shifts;
        end
    end
    
end
disp(['theta_max = ',num2str(theta_max/pi*180),'; phi_max = ',num2str(phi_max/pi*180)])
figure
hold on
surf(-x_plate,-y_plate,power_profile')
plot3(Max_point(1),Max_point(2),Max_point(3),'s', 'MarkerFaceColor','c', 'MarkerSize',10)

%%
iteration = 0;
initial = [Max_point(1),Max_point(2)];%
theta_initial = theta_max;
phi_initial = phi_max;
% initial = [0.5, -0.5];
alpha = 2;
P0 = Max_point(3);
grad_descent_arr = [];
while iteration < 100
    grad_descent_arr = [grad_descent_arr; [theta_initial/pi*180, phi_initial/pi*180, P0, initial]];
    iteration = iteration + 1;
    desc_step = lattice_step/(2*iteration);
    dx = initial + [1, 0] * desc_step;
    dy = initial + [0, 1] * desc_step;
    
    % calculation power profile for dx
    VS_dx = [dx, z_plate]/norm([dx, z_plate]);
    theta_dx = acos(-VS_dx(3));
    if VS_dx(3) == 1
        phi_dx = 0;
    else
        phi_dx = atan2(VS_dx(2),VS_dx(1));
    end
    Ant_dx = AF_func_turn(Ant_array, phi_dx);
    Shifts_dx = (Ant_dx(:,2) - min(Ant_dx(:,2)))*sin(theta_dx);
    P_dx = AF_func_powerprofile(Shifts_dx, impulse_signal, udt, c, w);
    plot3(dx(1), dx(2), P_dx, 's', 'MarkerFaceColor','y', 'MarkerSize',7)
    
    % calculation power profile for dx
    VS_dy = [dy, z_plate]/norm([dy, z_plate]);
    theta_dy = acos(-VS_dy(3));
    if VS_dy(3) == 1
        phi_dy = 0;
    else
        phi_dy = atan2(VS_dy(2),VS_dy(1));
    end
    Ant_dy = AF_func_turn(Ant_array, phi_dy);
    Shifts_dy = (Ant_dy(:,2) - min(Ant_dy(:,2)))*sin(theta_dy);
    P_dy = AF_func_powerprofile(Shifts_dy, impulse_signal, udt, c, w);
    plot3(dy(1), dy(2), P_dy, 's', 'MarkerFaceColor','g', 'MarkerSize',3)
    
    diff_dx = log10(P0) - log10(P_dx);
    diff_dy = log10(P0) - log10(P_dy);
    
%     keyboard
    
    initial = initial - alpha*[diff_dx, diff_dy];
    
    % calculation power profile for the gradient descent point
    VS_initial = [initial, z_plate]/norm([initial, z_plate]);
    theta_initial = acos(VS_initial(3));
    if VS_initial(3) == 1
        phi_initial = 0;
    else
        phi_initial = atan2(VS_initial(2),VS_initial(1));
    end
    Ant_initial = AF_func_turn(Ant_array, phi_initial);
    Shifts_initial = (Ant_initial(:,2) - min(Ant_initial(:,2)))*sin(theta_initial);
    
    P0 = AF_func_powerprofile(Shifts_initial, impulse_signal, udt, c, w);
    plot3(initial(1), initial(2), P0, 's', 'MarkerFaceColor','y', 'MarkerSize',3)
    
    
    
    disp(['iteration = ',num2str(iteration), '; theta = ',num2str(theta_initial/pi*180),'; phi = ',num2str(phi_initial/pi*180)])
    disp(['iteration = ',num2str(iteration), '; P = ', num2str(P0), '; [dx, dy] = ',num2str(alpha*[diff_dx, diff_dy])])
end
plot3(grad_descent_arr(:,4), grad_descent_arr(:,5), grad_descent_arr(:,3), '-o', 'MarkerFaceColor','r', 'MarkerSize',7)


% LoS angles: theta = 17.077; phi = -72.2415
% LoS angles: theta = 27.7592; phi = -145.8336
% LoS angles: theta = 52.7218; phi = -163.0913

figure
subplot(3,1,1), plot(grad_descent_arr(:,1)) 
hold on 
line([0 100],[17.07 17.07])
line([0 100],[27.7592 27.7592])
line([0 100],[52.7218 52.7218])
subplot(3,1,2), plot(grad_descent_arr(:,2))
hold on 
line([0 100],[-72.24 -72.24])
line([0 100],[-145.8336 -145.8336])
line([0 100],[-163.0913 -163.0913])
subplot(3,1,3), plot(grad_descent_arr(:,3))














