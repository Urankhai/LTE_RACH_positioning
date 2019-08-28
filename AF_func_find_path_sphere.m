function [phi_max, thetta_max, shift_max, Smooth_power_profile] = AF_func_find_path_sphere(impulse_signal, search, Ant_array, udt, c, w)
% delta_degree = search.delta_degree;
search.delta_degree = 1;
search.phi_arr = (-180:search.delta_degree:179);% * pi/180;
search.thetta_arr = (0:search.delta_degree:90);% * pi/180;

phi_arr = search.phi_arr;% * pi/180;
thetta_arr = search.thetta_arr;% * pi/180;


power_profile = zeros(length(phi_arr), length(thetta_arr));
Smooth_power_profile = zeros(length(phi_arr), length(thetta_arr));
keyboard

Ant_Num = size(Ant_array,1);

Fn = 1/udt;                                         % Nyuist frequency, frequency of cut
step_Fn = Fn/length(impulse_signal);                % the number of Hertz between frequency steps
f = 0:step_Fn:(length(impulse_signal)-1)*step_Fn;



variant = 2; % if 1 --- function AF_func_correction is used, else --- not used
max = 0;
for l = 1:length(phi_arr)
    if mod(l,10)==0
        disp(['l = ', num2str(l)])
    end
    phi = phi_arr(l)/180*pi;
    Ant = AF_func_turn(Ant_array, phi);
    
    for k = 1:length(thetta_arr)
        thetta = thetta_arr(k)/180*pi;
        
        Shifts = (Ant(:,2) - min(Ant(:,2)))*sin(thetta);
        
%                 if (phi_arr(l) == -135) && (thetta_arr(k) == 44)
%                     keyboard
%                 end
        
        ant_signals = reshape(impulse_signal, Ant_Num, []);
        if variant == 1
            [~, s_sum] = AF_func_correction(Shifts/c, udt, ant_signals, w);
            power_profile(l,k) = s_sum*s_sum';
            

            if power_profile(l,k) > max
                max = power_profile(l,k);
                phi_max = phi;
                thetta_max = thetta;
                shift_max = Shifts;
            end
            
        else
            
            if norm(Shifts/c/udt - floor(Shifts/c/udt)) < 1e-12
                Geom_shift = udt*floor(Shifts/c/udt);
            else
                Geom_shift = Shifts/c;
            end
            discret_part = udt*ceil(Geom_shift/udt); % dt*N, where N is integer
            nondcrt_part = discret_part - Geom_shift;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Realization of Discrete shift
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            S = fft(ant_signals,[],2);
            
            Phases = discret_part*f;         % The phases are oposite to the phases from AF_func_MIMO
            Exp = exp(1j*2*pi*Phases);         % Calculation of frequency offset
            
            
            
            %         keyboard
            % f(t-a)  <=>  exp(-1iwa) F(w)
            % In order to perform time shift, we need to go to frequency domain via FFT
            MIMO_freq = Exp.* S; % frequency shifting
            % Going back to time domain
            MIMO_time = ifft(MIMO_freq, [], 2);
            
            % Smooth
            Smooth_Phases = Shifts/c*f;         % The phases are oposite to the phases from AF_func_MIMO
            Smooth_Exp = exp(1j*2*pi*Smooth_Phases);         % Calculation of frequency offset
            
            
            
            %         keyboard
            % f(t-a)  <=>  exp(-1iwa) F(w)
            % In order to perform time shift, we need to go to frequency domain via FFT
            Smooth_MIMO_freq = Smooth_Exp.* S; % frequency shifting
            % Going back to time domain
            Smooth_MIMO_time = ifft(Smooth_MIMO_freq, [], 2);
            Smooth_sum_MIMO = sum(Smooth_MIMO_time);
            Smooth_power_profile(l,k) = Smooth_sum_MIMO*Smooth_sum_MIMO';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % End of Realization of Discrete shift
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Realization of Doppler Effect and nondiscrete shift
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Nondscrt_shift = exp(-1j*w*nondcrt_part');
            % sum signals from all antennas
            sum_MIMO = Nondscrt_shift*MIMO_time;
            power_profile(l,k) = sum_MIMO*sum_MIMO';
            
            % realization of shift in time domain
            
            
            if power_profile(l,k) > max
                max = power_profile(l,k);
                phi_max = phi;
                thetta_max = thetta;
                shift_max = Shifts;
            end
        end
    end
end
disp(['phi = ',num2str(phi_max/pi*180),'; thetta = ', num2str(thetta_max/pi*180)])
% figure
% surf(thetta_arr, phi_arr, power_profile)
% axis equal

