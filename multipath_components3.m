clear all;
close all;

% Global variables
global c F w lambda W Nfft

Normalization = 1;
K_grad = 180/pi;% Scaling coefficient
c = 3e8;
F = c*10;
lambda = c/F;   % in meters
w = 2*pi*F;
W = 30.72e6; % sampling rate
Nfft = 1024; % fft size in LTE PRACH

% LTE system parameters
w_LTE = W/1000;                        % the number of chips in 1ms
dt = 1/W;                      % the duration of one chip
N_chips = w_LTE/Nfft;

% Duration of PRACH is 1ms for a conventional format (R <= 14km)
UE_Num = 1;
UE_rach = zeros(UE_Num, w_LTE);
CZ_seq = zeros(839,UE_Num);
for i = 1:UE_Num
    [rach, seq] = AF_func_PRACH(i);
    UE_rach(i,:) = transpose(rach);
    CZ_seq(:,i) = seq;
    
    clear seq rach
end


% upsampling in order to simulate analog signal
up = 1; % the carrier's wavelength ~= 0.976m
udt = dt/up;
% UE.ups = kron(UE.s,ones(1,up));
UE_uprach = kron(UE_rach,ones(1,up));

% duration of rach signal after upconverter is 30720*up samples
samples = 0:30720*up-1;
carrier = exp(1j * w * udt * samples);
tx_sign = zeros(size(UE_uprach));
rx_sign = zeros(size(UE_uprach));



% Number of dominant paths
% we set the maximum time delay for 300 meters as 30 integer 30.72MHz
% chips; multipath components are inside of the duration of one chip
L_paths = 3;
disp(['Multipath with ', num2str(L_paths+1), ' dominant paths'])
time_delays = zeros(L_paths, UE_Num);
chip_delays = zeros(L_paths, UE_Num);
atten_coeff = zeros(L_paths, UE_Num);
 CD_general = zeros(1, UE_Num);
 TD_general = zeros(1, UE_Num);

for k = 1:UE_Num
    atten_coeff(:,k) = sort(rand(L_paths,1),'descend');
    chip_delays(:,k) = sort(randi(N_chips, L_paths,1));
    CD_general(k) = randi(300,1);
    disp(['General CD = ', num2str(CD_general(k)), ' chips; Multipath CDs = [', num2str(chip_delays(:,k)'), '] chips;'])
    
    time_delays(:,k) = chip_delays(:,k)/W;
    TD_general(k) = CD_general(k)/W;
    disp(['General TD = ', num2str(TD_general(k)), '; Multipath TDs = [', num2str(time_delays(:,k)'), '];'])
    
    figure
    stem([TD_general(k); TD_general(k)+time_delays(:,k)], [1; atten_coeff(:,k)])
    figure
    hold on
    plot(exp(1j*2*pi*[0:0.01:1]))
    plot(upsample(exp(1j*w*[TD_general(k); TD_general(k)+time_delays(:,k)]),2),'r-o')
    plot(upsample(exp(1j*w*TD_general(k)),2),'g-o','linewidth',2)

    tx_signal = [zeros(1,CD_general),carrier.*UE_uprach(k,:)];
    
    temp_signal = tx_signal(1:w_LTE);
    for paths = 1:L_paths
        temp_s = [zeros(1, chip_delays(paths, k)), tx_signal];
        temp_signal = temp_signal + atten_coeff(paths,k)*temp_s(1:w_LTE);
    end
    
    
    tx_sign(k,:) = temp_signal;
    rx_sign(k,:) = conj(carrier).*tx_sign(k,:);
    
    figure
    plot(abs(UE_uprach(k,:)))
    hold on
    plot(abs(rx_sign(k,:)))
    title(['General CD = ', num2str(CD_general(k)), ' chips; Multipath CDs = [', num2str(chip_delays(:,k)'), '] chips;'])
end

% PRACH Detection
s1 = rx_sign(1:up:end);
for UEi = 1:UE_Num
    [out_s1, ~] = AF_func_correlator(CZ_seq(:,UEi), s1);
    figure
    subplot(2,1,1), axis equal, plot3([1:2048],real(upsample(out_s1,2)),imag(upsample(out_s1,2)), '-o')
    subplot(2,1,2), plot(abs(out_s1), '-o')
end

% Reflection, penetration and other interactions with obstacles attenuate
% the signal's power. Even more, the attenuation coefficients caused by 
% interactions are always complex valued. This creates an ambiguity in time
% estimation, since it is complex valued, and the accuracy of time
% estimation can be delusioned up to 2 pi, which is in distance
% representation equals to the wavelength. From this observation, we make a
% conclusion that the accuracy of any localization using interacted signals
% cannot be better than the wavelength of the signal's carrier. This
% statement works for the case when the physical properties of materials
% are not accurate or unknown. 

% The ambiguity of time estimation can be explained by physical
% properties of materials with which signals are interacting. For the first
% stage, signals come to material and are absorbed by the material, then due to
% physical properties of the material some part of the signals' energy is
% emmiting from the material creating almost the attenuated copy of the
% incident signals. The absorbtion and the emission require some time, from
% which, obviously, the attenuation coefficients are complex valued. 







