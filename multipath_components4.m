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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of RACH signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rach, seq] = AF_func_PRACH2(1);
UE_rach = transpose(rach);
CZ_seq = seq;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End RACH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% upsampling in order to simulate analog signal
up = 1; % the carrier's wavelength ~= 0.976m
udt = dt/up;
% UE.ups = kron(UE.s,ones(1,up));
UE_uprach = kron(UE_rach,ones(1,up));

% duration of rach signal after upconverter is 30720*up samples
samples = 0:30720*up-1;
carrier = exp(1j * w * udt * samples);




% Number of dominant paths
% we set the maximum time delay for 300 meters as 30 integer 30.72MHz
% chips; multipath components are inside of the duration of one chip
L_paths = 0;
disp(['Multipath with ', num2str(L_paths+1), ' dominant paths'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of Multipath Environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atten_coeff = sort(rand(L_paths,1),'descend'); % coefficient of attenuation
chip_delays = sort(randi(N_chips, L_paths,1)); % integer number of delays

% Set up the LoS time delay
chip_iter = 0;
phase_diff = zeros(1,24);
for CD_general = 21:27;%randi(300,1); % LoS path's delay is an integer number
    chip_iter = chip_iter + 1;
disp(['General CD = ', num2str(CD_general), ' chips; Multipath CDs = [', num2str(chip_delays'), '] chips;'])

% Calculation time delays accordingly to the chip's delays
time_delays = chip_delays/W;
TD_general = CD_general/W;
disp(['General TD = ', num2str(TD_general), '; Multipath TDs = [', num2str(time_delays'), '];'])

figure(1) % Drwaing of multipath delays and the coefficients of attenuation
stem([TD_general; TD_general + time_delays], [1; atten_coeff])

figure(2) % Drawing of the phase rotation
hold on
plot(exp(1j*2*pi*[0:0.01:1]))
plot(upsample(exp(-1j*w*[TD_general; TD_general + time_delays]),2),'b-o')
plot(upsample(exp(-1j*w*TD_general),2),'-o','linewidth',2)
title(['Argument of the TD = ' , num2str(180*angle(exp(-1j*w*TD_general))/pi)])
arg_TD = 180*angle(exp(-1j*w*TD_general))/pi;
if arg_TD < 0 
    arg_TD = 360 + arg_TD;
end

% Set up the LoS's time delay
tx_signal = [zeros(1,CD_general),carrier.*UE_uprach];
% Set up the NLoSs' time delays
temp_signal = tx_signal(1:w_LTE);
for paths = 1:L_paths
    temp_s = [zeros(1, chip_delays(paths)), tx_signal];
    temp_signal = temp_signal + atten_coeff(paths)*temp_s(1:w_LTE);
end
tx_sign = temp_signal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Multipath Environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_sign = conj(carrier).*tx_sign;

figure(3)
plot(abs(UE_uprach))
hold on
plot(abs(rx_sign))
title(['General CD = ', num2str(CD_general), ' chips; Multipath CDs = [', num2str(chip_delays'), '] chips;'])
disp(['The phase shift if the carrier = ' , num2str(arg_TD)])


% PRACH Detection
s1 = rx_sign(1:up:end);

[out_s1, ~] = AF_func_correlator2(CZ_seq, s1);
max_Elm = max(out_s1);
if abs(180*angle(max_Elm)/pi - round(180*angle(max_Elm)/pi)) <1e-10
    arg_Elm = round(180*angle(max_Elm)/pi);
else
    arg_Elm = 180*angle(max_Elm)/pi;
end
if arg_Elm < 0
    arg_Elm = 360 + arg_Elm;
end

figure(4)
hold on
axis equal
plot(real(upsample(out_s1(1:30),1)),imag(upsample(out_s1(1:30),1)), '-o')
figure(5)
plot(abs(out_s1(1:30)), '-o')

disp(['       Argument of the max Elm = ', num2str(arg_Elm)])

asd = arg_TD - arg_Elm;
if  asd< 0
    asd = -asd;
end

disp(['Difference between the phase shift of the carrier and the phase of correlation = ', num2str((asd/180*pi))])
figure(6)
hold on
title(['Argument of the max Elm = ', num2str(arg_Elm)])
axis([-839 839 -839 839])
plot([0, max_Elm],'-o')

figure(2)
plot([0, max_Elm/norm(max_Elm)],'-o','linewidth', 2)
legend('Unit circle', 'Multipath delays', 'LoS delay', 'Corr phase')

phase_diff(chip_iter) = asd;
end
figure
plot(phase_diff)

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






