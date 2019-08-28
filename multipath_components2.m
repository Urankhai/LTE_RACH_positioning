clear all;
close all;

% Global variables
global c F w lambda W

Normalization = 1;
K_grad = 180/pi;% Scaling coefficient
c = 3e8;
F = c*10;
lambda = c/F;   % in meters
w = 2*pi*F;
W = 30.72e6;


% LTE system parameters
w_LTE = 30.72e3;                        % the number of chips in 1ms
dt = 1/w_LTE/1000;                      % the duration of one chip


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

time_delays = zeros(L_paths, UE_Num);
atten_coeff = zeros(L_paths, UE_Num);

for k = 1:UE_Num
    atten_coeff(:,k) = sort(rand(L_paths,1),'descend');
    time_delays(:,k) = sort(rand(L_paths,1))/W;
    TD = randi(300,1);
    
    figure
    stem(time_delays(:,k), atten_coeff(:,k))
    figure
    hold on
    plot(exp(1j*2*pi*[0:0.01:1]))
    plot(upsample(exp(1j*w*time_delays(:,k)),2),'-o')
    
    
    exponenta = atten_coeff(:,k).*exp(1j*w*time_delays(:, k));
    tx_signal = [zeros(1,TD),carrier.*UE_uprach(k,:)];
    
    
    
    tx_sign(k,:) = sum(exponenta*tx_signal(1:w_LTE),1);
    rx_sign(k,:) = conj(carrier).*tx_sign(k,:);
end

% PRACH Detection
s1 = rx_sign(1:up:end);
for UEi = 1:UE_Num
    [out_s1, ~] = AF_func_correlator(CZ_seq(:,UEi), s1);
    figure
    subplot(2,1,1), plot3([1:2048],real(upsample(out_s1,2)),imag(upsample(out_s1,2)), '-o')
    subplot(2,1,2), plot(abs(out_s1), '-o')
end








