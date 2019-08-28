clear all
close all

w_LTE = 30.72e6;

ts = 1/w_LTE;

w = 1250;
Nfft = 1024;
GF_below = 12;
GF_above = 13;
length_CP = 3168;
length_GT = 2976;
length_seq = 24576;

carriers = (1:864);
samples = 1:24576;

waveforms = exp( -1i*2*pi * (w*(carriers-1)') * (ts*(samples-1)) );

[rach, seq] = AF_func_PRACH(7);

fseq = [zeros(GF_above,1); fft(seq)/sqrt(839); zeros(GF_below,1)];

SS = fft(seq)/sqrt(839);

s = fseq'*waveforms;%/sqrt(864);
shift1 = [zeros(1,40),s];
shift2 = [zeros(1,120),s];
signal = (1*s + 1*shift1(samples) + 1*shift2(samples))/1;


% figure
% hold on
% plot(real(s),'c','linewidth',3)
% plot((real(sqrt(24576)*rach(3169:end-2976))),'r','linewidth',1)

% plot(real(signal),'g','linewidth',1)
% plot((real((signal(1)/rach(1))*rach)),'r','linewidth',1)

% figure
% hold on
% plot(s(1:10)','c','linewidth',3)
% plot(sqrt(24576)*rach(length_CP+1:length_CP+10),'r','linewidth',1)
% plot(signal(1:10),'g','linewidth',1)

search_taps = -10:10;
t_srch = ts/100;
t_los = search_taps * t_srch;
t_nlos1 = search_taps * t_srch + ts;
t_nlos2 = search_taps * t_srch + 2*ts;
t_nlos3 = search_taps * t_srch + 3*ts;

t_mp = ts*(1:1024);%[t_los, t_nlos1, t_nlos2, t_nlos3];


% Time domain
in_s_wocp = signal(1:24:end);
% Frequency domain
incom_seq = fft(in_s_wocp, Nfft)/(Nfft);
YY = flipud(incom_seq(174:1012)');
HH = YY./SS;

F = exp(-1i*2*pi * (w*(carriers-1)') * (t_mp))/sqrt(864);
% figure
% hold on
% for k = 1:length(t_mp)
%     plot(angle(F(:,k)))
% end
% plot(angle(HH))
FF = F(14:end-12,:);
% HH = FF*x
% x_est = (FF'*FF)\FF'*HH;


% [out_s1, ~] = AF_func_correlator(CZ_seq(:,UEi), s1);

length_p = size(F,2);
% p = zeros(length_p,1);
% p(30) = 1/3;
p = (randn(length_p, 1) + 1i*randn(length_p, 1))/sqrt(2)/3;
figure
hold on
plot(real(p))

convergence = 0;
gamma = 1/sqrt(max(eig(FF'*FF)));
alpha = 3;
disp(['alpha * gamma = ', num2str(gamma*alpha)])
epsilon = 0.01;
t = 0;
while convergence == 0
    keyboard
    p_old = p;
    p_temp = p - gamma*FF'*(FF*p - HH);
    p = AF_func_sparsify(p_temp, gamma*alpha);
    plot(real(p))
    
    if norm(p-p_old) < epsilon
        convergence = 1;
    else
        t = t + 1;
        disp(['step = ', num2str(t)])
    end
end





