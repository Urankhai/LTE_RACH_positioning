function [out_s, out_t] = AF_func_correlator(seq, in_s)

Nfft = 1024;
GF_below = 12;
GF_above = 13;
length_CP = 3168;
length_GT = 2976;

% keyboard

fseq = [zeros(GF_above,1); fft(seq)/sqrt(size(seq,1)); zeros(GF_below,1)];
% Time domain
local_seq = sqrt(Nfft)*ifft(fseq,Nfft);
% Frequency domain
fft_local = fft(local_seq)/sqrt(Nfft);

% keyboard
% Time domain
in_s_wocp = sqrt(24)*in_s(length_CP+1:24:end-length_GT);
% Frequency domain
incom_seq = fft(in_s_wocp, Nfft);



% figure(200)
% hold on
% plot(abs(local_seq),'g','linewidth',3)
% plot(abs(in_s_wocp),'r')

corr_seq = fft_local'.*incom_seq;
out_s = sqrt(Nfft)*ifft(corr_seq, Nfft);
% keyboard
[~, out_t] = max(abs(out_s));

