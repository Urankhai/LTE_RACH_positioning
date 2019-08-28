function [out_s, out_t] = AF_func_RACH_correlator(s, in_s)

Nfft = 1024;
length_CP = 3168;
length_GT = 2976;

% keyboard
% Time domain
s_wocp = sqrt(24)*s(length_CP+1:24:end-length_GT);
in_s_wocp = sqrt(24)*in_s(length_CP+1:24:end-length_GT);
% Frequency domain
seq = fft(s_wocp, Nfft);
incom_seq = fft(in_s_wocp, Nfft);



% figure(200)
% hold on
% plot(abs(local_seq),'g','linewidth',3)
% plot(abs(in_s_wocp),'r')

corr_seq = seq.*incom_seq;
out_s = sqrt(Nfft)*ifft(corr_seq, Nfft);
[~, out_t] = max(abs(out_s));

