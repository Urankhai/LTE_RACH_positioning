function [out_s, out_t] = AF_func_correlator(seq, in_s)

Nfft = 1024;
GF_below = 12;
GF_above = 13;
length_CP = 3168;
length_GT = 2976;

% keyboard

fseq = [zeros(GF_above,1); fft(seq)/sqrt(size(seq,1)); zeros(GF_below,1)];
local_seq = sqrt(Nfft)*ifft(fseq,Nfft);

% asd = [local_seq(end-29:end),local_seq(1:end-30)];
fft_local = fft(local_seq)/sqrt(Nfft);

in_s_wocp = sqrt(24)*in_s(length_CP+1:24:end-length_GT);
incom_seq = fft(in_s_wocp, Nfft)/sqrt(Nfft);

% [test_rach, test_seq] = AF_func_PRACH(1);
% test_preamble = test_rach(length_CP+1:24:end-length_GT);


figure(200)
hold on
plot(abs(local_seq),'g','linewidth',3)
plot(abs(in_s_wocp),'r')

figure(201)
hold on
plot(local_seq(1:10),'g','linewidth',3)
plot(in_s_wocp(1:10),'r')

corr_seq = fft_local'.*incom_seq;
out_s = Nfft*ifft(corr_seq, Nfft);
[~, out_t] = max(abs(out_s));
keyboard

