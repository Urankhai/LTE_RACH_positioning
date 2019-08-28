close all

sss1 = [zeros(1,1000),rach(1:end-1000)];
sss2 = [zeros(1,1001),rach(1:end-1001)];

rec_sss1_noise = sss1 + noise;
rec_sss1_fft = fft(rec_sss1_noise);
rec_sss1_fft_clear = zeros(size(rec_sss1_fft));
dw_freq = 19526;%20160;
up_freq = 20590;%21250;
rec_sss1_fft_clear(dw_freq:up_freq) = rec_sss1_fft(dw_freq:up_freq);

rec_sss1 = ifft(rec_sss1_fft_clear);

figure
hold on
plot(abs(fft(sss1)))
plot(abs(fft(noise)),'r')
plot(abs(fft(rec_sss1)),'g')

rec_sss2_noise = sss2 + noise;
rec_sss2_fft = fft(rec_sss2_noise);
rec_sss2_fft_clear = zeros(size(rec_sss2_fft));
dw_freq = 19526;%20160;
up_freq = 20590;%21250;
rec_sss2_fft_clear(dw_freq:up_freq) = rec_sss2_fft(dw_freq:up_freq);

rec_sss2 = ifft(rec_sss2_fft_clear);

figure
hold on
plot(abs(fft(sss2)))
plot(abs(fft(noise)),'r')
plot(abs(fft(rec_sss2)),'g')

rec_sss1_corr = rec_sss1.*conj(carrier);
rec_sss2_corr = rec_sss2.*conj(carrier);

% phases measurement
[out_sss1, ~] = AF_func_correlator(CZ_seq(:,1), rec_sss1_corr);
[out_sss2, ~] = AF_func_correlator(CZ_seq(:,1), rec_sss2_corr);

figure
hold on
plot(abs(out_sss1),'linewidth', 3)
plot(abs(out_sss2),'r')

figure
hold on
plot((out_sss1))
plot((out_sss2),'r')

angle11 = angle(max(out_sss1))
angle21 = angle(max(out_sss2))



figure
hold on
plot([0,max(out_sss1)],'-o')
plot([0,max(out_sss2)],'r-o')

correlation = zeros(1,length(rach));
for i =1:length(rach)
    qwe=circshift(rach, [1, i]);
    correlation(i) = (sss1+1000*noise)*qwe';
end
figure
plot(abs(correlation))

[A1,I1] = max(correlation);