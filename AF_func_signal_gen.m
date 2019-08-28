function s = AF_func_signal_gen(PUSCH)

Nfft = 2048;
first_cyclic = 160;
second_cyclic= 144;


N_smbls = size(PUSCH,2);
N_carri = size(PUSCH,1);
N_guard = 14;

fft_array = zeros(Nfft, N_smbls);

fft_array(N_guard + 1: N_guard + N_carri,:) = PUSCH;

pre_signal = ifft(fft_array, Nfft)*sqrt(Nfft);

s = zeros(1, N_smbls*Nfft + N_smbls/7*160 + N_smbls*6/7*144);

len = 0;
for l = 1:N_smbls
    if mod(l,7) == 1
        cyclic = pre_signal(end-first_cyclic+1:end,l);
    else
        cyclic = pre_signal(end-second_cyclic+1:end,l);
    end
    s(len + 1:len + length(cyclic) + 2048) = [reshape(cyclic, 1, []), reshape(pre_signal(:,l),1,[])];
    len = len + length(cyclic) + 2048;
end

% figure
% plot(abs(s(161:180)),'c','linewidth',3)
% hold on
% plot(abs(pre_signal(1:20,1)),'k')