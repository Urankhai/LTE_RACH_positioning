function s = AF_func_signal_gen_mult(N_smbls, fft_array)

Nfft = 2048;
first_cyclic = 160;
second_cyclic= 144;

length_GT = 2976; % the same Guard Time as for PRACH

pre_signal = ifft(fft_array, Nfft)*sqrt(Nfft);
s = zeros(1, N_smbls*Nfft + N_smbls/7*160 + N_smbls*6/7*144 + length_GT);

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

