clear all
close all
nfft = 1024;
% Duration of PRACH is 1ms for a conventional format (R <= 14km)
[rach, seq] = AF_func_PRACH(14); % since UE_Num is the number of reflectors

GF_below = 12;
GF_above = 13;
fseq = [zeros(GF_above,1); fft(seq)/sqrt(length(seq)); zeros(GF_below,1)]; % it takes 6 RBs
preamble = sqrt(nfft)*ifft(fseq, nfft);
fft_preamble = fft(preamble);


fid = fopen('zc_time.dat','w');
for k = 1:length(preamble)
    fwrite(fid,real(preamble(k)), 'float32');
    fwrite(fid,imag(preamble(k)), 'float32');
    disp(['step = ',num2str(k)])
end
fclose(fid);

fid = fopen('zc_freq.dat','w');
for k = 1:length(fft_preamble)
    fwrite(fid,real(fft_preamble(k)), 'float32');
    fwrite(fid,imag(fft_preamble(k)), 'float32');
    %     disp(['step = ',num2str(k)])
end
fclose(fid);

fid = fopen('zc_time_int16.dat','w');
for k = 1:length(preamble)
    fwrite(fid,real(preamble(k)), 'int16');
    fwrite(fid,imag(preamble(k)), 'int16');
    disp(['step = ',num2str(k)])
end
fclose(fid);

fid = fopen('zc_freq_int16.dat','w');
for k = 1:length(fft_preamble)
    fwrite(fid,real(fft_preamble(k)), 'int16');
    fwrite(fid,imag(fft_preamble(k)), 'int16');
    %     disp(['step = ',num2str(k)])
end
fclose(fid);

fid = fopen('zc_time.txt','w');
for k = 1:length(preamble)
    %     keyboard
    if k < length(preamble)
        if imag(preamble(k)) >= 0
            fprintf(fid,'%.11f +% .11fj,\t',[real(preamble(k)),imag(preamble(k))]);
        else
            fprintf(fid,'%.11f -% .11fj,\t',[real(preamble(k)),-imag(preamble(k))]);
        end
    else
        if imag(preamble(k)) >= 0
            fprintf(fid,'%.11f +% .11fj',[real(preamble(k)),imag(preamble(k))]);
        else
            fprintf(fid,'%.11f -% .11fj',[real(preamble(k)),-imag(preamble(k))]);
        end
    end
end
fclose(fid);

fid = fopen('zc_freq.txt','w');
for k = 1:length(fft_preamble)
%         keyboard
    if k < length(fft_preamble)
        if imag(fft_preamble(k)) >= 0
            fprintf(fid,'%.11f +% .11fj,\t',[real(fft_preamble(k)),imag(fft_preamble(k))]);
        else
            fprintf(fid,'%.11f -% .11fj,\t',[real(fft_preamble(k)),-imag(fft_preamble(k))]);
        end
    else
        if imag(fft_preamble(k)) >= 0
            fprintf(fid,'%.11f +% .11fj',[real(fft_preamble(k)),imag(fft_preamble(k))]);
        else
            fprintf(fid,'%.11f -% .11fj',[real(fft_preamble(k)),-imag(fft_preamble(k))]);
        end
    end
end
fclose(fid);


fileID = fopen('zc_time_int16.dat');
A = fread(fileID,[2,1024],'uint16');
figure
hold on
plot(A(1,:),A(2,:),'-oc','linewidth',3)
plot(preamble)





