clear all
temp_ID = 1;
[rach, seq] = AF_func_PRACH(temp_ID);

fid = fopen('ZCseq2.txt','w');
for k = 1:839%100:105
    if imag(seq(k)) >= 0
        asd = fprintf(fid, '%1.15f +%1.15f\r\n', real(seq(k)), imag(seq(k)));
    else
        asd = fprintf(fid, '%1.15f -%1.15f\r\n', real(seq(k)), -imag(seq(k)));
    end
end


fclose(fid);