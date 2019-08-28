clear all
close all

ref_sig = [];
rec_sig = [];

for l = 1:8
    
    name = ['measurements/reference-',num2str(l),'.txt'];
    fid=fopen(name,'r');
    
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        ref_sig = [ref_sig; str2num(tline)];
    end
    fclose(fid);
    
    name = ['measurements/measured-',num2str(l),'.txt'];
    fid=fopen(name,'r');
    
    fid_step = 0;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        rec_sig = [rec_sig; str2num(tline)];
        fid_step = fid_step + 1;
        if mod(fid_step,10000) == 0
            disp(num2str(fid_step))
        end
    end
    fclose(fid);
end

save('experiment5.mat', 'ref_sig', 'rec_sig')

% fileID = fopen('original.dat','r');
% A=fread(fileID,[9,2],'float32');
% fclose(fileID);