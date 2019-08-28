clear all
close all
Nsamples = 4000;


% name = 'measurements/reference1.txt';
% fid=fopen(name,'r');
% ref_sig = [];
% while 1
%     tline = fgetl(fid);
%     if ~ischar(tline), break, end
%     ref_sig = [ref_sig; str2num(tline)];
% end
% fclose(fid);
% 
% name = 'measurements/measured1.txt';
% fid=fopen(name,'r');
% rec_sig = [];
% while 1
%     tline = fgetl(fid);
%     if ~ischar(tline), break, end
%     rec_sig = [rec_sig; str2num(tline)];
% end
% fclose(fid);
% save('experiment1.mat', 'ref_sig', 'rec_sig')


load('experiment5.mat')

figure
hold on
plot(ref_sig(:,1),'-o')
plot(ref_sig(:,2),'-*r')
plot(rec_sig(:,1),'g')
plot(rec_sig(:,2),'m')
legend('real','imag','real','imag')



complex_ref = ref_sig(:,1) + 1i*ref_sig(:,2);
complex_rec = rec_sig(:,1) + 1i*rec_sig(:,2);

colors = {'b','g','r','c','m','k',':b',':g',':r',':c',':m',':k','-ob','-og','-or','-oc','-om','-ok'};
figure
hold on
plot(complex_ref(1:10),'linewidth',3)
plot(complex_rec(1:10),'g','linewidth',3)
plot(complex_ref(10:20))
plot(complex_rec(10:20),'g')

% 
% norm_ref = complex_ref./abs(complex_ref);
% norm_rec = complex_rec./abs(complex_rec);


dt = 1/1e6; % sampling rate of the USRPs
w = 50e3;
t = dt*(1:16*Nsamples);

Ant_Num = 8;
e_phi = zeros(1,Ant_Num);
est_error = zeros(1,Ant_Num);
dot_prod = zeros(1,Ant_Num);
phases = zeros(1,Ant_Num);

for i = 1:Ant_Num
    period = 1+(i-1)*Nsamples:i*Nsamples;
    
    ref = complex_ref(period);
    rec = complex_rec(period);
    
    % nonzero elements
    nnz_measurement = find(abs(rec));
    norm_ref = ref(nnz_measurement)./abs(ref(nnz_measurement));
    norm_rec = rec(nnz_measurement)./abs(rec(nnz_measurement));
    
    t_period = dt*(nnz_measurement-1);
    base_func = exp(1j*2*pi*w*t_period);
    
    z = norm_rec;
    e_phi(i) = (1/Nsamples)*(base_func'*z)/abs((1/Nsamples)*(base_func'*z));
    est_error(i) = norm(z - base_func*e_phi(i));
    
    
    figure 
    hold on
    plot(real(norm_ref(1:200)),'c','linewidth',3)
    plot(real(base_func(1:200)),'k')
    
    plot(real(norm_rec(1:200)),'g','linewidth',3)
    plot(real(base_func(1:200)*e_phi(i)),'r')
    corrected_rec = base_func*e_phi(i);
    
    dot_prod(i) = norm_ref'*corrected_rec/norm(norm_ref)/norm(corrected_rec);
    phases(i) = phase(dot_prod(i));
    
    
    
end
figure
plot(upsample(dot_prod,2))
figure
plot(phases)




% figure
% hold on
% plot(norm_ref(1:10),'linewidth',3)
% plot(norm_rec(1:10),'g','linewidth',3)
% plot(norm_ref(10:20))
% plot(norm_rec(10:20),'g')
% 
% figure
% hold on
% plot(t,real(norm_ref),'-o')
% plot(t,imag(norm_ref),'-*r')
% plot(t,real(norm_rec),'g')
% plot(t,imag(norm_rec),'m')
