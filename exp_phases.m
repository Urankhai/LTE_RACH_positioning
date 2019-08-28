clear all
close all
Nsamples = 4000;

load('experiment5.mat')

figure
hold on
plot(ref_sig(:,1),'-o')
% plot(ref_sig(:,2),'-*r')
plot(rec_sig(:,1),'g')
% plot(rec_sig(:,2),'m')
legend('reference','measured')



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

nnz_num = zeros(1,Ant_Num);

for i = 1:Ant_Num
    period = 1+(i-1)*Nsamples:i*Nsamples;
    
    ref = complex_ref(period);
    rec = complex_rec(period);
    
    % nonzero elements
    nnz_measurement = find(abs(rec));
    
    nnz_num(i) = length(nnz_measurement);
    
    norm_ref = ref(nnz_measurement)./abs(ref(nnz_measurement));
    norm_rec = rec(nnz_measurement)./abs(rec(nnz_measurement));
    
    t_period = dt*(nnz_measurement-1);
    base_func = exp(1j*2*pi*w*t_period);
    
    
    
    paths_num = 4;
    % z = H x + r; x = (H'*R^(-1)*H)^(-1)*H'*R^(-1)*z
%     H = base_func*ones(1, paths_num);
    z = norm_rec;
    
    sss = 0.5;
    R = sss*ones(length(z))+(1-sss)*eye(length(z));
    est = (1/Nsamples)*(base_func'*z);%(base_func'*R^(-1)*base_func)*(base_func'*R^(-1)*z);
    e_phi(i) = est/abs(est);
    est_error(i) = norm(z - base_func*e_phi(i));
    
    
    figure 
    hold on
    plot(real(norm_ref(1:200)),'c','linewidth',3)
    plot(real(base_func(1:200)),'k')
    
    plot(real(norm_rec(1:200)),'g','linewidth',3)
    plot(real(base_func(1:200)*e_phi(i)),'r')
    corrected_rec = base_func*e_phi(i);
    
    dot_prod(i) = norm_ref'*corrected_rec;%/norm(norm_ref)/norm(corrected_rec);
    phases(i) = angle(dot_prod(i));
    
    
    
end
figure
plot(upsample(e_phi,2))

figure
plot(upsample(dot_prod,2))
figure
plot(phases)

number_sycle = 1;
PPP = zeros(15, Ant_Num);
PPP(1,1) = phases(1);
steps = 1;
for k = 1:Ant_Num-1
    steps = steps + 1;
    if abs( phases(k) - phases(k+1) ) > pi
        number_sycle = number_sycle + 1;
        PPP(number_sycle, steps) = phases(k+1);
%         steps = 0;
    else
        PPP(number_sycle, steps) = phases(k+1);
    end
end
gaps_number = nnz(sum(PPP,2));

temp_phases = zeros(1, Ant_Num);
gap = 0;
shift = length(nonzeros(PPP(1,:)));
temp_phases(1:shift) = nonzeros(PPP(1,:));
for i = 1:gaps_number-1
    temp_u = nonzeros(PPP(i,:));
    temp_d = nonzeros(PPP(i+1,:));
    
    if temp_u(end) > temp_d(1)
        gap = gap + 1;
    else
        gap = gap - 1;
    end
    
    
    temp_phases(shift+1: shift + length(temp_d)) = temp_d + gap*2*pi;
    shift = shift + length(temp_d);
end


corrected_phases = fliplr(-temp_phases - min(-temp_phases));
steps_bar = 1600+(0:7)*4000;
value_bar = -6*10^(-4)*ones(1,8);
figure
subplot(2,1,1), plot(rec_sig(:,1),'g')
hold on
subplot(2,1,1), bar(steps_bar, value_bar)
subplot(2,1,2), plot(corrected_phases,'-o')

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
