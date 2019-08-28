% close all
% phase sorting
% figure(101)
% hold on
% title('Original phases')
% plot(-phases,'-o')


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
% plot(temp_phases,'-*')

% LS
z = temp_phases';
H = [(1:Ant_Num)', ones(Ant_Num,1)];
x = (H'*H)^(-1)*H'*z;
zz = H*x;

HH = [(1:Ant_Num).^2', (1:Ant_Num)', ones(Ant_Num,1)];
y = (HH'*HH)\HH'*z;
yy = HH*y;

% test_surface
% temp_phases = new_ph;

% keyboard

% figure(102)
% % hold on 
% % title('Corrected phases')
% subplot(5,1,1), plot(temp_phases,'-ro') 
% hold on 
% subplot(5,1,1), plot(zz,'-g*') 
% 
% subplot(5,1,2), plot(diff(temp_phases),'-ro')
% subplot(5,1,3), plot(diff(diff(temp_phases)),'-mo')
% 
% subplot(5,1,4), plot(temp_phases - zz','-ro','linewidth',3)
% % hold on
% % subplot(5,1,4), plot(zeros(size(temp_phases)),'b')
% 
% % subplot(5,1,4), plot(diff(temp_phases - zz'),'-ro')
% % subplot(5,1,4), plot(diff(temp_phases - zz'),'-ro')
% % subplot(5,1,5), plot(diff(diff(temp_phases - zz')),'-ro')
% 
% subplot(5,1,5), plot(corr_max)



