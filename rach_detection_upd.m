% close all

colors = {'b-o','r-*','g-d'};
% for users = 1
phases = zeros(1, Ant_Num);
corr_max = zeros(1, Ant_Num);
% SNR = 10;

for UEi = 1
    
    temp_max = zeros(1,Ant_Num);
    
    
    ampl_frwd = 0;
    ampl_cntr = 0;
    ampl_bcwd = 0;
    
    phases_frwd = zeros(1, Ant_Num);
    phases_cntr = zeros(1, Ant_Num);
    phases_bcwd = zeros(1, Ant_Num);
    
    distances = zeros(Ant_Num,2);
    
    for k = 1:Ant_Num
        
        %         temp_rach1 = reshape(rach_ant_user(k,UEi,:),1,[]);
        temp_rach1 = rach_ant(k,:);
        
        signal_power = norm(temp_rach1);
        %         temp_rach1 = temp_rach1 + noise;
        %         temp_rach1 = rach_ant(k,:);
        s1 = temp_rach1(1:up:end);
        [out_s1, ~] = AF_func_correlator(CZ_seq(:,UEi), s1);
        
        
        % upconverted check
        rach_org = UE.rach(1,:);
        test_seq = xcorr(s1,rach_org);
        [A1,I1] = max(test_seq);
        shift_pos = I1 - 30720;
        
        r_left = (shift_pos-1)*udt*c;
        r_right = shift_pos*udt*c;
        distances(k,:) = [r_left, r_right];
        
        
        %         keyboard
        %                         figure(10+k)
        %                         plot(abs(out_s1(1:100)), '-o')
        %                         keyboard
        temp_array = out_s1(1:100);
        diff_array = abs(diff(abs(temp_array)));
        [~, diff_num] = max(diff_array);
        
        
        
        
        [max_out_s, num_of_max] = max(temp_array);
        temp_max(k) = num_of_max;
        
        if k == 1
            center_pos = num_of_max;
        end
        
        if center_pos == 1
            ampl_cntr = ampl_cntr + abs(temp_array(center_pos));
            phases_cntr(k) = angle(temp_array(center_pos));
            
            ampl_frwd = 0;
            phases_frwd(k) = 0;
            
            ampl_bcwd = ampl_bcwd + abs(temp_array(center_pos+1));
            phases_bcwd(k) = angle(temp_array(center_pos+1));
        else
            ampl_cntr = ampl_cntr + abs(temp_array(center_pos));
            phases_cntr(k) = angle(temp_array(center_pos));
            
            ampl_frwd = ampl_frwd + abs(temp_array(center_pos-1));
            phases_frwd(k) = angle(temp_array(center_pos-1));
            
            ampl_bcwd = ampl_bcwd + abs(temp_array(center_pos+1));
            phases_bcwd(k) = angle(temp_array(center_pos+1));
        end
        
    end
    
    [max_ampl, pos] = max([ampl_frwd, ampl_cntr, ampl_bcwd]);
    
    if pos == 1
        phases = phases_frwd;
    elseif pos == 2
        phases = phases_cntr;
    else
        phases = phases_bcwd;
    end
    
    if nnz(diff(temp_max)) > 0
        disp('position change')
        %         keyboard
    end
    %     phases = phases_bcwd;
end
phases = -phases;
phase_sort
% keyboard

diff_phases = diff(temp_phases);
dd_phases = diff(diff_phases);

num_layouts = 0;
pos_layouts = [];
sign_layouts = [];
corrected_phases = temp_phases;
for l = 1:length(dd_phases)
    if (abs(dd_phases(l)) > 0.1)
        %                         keyboard
        num_layouts = num_layouts + 1;
        
        
%         if num_layouts > 2
%             keyboard
%         end
        
        if l > 1
            pos_layouts = [pos_layouts,l+2];
            sign_layouts = [sign_layouts, -sign(dd_phases(l))];
        else % l = 1
            if (abs(dd_phases(l+1)) > 0.05)
                pos_layouts = [pos_layouts,l+2];
                sign_layouts = [sign_layouts, -sign(dd_phases(l))];
            else
                pos_layouts = [pos_layouts, l+1];
                sign_layouts = [sign_layouts, sign(dd_phases(l))];
            end
            
        end
        delta_phase = 0.110446616721923;
        %             delta_phase = 0.010847358072824;
        if mod(num_layouts,2) == 1
            for k = pos_layouts(num_layouts):length(temp_phases)
                corrected_phases(k) = temp_phases(k) + delta_phase*sign_layouts(num_layouts);
                %                 keyboard
            end
        end
    end
end
% keyboard
% figure(103)
% % hold on
% % title('Corrected phases')
% subplot(4,1,1), plot(corrected_phases,'-ro')
% hold on
% subplot(4,1,1), plot(zz,'-g*')
% subplot(4,1,2), plot(diff(corrected_phases),'-ro')
% subplot(4,1,3), plot(diff(diff(corrected_phases)),'-mo')
% 
% subplot(4,1,4), plot(corrected_phases - yy','-ro','linewidth',3)


% figure(102)
% % hold on
% % title('Corrected phases')
% subplot(4,1,1), plot(temp_phases,'-ro')
% 
% subplot(4,1,2), plot(diff(temp_phases),'-ro')
% subplot(4,1,3), plot(diff(diff(temp_phases)),'-mo')
% 
% subplot(4,1,4), plot(temp_phases - yy','-ro','linewidth',3)
% hold on
% subplot(5,1,4), plot(zeros(size(temp_phases)),'b')

% subplot(5,1,4), plot(diff(temp_phases - zz'),'-ro')
% subplot(5,1,4), plot(diff(temp_phases - zz'),'-ro')
% subplot(5,1,5), plot(diff(diff(temp_phases - zz')),'-ro')

% subplot(5,1,5), plot(corr_max)

HH = [(1:Ant_Num).^2', (1:Ant_Num)', ones(Ant_Num,1)];
abc = (HH'*HH)\HH'*corrected_phases';
correct_phase_smooth = HH*abc;

% asd = temp_phases;
% delete
% temp_phase = zeros(1,44);
test_surface
% temp_phases = new_ph;
corrected_phases = correct_phase_smooth';
correct_pase = corrected_phases - min(corrected_phases);
circle = c*correct_pase/w;

[A2,I2] = min(circle);

% keyboard
if circle(1) == 0
    HH = Ant_Dist*(0:Ant_Num-1)';
    % inv_HH = (HH'*HH)^(-1)*HH';
    inv_HH = HH'/(sum(HH.^2));
    cos_alpha = inv_HH*circle';
    est_alpha = acos(cos_alpha);
    
    
    est_X = BS.ant(:,1) + cos_alpha*circle';
    est_Y = BS.ant(:,2) - sin(est_alpha)*circle';
    
    line_X = est_X;
    line_Y = BS.ant(:,2) - (BS.ant(1,1) - line_X)/tan(est_alpha);
elseif circle(end) == 0
    HH = -Ant_Dist*(fliplr(0:Ant_Num-1))';
    % inv_HH = (HH'*HH)^(-1)*HH';
    inv_HH = HH'/(sum(HH.^2));
    cos_alpha = inv_HH*circle';
    est_alpha = acos(cos_alpha);
    
    
    est_X = BS.ant(:,1) + cos_alpha*circle';
    est_Y = BS.ant(:,2) - sin(est_alpha)*circle';
    
    line_X = est_X;
    line_Y = BS.ant(:,2) + (line_X - BS.ant(end,1))/tan(est_alpha);
end


% figure
% hold on
% title('difference between line and the circle surface')
% subplot(4,1,1); plot(est_X, est_Y, 'c-o','linewidth',3); hold on; plot(est_X, line_Y, '-*')
% subplot(4,1,2); plot(est_X, est_Y - line_Y,'-o')
% subplot(4,1,3); plot(diff(est_Y)./diff(est_X),'-*')
% subplot(4,1,4); plot(diff(diff(est_Y))./diff(diff(est_X)),'-*')

% figure(100+Ant_Num)
% hold on
% plot(est_X,est_Y,'g-d','linewidth',2)
% plot(line_X, line_Y, 'k-p')


angle_real = 180*acos([1, 0]*(( UE.loc(UEi,:) - BS.loc )/norm( UE.loc(UEi,:) - BS.loc ))')/pi;



% ttt = 0:0.01:1;
% for rrr = 1:Ant_Num
%     x_small = circle(rrr)*cos(2*pi*ttt) + BS.ant(rrr,1);
%     y_small = circle(rrr)*sin(2*pi*ttt) + BS.ant(rrr,2);
%     plot(x_small,y_small)
% end


% plot(BS.ant(:,1), BS.ant(:,2) - circle', 'c-*', 'linewidth',3)
% figure(99)
% hold on
% title('Circle from phases')
% plot(ant_coord(:,1), circle,'-co')
% keyboard

% Bancroft's algorithm
collection_point = [250 299.9];
local_coordinates = BS.ant - ones(Ant_Num,1)*collection_point;

% The noise matrix W
W = ones(Ant_Num)*0.125+eye(Ant_Num)-0.125*eye(Ant_Num);%lambda/2/pi*(0.04)^2*eye(Ant_Num);


% good = [1, 4, 6, 8];

B = [local_coordinates, (circle)'];
KB = B'*W^(-1)*B;
% disp(['log10 of condition number of KB = ', num2str(log10(cond(KB)))])
invB = (KB)^(-1)*B'*W^(-1);
M = [1 0 0 ; 0 1 0 ; 0 0 -1];
I = ones(Ant_Num, 1);
a = zeros(Ant_Num, 1);

for k = 1:Ant_Num
    a(k) = 0.5* B(k,:)*M*B(k,:)';
end

% c2 * Lambda^2 + c1 * Lambda + c0 = 0
c0 = (invB*a)' * M * invB*a;
c1 = 2*((invB*I)' * M * invB*a - 1);
c2 = (invB*I)' * M * invB*I;


Lambda1 = (-c1 + sqrt(c1^2 - 4*c2*c0))/(2*c2);
solution1 = M * invB * (a + Lambda1*I);


Lambda2 = (-c1 - sqrt(c1^2 - 4*c2*c0))/(2*c2);
solution2 = M * invB * (a + Lambda2*I);

if solution1(2) + collection_point(2) > BS.loc(2)
    UE_location = solution2(1:2)' + collection_point;
    UE_radius = -solution2(3);
else
    UE_location = solution1(1:2)' + collection_point;
    UE_radius = -solution1(3);
end

theta_estm = acos([1, 0]*(( UE_location - BS.loc )/norm( UE_location - BS.loc ))');
plecho = ant_coord(1,1) - BS.loc(1,1);

if (UE_radius > distances(I2,1)) && (UE_radius < distances(I2,2))
    disp('estimation is probably correct')
else
    disp('estimation is incorrect')
%     keyboard

    r_mean = (distances(I2,1) + distances(I2,2))/2;
    
    if theta_estm <= pi/2
%         keyboard
        r_bs = plecho*sin(pi/2 - theta_estm);
        r_sum = r_bs + r_mean;
        pUE_mean = [BS.loc(1,1) + r_sum*sin( pi/2 - theta_estm ),BS.loc(1,2) - r_sum*cos( pi/2 - theta_estm )];
    else
%         keyboard
        r_bs = plecho*sin(theta_estm - pi/2);
        r_sum = r_bs + r_mean;
        pUE_mean = [BS.loc(1,1) + r_sum*sin( pi/2 - theta_estm ),BS.loc(1,2) - r_sum*cos( pi/2 - theta_estm )];
    end
    UE_location = pUE_mean;
    UE_radius = norm(pUE_mean - BS.loc);
end



% if UE_radius

% keyboard



% if UE_radius <= 0
%     disp('need to check')
%     keyboard
% end


disp(['     Real UE position = [', num2str(UE.loc(UEi, 1), '%10.2f'),', ', num2str(UE.loc(UEi, 2), '%10.2f'),'];'])
disp(['Estimated UE location = [', num2str(UE_location(1), '%10.2f'),', ', num2str(UE_location(2), '%10.2f'), ']; Rasius of the circle = ', num2str(UE_radius)])

 figure(100+Ant_Num)
 plot([BS.loc(1), UE_location(1)],[BS.loc(2), UE_location(2)],'-s', 'MarkerEdgeColor','k', 'MarkerFaceColor','c','MarkerSize',6)

nonlin_search


