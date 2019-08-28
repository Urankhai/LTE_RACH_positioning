clear all
close all

CASE_FACTOR = 1; % if 1 then multipath, in another case single path.

phys_data
% Human penetration is human_pntr

ant_vect = 32;%[8, 16];%, 24, 32, 48, 64];
std_diff = zeros(size(ant_vect));
std_angle_diff = zeros(size(ant_vect));
PTx = 100; %mW
SNR14_dBm = 10*log10(PTx*(lambda/(4*pi*14000))^2); % dB
SNR_global = [];

colors_main = {'k-o', 'k--o', 'k-.o', 'k:o', 'k-'};
lwidth = [3,3,3,3,1];

for scenario_ind = 1:1
    for ant_step = 1:length(ant_vect)
        Ant_Num = ant_vect(ant_step);
        UE_Num = 1;     % in this scenarion we consider up to three UEs
        
        
        N_iter = 10;
        r_estm = zeros(1,N_iter);
        r_real = zeros(1,N_iter);
        x_diff = zeros(1,N_iter);
        y_diff = zeros(1,N_iter);
        r_diff = zeros(1,N_iter);
        r_orgn = zeros(UE_Num,N_iter);
        
        angle_diff = zeros(1,N_iter);
        angle_orgn = zeros(UE_Num,N_iter);
        
%         cond_KB = zeros(1,N_iter);
        for iii = 1:N_iter
%             close all
            
            Scenario_plain        % reflection planes, coordinates of UE and BS are defined in scenario
            r_orgn(:,iii) = radius;
            angle_orgn(:,iii) = angles;
            
            lte_tx          % signals "s" and rach preambles "rach" come from this script
            % upsampling in order to simulate analog signal
            up = 1; % the carrier's wavelength ~= 0.976m
            udt = dt/up;
            % UE.ups = kron(UE.s,ones(1,up));
            UE.uprach = kron(UE.rach,ones(1,up));
            
            % duration of rach signal after upconverter is 30720*up samples
            samples = 0:30720*up-1;
            carrier = exp(1j * w * udt * samples);
            tx_sign = zeros(size(UE.uprach));
            
            for k = 1:UE_Num
                tx_sign(k,:) = carrier.*UE.uprach(k,:);
            end
            
            Antenna         % from this point we know coordinates of the antenna elements
            
            
            rach_ant_user = zeros(Ant_Num, UE_Num, length(carrier));
            rach_ant = zeros(Ant_Num, length(carrier));
            rtd = zeros(1,Ant_Num);
            dsh = zeros(1,Ant_Num);
            nsh = zeros(1,Ant_Num);
            esh = zeros(1,Ant_Num);
            for ant_i = 1:Ant_Num
                Rx = BS.ant(ant_i,:);       % coordinates of the antenna element
                Rp = BS.pol;            % polarization of the antenna element
                
                temp_ant = 0;
                
                for k = 1:UE_Num
                    Tx = UE.loc(k,:);   % coordinates of UE
                    Tp = UE.pol(k,:);   % polarization of UE
                    
                    rach = tx_sign(k,:);
                    
                    [LoS_att, T_LoS] = AF_func_LoS([Rx, 0], Rp , [Tx, 0], Tp, c, lambda);
                    [LoS_UL_r, discret_shift, nondcrt_part, Nondscrt_shift] = AF_func_MIMO_rec(w, rach, udt, T_LoS);
                    %                 keyboard
                    %         rec_rach = LoS_att * LoS_UL_r;
                    dist_temp = norm(Rx - Tx);
                    SNR_dBm = 10*log10(PTx*(lambda/(4*pi*dist_temp))^2); % dB
                    SNR = SNR_dBm - SNR14_dBm;
                    SNR_global = [SNR_global, SNR];
                    
                    sigma = std(LoS_UL_r)/sqrt(2)/(10^(SNR/20));
                    noise = sigma*(randn(size(LoS_UL_r)) + 1i*randn(size(LoS_UL_r)));
                    
                    rec_rach_noise = LoS_UL_r + noise;
                    rec_rach_fft = fft(rec_rach_noise);
                    rec_rach_fft_clear = zeros(size(rec_rach_fft));
                    dw_freq = 19526;%20160;
                    up_freq = 20590;%21250;
                    rec_rach_fft_clear(dw_freq:up_freq) = rec_rach_fft(dw_freq:up_freq);
                    
                    rec_rach = ifft(rec_rach_fft_clear);
                    
                    rach_ant_user(ant_i, k, :) = rec_rach.*conj(carrier);
                    temp_ant = temp_ant + rec_rach;
                    
                    %                 disp(['T_LoS = ', num2str(ceil(T_LoS/udt))])
                    
                    if k == 1
                        rtd(ant_i) = T_LoS;
                        dsh(ant_i) = discret_shift;
                        nsh(ant_i) = nondcrt_part;
                        esh(ant_i) = Nondscrt_shift;
                    end
                end
                
                rach_ant(ant_i, :) = temp_ant.*conj(carrier);
                
            end
            asd = reshape(rach_ant_user(:,1,:),[],30720);
            
            nnz_u = nnz(asd(1,1:500));
            nnz_d = nnz(asd(end,1:500));
            delta_phase = 0.110446616721923;
            
            %         if  nnz_u ~= nnz_d
            %             disp('layout detected')
            %             keyboard
            %
            % %             nnz_ant = asd(:,nnz_u);
            % %             layout_ant =
            %
            %         end
            
            rach_detection
            
            
            r_estm(iii) = sqrt((BS.loc(1) - UE_location(1))^2 + (BS.loc(2) - UE_location(2))^2);
            r_real(iii) = radius(1);
            r_diff(iii) = r_estm(iii) - r_real(iii);
            
            
            angle_real = 180*acos([1, 0]*(( UE.loc(UEi,:) - BS.loc )/norm( UE.loc(UEi,:) - BS.loc ))')/pi;
            angle_estm = 180*acos([1, 0]*(( UE_location - BS.loc )/norm( UE_location - BS.loc ))')/pi;
            angle_diff(iii) = (angle_real - angle_estm);
            
            
            
            %         disp(['x diff = ', num2str(x_diff(iii)), '; y diff = ', num2str(y_diff(iii))])
            disp(['Radius error = ', num2str(r_diff(iii)),'; angle error = ', num2str(angle_diff(iii))])
            
            
                    x_diff(iii) = UE.loc(UEi, 1) - UE_location(1);
                    y_diff(iii) = UE.loc(UEi, 2) - UE_location(2);
                    r_diff(iii) = sqrt(x_diff(iii)^2 + y_diff(iii)^2);
%             cond_KB(iii) = log10(cond(KB));
            
            if abs(r_diff(iii)) > 20
                disp('error is too big')
                %             keyboard
            end
            
            %         close Figure 101 Figure 102
            %         keyboard
            if c1^2 - 4*c2*c0 < 0
                disp('wrong equations')
                keyboard
                iii = iii - 1;
            end
            %         close all
        end
        
        radang_temp = [r_diff;angle_diff];
        
        std_r_temp = std(r_diff);
        num_exceed = nnz(abs(r_diff - mean(r_diff)) > 4*std_r_temp);
%         figure
%         hold on
        while num_exceed >= 1
            %         plot(radang_temp(1,:))
            %         keyboard
            [~, pos_temp] = max(radang_temp(1,:));
            radang_temp = [radang_temp(:, 1:pos_temp-1),radang_temp(:, pos_temp+1:end)];
            std_r_temp = std(radang_temp(1,:));
            
            num_exceed = nnz(abs(radang_temp(1,:) - mean(radang_temp(1,:))) > 3.5*std_r_temp); % 3.5 sigma = 99.95
        end
        
        disp(['Iterations = ',num2str(N_iter), '; number left = ', num2str(length(radang_temp))])
        r_stat = radang_temp(1,:);
        angle_stat = radang_temp(2,:);
        
%         figure(1000+ant_step)
%         hold on
%         title('coordinates error')
%         plot(r_stat,'-o')
        
%         figure(2000+ant_step)
%         hold on
%         title('angle error')
%         plot(angle_stat,'-o')
        
        %     figure(3000+ant_step)
        %     hold on
        %     title('normalizaed errors')
        %     plot(r_diff/min(r_diff),'k-o')
        %     plot(angle_diff/max(angle_diff),'r-*')
        %     legend('distance', 'angle')
        
        std_diff(ant_step) = std(r_stat);%(sum(r_diff)/N_iter);
        std_angle_diff(ant_step) = std(angle_stat);%(sum(angle_diff)/N_iter);
        
        %     disp([num2str(N_iter),' iterations: average r_diff = ', num2str(avg_diff(ant_step))])
        %     disp([num2str(N_iter),' iterations: average angle error = ', num2str(avg_angle_diff(ant_step))])
        
        
        disp(['mean distance error = ', num2str(mean(r_stat)), '; std of distance error = ', num2str(std(r_stat))])
        disp(['   mean angle error = ', num2str(mean(angle_stat)), ';    std of angle error = ', num2str(std(angle_stat))])
        
        %     figure
        %     hold on
        %     title('Relation matrix condition with coordinate errors')
        %     grid on
        %     bar(cond_KB, log10(abs(r_stat)))
        
%         figure
%         hist(r_stat)
% %             
        
        
    end
%     keyboard
    
    
    figure(5000)
    hold on
    title('Average error of coordinates estimation')
    plot(ant_vect, std_diff,colors_main{scenario_ind},'linewidth',lwidth(scenario_ind))


    
    figure(6000)
    hold on
    title('Average error of angle etimation')
    plot(ant_vect, std_angle_diff, colors_main{scenario_ind}, 'linewidth',lwidth(scenario_ind))

    
end
figure(5000)
legend('30-100','100-200','200-300','300-500','500-1000')
figure(6000)
legend('30-100','100-200','200-300','300-500','500-1000')

figure
hold on
title('SNR in dBm')
plot(SNR_global)