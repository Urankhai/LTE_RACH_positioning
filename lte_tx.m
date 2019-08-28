% LTE system parameters
w_LTE = 30.72e3;                        % the number of chips in 1ms
dt = 1/w_LTE/1000;                      % the duration of one chip

UE_Num = 1;
% Duration of PRACH is 1ms for a conventional format (R <= 14km)
UE.rach = zeros(UE_Num, w_LTE);
CZ_seq = zeros(839,UE_Num);
for i = 1:UE_Num
    [rach, seq] = AF_func_PRACH(1); % since UE_Num is the number of reflectors
    UE.rach(i,:) = transpose(rach);
    CZ_seq(:,i) = seq;
    
%     clear seq rach
end

% 
% RB_max = 100;%169;% is the maximum number of RBs
% N_sub = 12;
% N_frm = 1;
% N_sbfrm = N_frm*2;
% frm_elm = 12;
% % N_smbls = N_frm*frm_elm + N_sbfrm;% number of RS symbols is equal to the number of subframes
% N_guard = 5;
% 
% Nfft = 2048;
% 
% set_columns = 1:N_frm*(frm_elm + 2);
% set_dmrs = 4:7:N_frm*(frm_elm + 2);
% 
% set_inf = setdiff(set_columns,set_dmrs);
% 
% % Modulation and Demodulation setup
% M = 16;                         % The order of Modulation constellation
% hMod = modem.qammod(M);         % Create a 16-QAM modulator
% hMod.InputType = 'Bit';         % Accept bits as inputs
% hMod.SymbolOrder = 'Gray';      % Accept bits as inputs
% hDemod = modem.qamdemod(hMod);  % Create a 16-QAM based on the modulator
% 
% 
% 
% % s = AF_func_signal_gen(PUSCH);
% if UE_Num == 3
%     RBs = [30,30,40];
% elseif UE_Num == 2
%     RBs = [45, 45];
% else
%     RBs = 100;
% end
% 
% if sum(RBs) > 100
%     disp('choose another resource allocation RBs')
%     return
% end
% 
% start_point = N_guard;
% length_GT = 2976; % the same Guard Time as for PRACH
% UE.s = zeros(UE_Num, N_frm*(w_LTE + length_GT));  % Generated signals
% UE.RL = zeros(UE_Num,2);            % Resource Location
% for i = 1:UE_Num
%     
%     
%     [INF, INF_MOD, PUSCH, RS] = AF_func_UL_signal(RBs(i), N_frm, hMod);
%     N_carri = size(PUSCH,1);
%     N_smbls = size(PUSCH,2);
%     
%     fft_array = zeros(Nfft, N_smbls);
%     finish_point = start_point + N_carri;
%     fft_array(start_point+1: finish_point,:) = PUSCH;
%     
%     % signal generation for multiuser mode
%     s = AF_func_signal_gen_mult(N_smbls, fft_array);
%     
%     UE.RL(i,:) = [start_point+1, finish_point];
%     start_point = finish_point + N_guard;  
%     UE.s(i,:) = s;
%     
%     clear s
% end
