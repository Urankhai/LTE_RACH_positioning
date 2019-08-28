function [INF, INF_MOD, PUSCH, ZC] = AF_func_UL_signal(Ny_RB, Nx_frm, hMod)

% time domain
Nx_sbfrm = Nx_frm*2; % number of subframes
Nx_smbls = 12; % number of symbols per frame
Nx = Nx_frm*Nx_smbls; % Number of OFDM symbols
% frequency domain
Ny_sub = 12; % number of subcarriers per RB
Ny = Ny_RB*Ny_sub;    % Number of all subcarriers


N_elm =  Nx * Ny;     % number of elements in Radio Recourse

% generation of Zadoff Chu sequencies
ZC = zeros(Ny_RB*Ny_sub, Nx_sbfrm);
root = [3,5,11,17,23,31,41];
for k = 1:Nx_sbfrm
%     keyboard
    asd = lteZadoffChuSeq(root(mod(k - 1,length(root)) + 1),Ny_RB*Ny_sub+1);
    ZC(:,k) = asd(1:Ny_RB*Ny_sub);
end


INF = randi([0 1],N_elm*sqrt(hMod.M),1);         % bit stream
% INF_DEC = bi2de(reshape(INF,[],log2(hMod.M)));   % bit to decimal convertion
% INF_MOD = qammod(INF_DEC, hMod.M);%, 0,'gray');               % modulation
INF_MOD = modulate(hMod,INF);
% keyboard
% INF_MOD = INF_MOD;%/mean(abs(INF_MOD));


INF_MAP = zeros(Ny, Nx);

for k = 1:Nx_frm
    for l = 1:Ny_RB
        temp_positions = ((k-1)*Ny_RB+(l-1))*Ny_sub*Nx_smbls+1 : ((k-1)*Ny_RB+l)*Ny_sub*Nx_smbls;
        INF_MAP((l-1)*Ny_sub+1 : l*Ny_sub, (k-1)*Nx_smbls+1 : k*Nx_smbls) = reshape(INF_MOD(temp_positions), Ny_sub, Nx_smbls);
        %reshape(INF_MOD((l-1)*Ny_sub*Nx_smbls + (k-1)*Ny_RB*Ny_sub*Nx_smbls + 1:l*Ny_sub*Nx_smbls + (k-1)*Ny_RB*Ny_sub*Nx_smbls), Ny_sub, Nx_smbls);
    end
end

set_columns = 1:Nx_frm*(Nx_smbls + 2);
set_dmrs = 4:7:Nx_frm*(Nx_smbls + 2);
set_inf = setdiff(set_columns,set_dmrs);
PUSCH = zeros(Ny_RB*Ny_sub, Nx_frm*(Nx_smbls + 2));

PUSCH(:,set_dmrs) = ZC;
PUSCH(:,set_inf) = fft(INF_MAP, Ny_RB*Ny_sub)/sqrt(Ny);