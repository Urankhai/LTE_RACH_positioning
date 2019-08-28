function [LoS_att, T_LoS] = AF_func_LoS(BS, p_BS , UE, p_UE, c, lambda)

% Path loss (simple approach)
LoS_dist = norm(BS - UE);
if LoS_dist < 1
    PL_LoS = (lambda/(4*pi))^2;
else
    PL_LoS = (lambda/(4*pi*LoS_dist))^2;
end
% Polarization change
if LoS_dist == 0
    LoS_vect = [1 0 0];
else
    LoS_vect = (BS - UE)/LoS_dist;
end
[~, LoS_pol] = AF_func_polarization(p_UE, LoS_vect);
G_LoS = abs(LoS_pol*p_BS');

% Output data
T_LoS = LoS_dist/c;
LoS_att = G_LoS*PL_LoS;