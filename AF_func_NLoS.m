function [NLoS_att, T_NLoS] = AF_func_NLoS(c, lambda, NLoS_dist, cos_theta, er, p_BS, p_UE, temp_basis)

T_NLoS = NLoS_dist/c;
PL_NLoS = (lambda/(4*pi*NLoS_dist))^2;

[G_perp, G_parr] = AF_func_Fresnel(cos_theta, er);

% calculation of perp and parr components of p_UE
pUE_perp = temp_basis(2,:)*p_UE';
pUE_parr = norm(p_UE - pUE_perp*temp_basis(2,:)); 
% polarization change
NLoS_pol = G_perp * pUE_perp * temp_basis(2,:) + G_parr * pUE_parr * temp_basis(3,:);
G_NLoS = abs(NLoS_pol*p_BS');

NLoS_att = G_NLoS*PL_NLoS;
