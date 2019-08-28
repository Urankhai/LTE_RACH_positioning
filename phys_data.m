% Global variables
global c F w lambda human_pntr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Normalization = 1;
K_grad = 180/pi; % Scaling coefficient
c = 3e8;
F = 2.6e9;c*10;
GHz = F/3;
lambda = c/F;
e0 = 1/(4*pi * 1e-7)/c^2; %Vacuum permittivity
k0 = 2*pi/lambda;
w = 2*pi*F;
%--------------------------------------------------------------------------
% Coated glass
%--------------------------------------------------------------------------
glass_a = 7;
glass_b = 0;
glass_c = 0.25;
glass_d = 0;

glass_e1 = glass_a*(F/GHz)^glass_b;
glass_sigma = glass_c*(F/GHz)^glass_d;
glass_err = 17.98*glass_sigma/(F/GHz);
glass_er = glass_e1 - 1j*glass_err;
%--------------------------------------------------------------------------
% Concrette
%--------------------------------------------------------------------------
concrette_a = 5.31;
concrette_b = 0;
concrette_c = 0.0326;
concrette_d = 0.8095;

concrette_e1 = concrette_a*(F/GHz)^concrette_b;
concrette_sigma = concrette_c*(F/GHz)^concrette_d;
concrette_err = 17.98*concrette_sigma/(F/GHz);
concrette_er = concrette_e1 - 1j*concrette_err;

%--------------------------------------------------------------------------
% Human body
%--------------------------------------------------------------------------
human_thickness = 0.6; % the average human's thickness is equal to 60 cm
human_a = 2.94;
human_b = 0;
human_c = 0.0116;
human_d = 0.7076;

human_e1 = human_a*(F/GHz)^human_b;
human_sigma = human_c*(F/GHz)^human_d;
human_err = 17.98*human_sigma/(F/GHz);
human_er = human_e1 - 1j*human_err;

% we will supose that human body is penetrated by signal, thus, the angle
% of incident is equal to 90 degree
human_appe = sqrt(human_er - 1);
human_perp = (1 - human_appe)/(1 + human_appe);
human_parr = (human_er - human_appe)/(human_er + human_appe);
human_delt = k0*sqrt(human_er - 1);
T_perp = ((1-human_perp^2)*exp(-1j*(human_delt - k0*human_thickness)))/(1 - human_perp^2* exp(-1j*2*human_delt));
T_parr = ((1-human_parr^2)*exp(-1j*(human_delt - k0*human_thickness)))/(1 - human_parr^2* exp(-1j*2*human_delt));

human_pntr = sqrt(T_perp*T_perp' + T_parr*T_parr');
% %--------------------------------------------------------------------------
% % Concrete
% %--------------------------------------------------------------------------
% concrete_perp = 0.6; % according to METIS for 2.6 GHz
% concrete_parr = 0.3; % according to METIS for 2.6 GHz