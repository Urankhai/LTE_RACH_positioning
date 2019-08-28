clear all
close all

CASE_FACTOR = 1; % if 1 then multipath, in another case single path.

phys_data
% Human penetration is human_pntr

UE_Num = 3;     % in this scenarion we consider up to three UEs
Scenario        % reflection planes, coordinates of UE and BS are defined in scenario
lte_tx          % signals "s" and rach preambles "rach" come from this script

% upsampling in order to simulate analog signal
up = 10; % the carrier's wavelength ~= 0.976m
udt = dt/up;
UE.ups = kron(UE.s,ones(1,up));
UE.uprach = kron(UE.rach,ones(1,up));

step_grid = 1:size(UE.uprach,2);
time_grid = step_grid*udt;
carrier = exp(1i * w * time_grid);


Antenna         % from this point we know coordinates of the antenna elements


rach_signal = carrier.*UE.uprach(1,:);