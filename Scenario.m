[SS, Rx, Oa, sigma] = AF_data_KIOTorres();
planes = SS(1:3,:);%[[0,0,1,0]; SS(1:2,:)]; % add flat ground


rx_p = [0, 0, 1];                     % Polarization vector of BS
rx_p = rx_p/norm(rx_p);
BS.loc = Rx;                          % BS's position
BS.pol = rx_p;


Ue1 = [150,-20,0;150,-20,1.9];
Ue2 = [100,-10,0;100,-10,1.5];
Ue3 = [50,20,0;50,20,0.5];
pUe = [Ue1; Ue2; Ue3];


UE.loc = zeros(UE_Num, 3);
UE.pol = zeros(UE_Num, 3);
figure(100)
for i = 1:UE_Num
    Ue = pUe(i*2-1:i*2, :);
    plot3(Ue(:,1),Ue(:,2),Ue(:,3),'linewidth',3)
    plot3(Ue(2,1),Ue(2,2),Ue(2,3),'s', 'MarkerFaceColor','r','MarkerSize',3)
    
    
    % Users location 
    UE.loc(i,:) = Ue(2,:);
    LoS_vect = [Rx;Rx - (Rx - Ue(2,:))/norm((Rx - Ue(2,:)))];
    plot3(LoS_vect(:,1),LoS_vect(:,2),LoS_vect(:,3))
    
    LoS_dec = -Oa*(Rx - Ue(2,:))'/norm((Rx - Ue(2,:)));
    
    theta = 180*acos(LoS_dec(1))/pi;
    phi = 180*atan2(LoS_dec(3), LoS_dec(2))/pi;
    disp(['LoS angles: theta = ',num2str(theta), '; phi = ',num2str(phi)])
    
    % Users polarization
    tx_p = [0, 0, 1];                     % Polarization vector of UE
    UE.pol(i,:) = tx_p/norm(tx_p);
end

% LoS angles: theta = 17.077; phi = -72.2415
% LoS angles: theta = 27.7592; phi = -145.8336
% LoS angles: theta = 52.7218; phi = -163.0913