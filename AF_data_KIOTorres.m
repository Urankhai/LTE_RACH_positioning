function [planes, BS, Oa, sigma] = AF_data_KIOTorres()
% Parameters of the buildings from Google Maps
bs1 = [-3.689503 40.467352];
bs2 = [-3.688501 40.467197];

Tw1 = [-3.690289 40.467219;
    -3.689897 40.467156;
    -3.689981 40.466841;
    -3.690393 40.466911;
    -3.690289 40.467219];

Tw2 = [-3.688237 40.466895;
    -3.687831 40.466826;
    -3.687954 40.466517;
    -3.688319 40.466581;
    -3.688237 40.466895];

h = 115; % The heights of the buildings in m
a = 15;  % deg, inclination angle
tr = h*tan(a/180*pi);
dtr = h*tan(1/180*pi);
% We use assumption that the ground is horizontal flat
flat_ground = [-20 -40 0; -20 120 0; 240 120 0; 240 -40 0];
% figure
% plot(Tw1(:,1), Tw1(:,2),'-o')
% hold on
% plot(Tw2(:,1), Tw2(:,2),'-o')
% plot(bs1(:,1), bs1(:,2),'*')
% plot(bs2(:,1), bs2(:,2),'*')
% axis equal

L0 = -3; % deg, this information is needed for GaussKrueger
[xbs1, ybs1] = geo2gausskrueger( bs1(1), bs1(2), L0 );
[xbs2, ybs2] = geo2gausskrueger( bs2(1), bs2(2), L0 );
xybs1 = [xbs1, ybs1];
xybs2 = [xbs2, ybs2];

xyTw1 = zeros(size(Tw1));
xyTw2 = zeros(size(Tw2));
for i = 1:length(Tw1)
    [xTw1, yTw1] = geo2gausskrueger( Tw1(i,1), Tw1(i,2), L0 );
    [xTw2, yTw2] = geo2gausskrueger( Tw2(i,1), Tw2(i,2), L0 );
    
    xyTw1(i,:) = [xTw1, yTw1];
    xyTw2(i,:) = [xTw2, yTw2];
end


% figure
% plot(xyTw1(:,1), xyTw1(:,2),'-o')
% hold on
% plot(xyTw2(:,1), xyTw2(:,2),'-o')
% plot(xybs1(:,1), xybs1(:,2),'*')
% plot(xybs2(:,1), xybs2(:,2),'*')
% axis equal

whole = [xybs1;xybs2;xyTw1;xyTw2];
xmin = min(whole(:,1));
ymin = min(whole(:,2));
normalized_coordinates = whole - ones(length(whole),1)*[xmin, ymin];

% figure
% plot(xyTw1(:,1)-xmin, xyTw1(:,2)-ymin,'-o')
% hold on
% plot(xyTw2(:,1)-xmin, xyTw2(:,2)-ymin,'-o')
% plot(xybs1(:,1)-xmin, xybs1(:,2)-ymin,'*')
% plot(xybs2(:,1)-xmin, xybs2(:,2)-ymin,'*')
% axis equal

% Tower 1, ground floor
b1 = normalized_coordinates(3:6,:);
w1 = b1(2,:)-b1(3,:);  % wall, which is opposit to the Tower 2
tr1 = [w1(2), -w1(1)]/norm(w1);
% w1u = b1(2,:)-b1(1,:); % wall up
% w1d = b1(3,:)-b1(4,:); % wall down

% Tower 2, ground floor
b2 = normalized_coordinates(8:11,:);
w2 = b2(1,:)-b2(4,:);  % wall, which is opposit to the Tower 1
tr2 = [-w2(2), w2(1)]/norm(w2);
% w2u = b2(1,:)-b2(2,:); % wall up
% w2d = b2(4,:)-b2(3,:); % wall down

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of standard deviation of the surface height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume the length of w1 = w2 ~= 36
% The total number of windows cells = 30, thus the width of each window =
% 36/30 = 1.2 m.
% Silver surface take 6 packs and the height is 1 m, windows connectors take in total 2 packs and height is 0.02m,
% red line takes 1/5 of the window's length and the height is 0.1m.
expct = 1/24*(0.02*1 + 1/5*0.05);
sigma = sqrt(1/24*(1*(0.02 - expct)^2 + 1/5*(0.05-expct)));
% keyboard


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tower 1, top piece
tp1 = b1 + ones(length(b1),1)*tr*tr1;
vp1 = b1;
dp1 = b1 + ones(length(b1),1)*dtr*tr1;
% Tower 2, top piece
tp2 = b2 + ones(length(b1),1)*tr*tr2;
vp2 = b2;
dp2 = b2 + ones(length(b1),1)*dtr*tr2;

% figure(2)
% plot(tp1(:,1),tp1(:,2),'r-d')
% plot(tp2(:,1),tp2(:,2),'r-d')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Model, Tower 1
base1 = [b1,zeros(length(b1),1)];
top1 = [tp1,h*ones(length(b1),1)];
Tw1 = [base1; top1];
Fs1 = [1 2 6 5; 2 6 7 3; 8 7 3 4; 5 8 4 1; 1 2 3 4; 5 6 7 8];
% Inclined plane calculation
v1 = top1(3,:) - base1(3,:);
n1 = cross([w1,0],v1)/norm(cross([w1,0],v1));   % normal vector
D1 = -base1(3,:)*n1';                           % the element D in equation Ax+By+Cz+D=0
S1 = [n1,D1, base1(2,:),base1(3,:), top1(3,:),top1(2,:)];

vector_n1 = [top1(3,:);top1(3,:)+n1];

% vertical plane calculation
vbase1 = [b1,zeros(length(b1),1)];
vtop1 = [vp1,h*ones(length(b1),1)];
vTw1 = [vbase1; vtop1];
vn1 = [tr1,0];
vD1 = -base1(3,:)*vn1';                          % the element D in equation Ax+By+Cz+D=0
vS1 = [vn1,vD1, vbase1(2,:),vbase1(3,:), vtop1(3,:),vtop1(2,:)];
vector_vn1 = [top1(3,:);top1(3,:)+vn1];

% displaced plane calculation
dbase1 = [b1,zeros(length(b1),1)];
dtop1 = [dp1,h*ones(length(b1),1)];
dTw1 = [dbase1; dtop1];

dv1 = dtop1(3,:) - dbase1(3,:);
dn1 = cross([w1,0],dv1)/norm(cross([w1,0],dv1));
dD1 = -base1(3,:)*dn1';                          % the element D in equation Ax+By+Cz+D=0
dS1 = [dn1,dD1, dbase1(2,:),dbase1(3,:), dtop1(3,:),dtop1(2,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Model, Tower 2
base2 = [b2,zeros(length(b1),1)];
top2 = [tp2,h*ones(length(b1),1)];
Tw2 = [base2; top2];
Fs2 = Fs1;
% Inclined plane calculation
v2 = top2(4,:) - base2(4,:);
n2 = cross(v2,[w2,0])/norm(cross(v2,[w2,0]));   % normal vector
D2 = -base2(4,:)*n2';                           % the element D in equation Ax+By+Cz+D=0
S2 = [n2,D2, base2(1,:),base2(4,:), top2(4,:),top2(1,:)];

vector_n2 = [top2(4,:);top2(4,:)+n2];

% vertical plane calculation
vbase2 = [b2,zeros(length(b1),1)];
vtop2 = [vp2,h*ones(length(b1),1)];
vTw2 = [vbase2; vtop2];

vn2 = [tr2,0];
vD2 = -base2(4,:)*vn2';                          % the element D in equation Ax+By+Cz+D=0
vS2 = [vn2,vD2, vbase2(1,:),vbase2(4,:), vtop2(4,:),vtop2(1,:)];
vector_vn2 = [top2(4,:);top2(4,:)+vn2];

% displaced plane calculation
dbase2 = [b2,zeros(length(b1),1)];
dtop2 = [dp2,h*ones(length(b1),1)];
dTw2 = [dbase2; dtop2];

dv2 = dtop2(4,:) - dbase2(4,:);
dn2 = cross(dv2,[w2,0])/norm(cross(dv2,[w2,0]));   % normal vector
dD2 = -dbase2(4,:)*dn2';                           % the element D in equation Ax+By+Cz+D=0

% dn2 = [tr2,-dh]/norm([tr2,-dh]);
% dD2 = -base2(4,:)*dn2';                          % the element D in equation Ax+By+Cz+D=0
dS2 = [dn2,dD2, dbase2(1,:),dbase2(4,:), dtop2(4,:),dtop2(1,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinates of base stations
ah = 30; % the height of the BS's antenna
% base_station1 = [normalized_coordinates(1,:),0; normalized_coordinates(1,:),ah];
base_station2 = [normalized_coordinates(2,:),0; normalized_coordinates(2,:),ah];
% orthonotmal trihedron
az = [0,0,1];                   % vertical vector coincides with ez
ax = [-w2, 0]/norm([-w2, 0]); % the normal vector of the BS's antenna
ay = cross(az, ax);
Oa = [ax;ay;az];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% drawings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(101)
% hold on
% grid on
% axis equal
% colormap winter
% view(3);
% patch('Vertices',vTw1,'Faces',Fs1,'FaceVertexCData',[.3 .7 .7],'FaceColor','flat')
% patch('Vertices',vTw2,'Faces',Fs2,'FaceVertexCData',[.4 .8 .8],'FaceColor','flat')
% patch(flat_ground(:,1), flat_ground(:,2), flat_ground(:,3),[.8 .8 .8])
% 
% figure(102)
% hold on
% grid on
% axis equal
% colormap winter
% view(3);
% patch('Vertices',dTw1,'Faces',Fs1,'FaceVertexCData',[.3 .7 .7],'FaceColor','flat')
% patch('Vertices',dTw2,'Faces',Fs2,'FaceVertexCData',[.4 .8 .8],'FaceColor','flat')
% patch(flat_ground(:,1), flat_ground(:,2), flat_ground(:,3),[.8 .8 .8])

figure(100)
hold on
grid on
axis equal
plot3(base1(:,1),base1(:,2),base1(:,3))
plot3(top1(:,1),top1(:,2),top1(:,3))
colormap winter
view(3);
patch('Vertices',Tw1,'Faces',Fs1,'FaceVertexCData',[.3 .7 .7],'FaceColor','flat')


plot3(base2(:,1),base2(:,2),base2(:,3))
plot3(top2(:,1),top2(:,2),top2(:,3))
patch('Vertices',Tw2,'Faces',Fs2,'FaceVertexCData',[.4 .8 .8],'FaceColor','flat')
patch(flat_ground(:,1), flat_ground(:,2), flat_ground(:,3),[.8 .8 .8])

% keyboard

plot3(base_station2(:,1),base_station2(:,2),base_station2(:,3),'linewidth',3)
plot3(base_station2(2,1),base_station2(2,2),base_station2(2,3),'.r')


% inclined plane
plot3(vector_n1(:,1),vector_n1(:,2), vector_n1(:,3),'r')
plot3(vector_n2(:,1),vector_n2(:,2), vector_n2(:,3),'r')
% vertical plane
plot3(vector_vn1(:,1),vector_vn1(:,2), vector_vn1(:,3),'g')
plot3(vector_vn2(:,1),vector_vn2(:,2), vector_vn2(:,3),'g')

S_ground = [0, 0, 1, 0, flat_ground(1,:), flat_ground(2,:), flat_ground(3,:), flat_ground(4,:)];
planes = [S_ground; S1;S2;vS1;vS2;dS1;dS2];
BS = base_station2(2,:);
% ax = (UE - BS)/norm(UE - BS); % the normal vector of the BS's antenna
% ay = cross(az, ax);
% Oa = [ax;ay;az];

% figure(101)
% hold on
% grid on
% axis equal
% vect_ax = [BS;BS+5*ax];
% vect_ay = [BS;BS+5*ay];
% vect_az = [BS;BS+5*az];
% 
% plot3(vect_ax(:,1),vect_ax(:,2),vect_ax(:,3),'c','linewidth',2)
% plot3(vect_ax(2,1),vect_ax(2,2),vect_ax(2,3),'d', 'MarkerFaceColor','c','MarkerSize',3)
% % arrow(BS,BS+5*ax,1,'EdgeColor','c','FaceColor','c','linewidth',2)
% text(BS(1)+5*ax(1),BS(2)+5*ax(2),BS(3)+5*ax(3)+0.3,'aO_x');
% 
% plot3(vect_ay(:,1),vect_ay(:,2),vect_ay(:,3),'m','linewidth',2)
% plot3(vect_ay(2,1),vect_ay(2,2),vect_ay(2,3),'d', 'MarkerFaceColor','m','MarkerSize',3)
% % arrow(BS,BS+5*ay,1,'EdgeColor','m','FaceColor','m','linewidth',2)
% text(BS(1)+5*ay(1),BS(2)+5*ay(2),BS(3)+5*ay(3)+0.3,'aO_y');
% 
% plot3(vect_az(:,1),vect_az(:,2),vect_az(:,3),'b','linewidth',2)
% plot3(vect_az(2,1),vect_az(2,2),vect_az(2,3),'d', 'MarkerFaceColor','b','MarkerSize',3)
% % arrow(BS,BS+5*az,1,'EdgeColor','b','FaceColor','b','linewidth',2)
% text(BS(1)+5*az(1),BS(2)+5*az(2),BS(3)+5*az(3)+0.3,'aO_z');
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End drawings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
