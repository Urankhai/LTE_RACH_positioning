close all

w = 2*pi*3.2e9; % Frequency of the carrier
c = 299792458;  % Speed of light


BS_ant = [0.05*(1:Ant_Num)',zeros(Ant_Num,1)];
BS_loc = [0.05, 0];

circle = c*corrected_phases/w;

ttt = 0:0.01:1;
figure
hold on
plot(BS_ant(:,1),BS_ant(:,2),'*')
for rrr = 1:Ant_Num
    x_small = circle(rrr)*cos(2*pi*ttt) + BS_ant(rrr,1);
    y_small = circle(rrr)*sin(2*pi*ttt) + BS_ant(rrr,2);
    plot(x_small,y_small)
end

collection_point = [0.2 -0.1];
temp_index = 1:8;
local_coordinates = BS_ant(temp_index,:) - ones(length(temp_index),1)*collection_point;

sigma = 4*10^(-4);
W = sigma*eye(length(temp_index)); %;ones(length(temp_index))*sigma+eye(length(temp_index))-sigma*eye(length(temp_index));

B = [local_coordinates, (circle(temp_index))'];
KB = B'*W^(-1)*B;
disp(['log10 of condition number of KB = ', num2str(log10(cond(KB)))])
invB = (KB)^(-1)*B'*W^(-1);
M = [1 0 0 ; 0 1 0 ; 0 0 -1];
I = ones(length(temp_index), 1);
a = zeros(length(temp_index), 1);

for k = 1:length(temp_index)
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

plot(solution1(1)+collection_point(1),solution1(2)+collection_point(2),'*r')
plot(solution2(1)+collection_point(1),solution2(2)+collection_point(2),'*g')

t=0:0.001*pi:2*pi;
est_x1=solution1(3)*cos(t)+(solution1(1)+collection_point(1));
est_y1=solution1(3)*sin(t)+(solution1(2)+collection_point(2));
plot(est_x1,est_y1);

est_x2=solution2(3)*cos(t)+(solution2(1)+collection_point(1));
est_y2=solution2(3)*sin(t)+(solution2(2)+collection_point(2));
plot(est_x2,est_y2);


if solution1(2) + collection_point(2) > 0
    UE_location = solution1(1:2)' + collection_point;
    UE_radius = -solution1(3);
else
    disp('need to be investigated')
end

figure
hold on
est_x=UE_radius*cos(t)+UE_location(1);
est_y=UE_radius*sin(t)+UE_location(2);
plot(est_x,est_y);
plot(UE_location(1),UE_location(2),'s','MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',10)


plot(BS_ant(:,1),BS_ant(:,2),'*')
for rrr = 1:Ant_Num
    x_small = circle(rrr)*cos(2*pi*ttt) + BS_ant(rrr,1);
    y_small = circle(rrr)*sin(2*pi*ttt) + BS_ant(rrr,2);
    plot(x_small,y_small)
end


disp(['User location is [', num2str(UE_location(1)),', ', num2str(UE_location(2)),'], radius of the circle = ', num2str(UE_radius)])