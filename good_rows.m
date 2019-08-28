good = [1, 6, 8];

B = [local_coordinates(good,:), (circle(good))'];
KB = B'*B;
disp(['log10 of condition number of KB = ', num2str(log10(cond(KB)))])
invB = (KB)^(-1)*B';
M = [1 0 0 ; 0 1 0 ; 0 0 -1];
I = ones(length(good), 1);
a = zeros(length(good), 1);

for k = 1:length(good)
    a(k) = 0.5* B(k,:)*M*B(k,:)';
end

% c2 * Lambda^2 + c1 * Lambda + c0 = 0
c0 = (invB*a)' * M * invB*a;
c1 = 2*((invB*I)' * M * invB*a - 1);
c2 = (invB*I)' * M * invB*I;

if c1^2 - 4*c2*c0 < 0
    disp('need to check')
    keyboard
end

Lambda1 = (-c1 + sqrt(c1^2 - 4*c2*c0))/(2*c2);
solution1 = M * invB * (a + Lambda1*I);


Lambda2 = (-c1 - sqrt(c1^2 - 4*c2*c0))/(2*c2);
solution2 = M * invB * (a + Lambda2*I);

if solution1(2) + collection_point(2) > BS.loc(2)
    UE_location = solution2(1:2)' + collection_point;
    UE_radius = -solution2(3);
else
    UE_location = solution1(1:2)' + collection_point;
    UE_radius = -solution1(3);
end

% if UE_radius <= 0
%     disp('need to check')
%     keyboard
% end


disp(['     Real UE position = [', num2str(UE.loc(UEi, 1), '%10.2f'),', ', num2str(UE.loc(UEi, 2), '%10.2f'),'];'])
disp(['Estimated UE location = [', num2str(UE_location(1), '%10.2f'),', ', num2str(UE_location(2), '%10.2f'), ']; Rasius of the circle = ', num2str(UE_radius)])

figure(100)
plot(UE_location(1),UE_location(2),'s', 'MarkerEdgeColor','k', 'MarkerFaceColor','c','MarkerSize',6)



t=0:0.001*pi:2*pi; 
est_x=UE_radius*cos(t)+UE_location(1); 
est_y=UE_radius*sin(t)+UE_location(2); 
plot(est_x,est_y); 