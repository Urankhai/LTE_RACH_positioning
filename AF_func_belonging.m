function answer = AF_func_belonging(K, plane)
epsilon = 1e-10;
p = zeros(4,3);
for k = 1:4
    p(k,:) = plane(5+(k-1)*3: 4+k*3);
end

S_quad = AF_func_triarea(p(1,:),p(2,:),p(3,:)) + AF_func_triarea(p(1,:),p(4,:),p(3,:));
S = AF_func_triarea(p(1,:),p(2,:),K) + AF_func_triarea(p(2,:),p(3,:),K) + AF_func_triarea(p(3,:),p(4,:),K)+ AF_func_triarea(p(4,:),p(1,:),K);

% keyboard

if abs(S-S_quad) < epsilon
    answer = 1;
else
    answer = 0;
% %     keyboard
%     figure
%     hold on
%     plot3(K(1),K(2),K(3),'ro')
%     patch(p(:,1),p(:,2),p(:,3),[.8 .8 .8])
end
