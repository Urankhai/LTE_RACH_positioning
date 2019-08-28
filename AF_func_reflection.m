function [NLoS_dep, K, NLoS_arr, NLoS_dist, temp_cos, temp_M, belonging]= AF_func_reflection(BS, UE, plane)

% Given plane
n1 = plane(1:3);
D1 = plane(4);
if n1(3) ~= 0
    p1 = 0;
    p2 = 0;
    p3 = -D1/n1(3);
elseif n1(2) ~= 0
    p1 = 0;
    p3 = 0;
    p2 = -D1/n1(2);
else
    p2 = 0;
    p3 = 0;
    p1 = -D1/n1(1);
end
a_point_on_plane = [p1,p2,p3];

% searching the mirror point from another side of the plane
v1 = UE - a_point_on_plane;
v2 = (v1*n1')*n1;
S_UE = v1 - v2 + a_point_on_plane;
mirror_UE = S_UE - v2;

% searching of reflection point K
temp_vect = BS - mirror_UE;
temp_dist = norm(temp_vect);
temp_ort  = temp_vect/temp_dist;
temp_cos  = abs(temp_ort*n1');

n2 = cross(temp_ort, n1);
n2 = n2/norm(n2); 
D2 = -BS*n2';

n3 = cross(n2, temp_ort);
n3 = n3/norm(n3);
D3 = -BS*n3';

temp_M = [n1; n2; n3];
temp_D = [D1; D2; D3];
K = (temp_M\-temp_D)';

belonging = AF_func_belonging(K, plane);

NLoS_dep = (K - UE)/norm(K - UE);
NLoS_arr = temp_ort;
NLoS_dist = norm(K - UE) + norm(BS - K);
