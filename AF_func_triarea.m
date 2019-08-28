function S = AF_func_triarea(p1,p2,p3)

v1 = p1 - p2;
v2 = p3 - p2;

cos_alpha = v1*v2'/(norm(v1)*norm(v2));
sin_alpha = sqrt(1-cos_alpha^2);

h = norm(v1)*sin_alpha;
b = norm(v2);
S = h*b/2;