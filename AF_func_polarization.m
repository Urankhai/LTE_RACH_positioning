function [att, pol] = AF_func_polarization(p, vect)
p = p/norm(p);
vect = vect/norm(vect);

n1 = cross(vect, p);
n1 = n1/norm(n1);

n2 = cross(n1, vect);

att = p*n2';
pol = att*n2;