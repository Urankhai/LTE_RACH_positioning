function [G_perp, G_parr] = AF_func_Fresnel(cos_theta, er)

root_er = sqrt(er - 1 + cos_theta^2);
G_perp = (cos_theta - root_er)/(cos_theta + root_er);
G_parr = (er*cos_theta - root_er)/(er*cos_theta + root_er);
