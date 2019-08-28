% calculation of reflection
figure(100+Ant_Num)
est_p = UE_location;
% plot(est_p(1),est_p(2),')

sp1 = [est_p(1),est_p(2) + a1]*[-1;1]/norm([est_p(1),est_p(2) + a1])/sqrt(2);
dd1 = norm([est_p(1),est_p(2) + a1])*sin(acos(sp1));
ref1 = est_p + sqrt(2)*dd1*[-1, -1];


sp2 = [est_p(1),est_p(2) + a2]*[1;1]/norm([est_p(1),est_p(2) + a2])/sqrt(2);
dd2 = norm([est_p(1),est_p(2) + a2])*sin(acos(sp2));
ref2 = est_p + sqrt(2)*dd2*[1, -1];


plot(ref1(1),ref1(2), 'd')
plot(ref2(1),ref2(2), 'd')

keyboard
