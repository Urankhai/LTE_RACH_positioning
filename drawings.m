figure

plot(abs(reshape(UE1_sign(1,1,:),1,[])))
hold on
plot(abs(reshape(UE1_sign(1,2,:),1,[])),'r')
plot(abs(reshape(UE1_sign(1,3,:),1,[])),'g')
plot(abs(reshape(UE1_sign(1,4,:),1,[])),'c')

figure

plot(abs(reshape(UE2_rach(1,1,:),1,[])))
hold on
plot(abs(reshape(UE2_rach(1,2,:),1,[])),'r')
plot(abs(reshape(UE2_rach(1,3,:),1,[])),'g')
plot(abs(reshape(UE2_rach(1,4,:),1,[])),'c')


figure

plot(abs(reshape(UE3_rach(1,1,:),1,[])))
hold on
plot(abs(reshape(UE3_rach(1,2,:),1,[])),'r')
plot(abs(reshape(UE3_rach(1,3,:),1,[])),'g')
plot(abs(reshape(UE3_rach(1,4,:),1,[])),'c')