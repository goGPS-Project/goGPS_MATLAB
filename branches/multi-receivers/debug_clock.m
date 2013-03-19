 close all
figure
plot(dtR(:,1,1))


figure
plot(dtR(:,1,1)-(time_R(:,:,1)-time_GPS),'*b')
hold on
plot(((time_R(:,:,1)-time_GPS))*(( max(abs((time_R(:,:,1)-time_GPS)))/max(abs((dtR(:,1,1)-(time_R(:,:,1)-time_GPS))))))^-1,'-r')
