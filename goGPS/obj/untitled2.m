t=[1000:1:50000];
[X,V]=core_sky.polyInterpolate(t',2);
Vnum=diff(X,1,1);
t=[1000.5:1:49999.5];
[X,V]=core_sky.polyInterpolate(t',2);
scala=1
figure
subplot(3,1,1)
plot(Vnum(:,:,1)-V(:,:,1))
hold on
%plot(V(:,:,1)/scala)
legend('num','an')
subplot(3,1,2)
plot(Vnum(:,:,2)-V(:,:,2))
hold on
%plot(V(:,:,2)/scala)
subplot(3,1,3)
plot(Vnum(:,:,3)-V(:,:,3))
hold on
%plot(V(:,:,3)/scala)