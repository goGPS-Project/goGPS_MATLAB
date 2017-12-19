clear
close all
t = zeros(1,100);
idx = [5 50:100];
t(idx) = 1;
t = t +randn(1,100)*0.2;
t = movmedian(t,5);
plot(t)
dt = Core_Pre_Processing.diffAndPred(t',1);
dtr = flip(diff(flip(t),1))
figure
plot(abs(dtr(1:end-1)) + abs(dtr(2:end)))


filter = 
