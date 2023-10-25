clc;
clear all;
%% 已知条件
z=[-5 -2 0]';
p=[-1 -3 -4+i*3 -4-i*3]';
k=5;
t=0:0.001:8;
%% 传递函数
% G=zpk(z,p,k)
[num,den]=zp2tf(z,p,k);
G=tf(num,den)
%% 求单位阶跃响应
% step(num,den)
% step(G)
figure(1);
subplot(121);
y=step(t,G);
plot(t, y,'g-','linewidth',1);
title('单位阶跃响应');
xlabel('T/s');ylabel('幅度');
grid on;box on;

%% 求单位冲击响应
% y = impulse(num,den);
y = impulse(G,0:0.001:8);
subplot(122);
plot(t, y,'r-','linewidth',1);
grid on;box on;
title('单位冲击响应');
xlabel('T/s');ylabel('幅度');