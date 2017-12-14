function polarScatter(az,decl,size,color,flag)
%%% INTERNAL PARAMETER
scale = 1;
%%%

decl_n = decl/(pi/2)*scale;
az = az + pi/2;
x = cos(az).*decl_n;
y = sin(az).*decl_n;
if nargin > 4 && strcmp(flag,'filled')
    scatter(x,y,size,color,'filled')
else
    scatter(x,y,size,color)
end
hold on
%plot parallel
az_l = [0:pi/200:2*pi];
d_step = 15/180*pi;
decl_s = ([0:d_step:pi/2]/(pi/2))*scale;
for d = decl_s
    x = cos(az_l).*d;
    y = sin(az_l).*d;
    text(cos(80/180*pi)*d,sin(80/180*pi)*d,sprintf('%d',round(d*90)),'HorizontalAlignment','center');
    plot(x,y,'color',[0.6 0.6 0.6]);
    
end
%plot meridian
az_step = 30/180 *pi;
az_s = [0:az_step:2*pi];
decl_l = ([0 1])*scale;
for a = az_s
    x = cos(a).*decl_l;
    y = sin(a).*decl_l;
    if abs(a-2*pi) > 0.0001
        text(cos(a)*1.1,sin(a)*1.1,sprintf('%d',round(a/pi*180)),'HorizontalAlignment','center');
    end
    plot(x,y,'color',[0.6 0.6 0.6]);
    
end
axis equal
% xlim([-2 2])
% ylim([-2 2])
axis off
set(gcf,'color','w');
hold off
end