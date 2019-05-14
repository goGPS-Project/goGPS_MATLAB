function g=plotEarth(b)
g = hgtransform;
if nargin==0
b=6378137.0; %% semimajor axis ellipsoid [m]
end
coast=loadjson('coast.geojson');
for i=1:length(coast.features);
    toplot=LL2V(coast.features(i).geometry.coordinates/180*pi)*b;
    plot3(toplot(:,1),toplot(:,2),toplot(:,3),'Color',[0.4 0.4 0.4],'Parent',g)
    hold on

end
color_sec=[0.6 0.6 0.6];
second_lw=0.1;
r=b;
p=10; % parallel each 20 degree
phi=linspace(0,pi,p);
phi(1)=[];
phi(end)=[];
m=19; % meridian each 20 degree
lambda=linspace(0,2*pi,m);
lambda(1)=[];
lambda(end)=[];
teta1=-pi/2+pi/(p-1):0.01:pi/2-pi/(p-1);
teta2=[pi/2+pi/(p-1):0.01:pi,-pi:0.01:-pi/2-pi/(p-1)];
% teta(and(teta<pi/2+pi/(p-1),teta>pi/2-pi/(p-1)))=[];
% teta(and(teta<-pi/2+pi/(p-1),teta>-pi/2-pi/(p-1)))=[];
x1=r*cos(teta1);
y1=r*sin(teta1);
x2=r*cos(teta2);
y2=r*sin(teta2);
for i=1:m-2
    plot3(cos(lambda(i))*x1,sin(lambda(i))*x1,y1,'Color',color_sec,'LineWidth',second_lw,'LineStyle','-','Parent',g)
    hold on
    plot3(cos(lambda(i))*x2,sin(lambda(i))*x2,y2,'Color',color_sec,'LineWidth',second_lw,'LineStyle','-','Parent',g)
    hold on
end
r=b;
teta=-pi:0.01:pi;
x=r*cos(teta);
y=r*sin(teta);
for i=1:p-2
    plot3(x*sin(phi(i)),y*sin(phi(i)),r*ones(1,numel(x))*cos(phi(i)),'Color',color_sec,'LineWidth',second_lw,'LineStyle','-','Parent',g) 
 hold on
end
axis equal
axis off
set(gcf,'color','w');
[x,y,z] = sphere;
surf(x*(b-1),y*(b-1),z*(b-1),'FaceColor', [1 1 1],'EdgeColor','none','Parent',g)
end
function V=LL2V(lonlat)
V=[cos(lonlat(:,2)).*cos(lonlat(:,1)) cos(lonlat(:,2)).*sin(lonlat(:,1)) sin(lonlat(:,2))];
end
