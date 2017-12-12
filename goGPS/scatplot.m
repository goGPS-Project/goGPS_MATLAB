% Draws a modified "polar" plot. The N data sets are assigned to a 360/N degree sector, 
% and each point within a data set is plotted inside the appropriate sector. This is done 
% by computing the point's polar coordinate, and mapping the resulting angle component 
% on to the appropriate 360/N degree sector, while retaining the point's magnitude. 
%
% Syntax: scatplot(x_1,y_1,S_1,x_2,y_2,S_2,...)
% x_n and y_n are the nth data sets, and S_n is the corresponding symbol string.
%
% Optional parameter syntax: scatplot(x_1,y_1,S_1,x_2,y_2,S_2,...,S)
% S is one of the following strings:
%       'none' does not draw any marker
%       'c'    only draws the circle marker
%       'l'    only draws the line marker
%       'cl'   draws the circle and line marker (default)
%       'lc'   is the same as 'cl'
%
% Example usage: 
% a=randn(2,1000);
% b=randn(2,1000)*2;
% c=randn(2,1000)*3;
% scatplot(a(1,:),a(2,:),'rd',b(1,:),b(2,:),'bs',c(1,:),c(2,:),'go');
%
%
%
% Hasan Mir
% Email: hmir@mit.edu
% 12/2/2005


%------------------------------------------------------------------%
function scatplot(varargin)

if mod(nargin,3)~=0
    
    N=(nargin-1)/3;
    
else 
    
    N=nargin/3;
    
end

sect_width=2*pi/N;    
offset_angle=0:sect_width:2*pi-sect_width;

R=0;
m=1;

for n=1:N
    
    x=varargin{m};
    y=varargin{m+1};
    
    phi1=offset_angle(n);
    phi2=sect_width;
    theta1=atan2(y,x)+pi;
    theta2=phi1+phi2*theta1/(2*pi);
    
    r=sqrt(x.^2+y.^2);
    z=r.*exp(j*theta2);
    x=real(z); 
    y=imag(z);
    
    if max(r)>R
        
        R=max(r);
        
    end
    
    plot(x,y,varargin{m+2})
    hold on;
    
    m=m+3;
    
end

if mod(nargin,3)==0
    
    drawline(R,offset_angle)
    drawcircle(R)
    
else
    
    switch varargin{end}
        
        case 'none'
            
        case 'l'
            
            drawline(R,offset_angle)
            
        case 'c'
            
            drawcircle(R)
            
        case {'lc','cl'}
            
            drawline(R,offset_angle)
            drawcircle(R)
            
    end
    
end

hold off; axis equal
%------------------------------------------------------------------%


%------------------------------------------------------------------%
function drawline(R,offset_angle)

for n=1:length(offset_angle)
    
    plot(real([0 R]*exp(j*offset_angle(n))),imag([0 R]*exp(j*offset_angle(n))),'k-')
    
end
%------------------------------------------------------------------%


%------------------------------------------------------------------%
function drawcircle(R)

r=linspace(0,R,5);
w=0:.01:2*pi;

for n=2:length(r)
    
    plot(real(r(n)*exp(j*w)),imag(r(n)*exp(j*w)),'k--')
    
end
%------------------------------------------------------------------%