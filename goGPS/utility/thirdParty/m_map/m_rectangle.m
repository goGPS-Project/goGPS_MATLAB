function h=m_rectangle(long,lat,width,height,varargin)
%M_RECTANGLE Create rectangle, rounded-rectangle, or ellipse
%   M_RECTANGLE adds a default rectangle to the current map axes.
%   
%   M_RECTANGLE(long,lat,width,height) creates a rectangle in 2-D coordinates.
%   The long and lat elements determine the location and the width and height
%   elements determine the size. The function plots into the current axes
%   without clearing existing content from the axes.
%   M_RECTANGLE(long,lat,width,height,flag) creates different rectangle. If
%   flag=0(default),it is a normal rectangle. If flag=1, it is a curve box
%   determined to the projection.
%   
%   M_RECTANGLE(...,'Curvature',cur) adds curvature to the sides
%   of the rectangle. For different curvatures along the horizontal and
%   vertical sides, specify cur as a two-element vector of the form
%   [horizontal vertical]. For the same length of curvature along all
%   sides, specify cur as a scalar value. Specify values between 0 (no
%   curvature) and 1 (maximum curvature). Use [1 1] to create an ellipse or
%   circle.
%
%   M_RECTANGLE(...,Name,Value) specifies rectangle properties using one or
%   more Name,Value pair arguments. It is the same as function rectangle.
%   
%   M_RECTANGLE(container,...) creates the rectangle in the axes, group, or
%   transform specified by container, instead of in the current axes.
%   
%   R = M_RECTANGLE(...) returns the rectangle object created.
%
%   Execute GET(R), where R is a rectangle object, to see a list of
%   rectangle object properties and their current values.
%   Execute SET(R) to see a list of rectangle object properties and legal
%   property values.
%


% Shi Weiheng (tfoterye@gmail.com) 25/Mar/18
%
% Based on m_quiver written by Prof Rich Pawlowicz.
% Comments mainly come from function rectangle.
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

if (varargin{1}==0 || varargin{1}==1)
    flag=varargin{1};
    varargin1={varargin{2:end}};
else
    flag=0;
    varargin1=varargin;
end


[X,Y]=m_ll2xy(long,lat,'clip','off');
[X1,Y1]=m_ll2xy(long+width,lat+height,'clip','off');
[X2,Y2]=m_ll2xy(long,lat+height,'clip','off');
[X3,Y3]=m_ll2xy(long+width,lat,'clip','off');

[small_x1,small_y1]=m_ll2xy(long:width/100:long+width,lat*ones(1,101),'clip','off');
[small_x2,small_y2]=m_ll2xy(long:width/100:long+width,(lat+height)*ones(1,101),'clip','off');

if flag==0
h=rectangle('Position',[X,Y,X1-X,Y1-Y],varargin1{:});
end
if flag==1
   h=line( [X,small_x1,X1,small_x2(end:-1:1),X],[Y,small_y1,Y1,small_y2(end:-1:1),Y],varargin1{:});
end
set(h,'tag','m_rectangle');
if nargout==0
 clear h
end
