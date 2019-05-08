function [cs,h]=m_contourf(long,lat,data,varargin)
%  M_CONTOURF Adds filled contours to a map
%    M_CONTOURF(LONG,LAT,DATA,...) is the same as M_CONTOUR except
%    that contours are filled. Areas of data above a given level are
%    filled, areas below are left blank or are filled by a lower level.
%    NaN's in the data leave holes in the filled plot/
%
%    [CS,H] = M_CONTOURF(...) returns contour matrix C as described in 
%    CONTOURC and a vector H of handles to PATCH objects (for use by
%    CLABEL).
%
%    See also M_CONTOUR, CONTOURF

% Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 19/02/98 - type - should have been 'clip','patch', rather than 'off'.
%  9/12/98 - handle all-NaN plots without letting contour crash.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% Apr/06  - workaround for v7 bug in contourf.


global MAP_PROJECTION 

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end
 
if min(size(long))==1 && min(size(lat))==1
 [long,lat]=meshgrid(long,lat);
end
 
[X,Y]=m_ll2xy(long,lat,'clip','on');  %First find the points outside
 
i=isnan(X);      % For these we set the *data* to NaN...
data(i)=NaN;
 
                 % And then recompute positions without clipping. THis
                 % is necessary otherwise contouring fails (X/Y with NaN
                 % is a no-no. Note that this only clips properly down
                 % columns of long/lat - not across rows. In general this
                 % means patches may nto line up properly a right/left edges.
if any(i(:)), [X,Y]=m_ll2xy(long,lat,'clip','patch'); end  

if any(~i(:))

 % Bug in contourf call - Solution Number: 1-1W36E8
 % (that involved whether or not the X/Y matrices were handled
 % correctly). In later versions a problem comes up with the renderer
 % that could be solved here, but is better handled in m_grid (see comments
 % in code there).
  if any(strfind(version,'(R14) Service Pack 3')) || ...
     any(strfind(version,'7.2.0.294 (R2006a)')) || ...
     any(strfind(version,'7.2.0.294 (R2006a)')) %%| ....
   %%  any(strfind(version,'7.4.0.336 (R2007a)')),
   [cs,h]=contourf('v6',X,Y,data,varargin{:});
  else
   [cs,h]=contourf(X,Y,data,varargin{:});
  end
   
 set(h,'tag','m_contourf');
else
  cs=[];h=[];
end

if nargout==0
 clear cs h
end
