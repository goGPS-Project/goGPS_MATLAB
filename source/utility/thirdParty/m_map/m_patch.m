function h=m_patch(long,lat,C,varargin)
% M_PATCH Create patches on a map
%    M_PATCH(LONG,LAT,C) is a drop-in replacement for PATCH that uses 
%    longitude/latitude coordinates to draw a patch on the current map. 
%    See PATCH for more details about the way in which patch colours and 
%    properties should be specified.
%
%    Currently you cannot specify C to be other than a string or 1x3 RGB
%    vector.
%
%    See also M_LINE, M_LL2XY

% Rich Pawlowicz (rich@ocgy.ubc.ca) 3/Sep/98
% 
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

%  10/Mar/99 - changed order of calls ('c' not handled correctly in mu_coast otherwise)
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)

[m,n]=size(long);

if m==1 && n>1
  h=mu_coast('vector',[long' lat';long(1) lat(1)],'patch',C,'tag','m_patch',varargin{:});
elseif m>1 && n==1
  h=mu_coast('vector',[long lat;long(1) lat(1)],'patch',C,'tag','m_patch',varargin{:});
else
  h=mu_coast('vector',[reshape([long;long(1,:);NaN+ones(1,n)],(m+2)*n,1),...
                     reshape([lat;lat(1,:);NaN+ones(1,n)],(m+2)*n,1)],'patch',C,'tag','m_patch',varargin{:});
end

if nargout==0
 clear h
end
