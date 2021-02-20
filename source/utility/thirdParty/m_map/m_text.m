function h=m_text(long,lat,varargin)
% M_TEXT Text Annotation
%    M_TEXT(LONG,LAT,'string') adds the text in the quotes to location 
%    (LONG,LAT) on the currently defined map projection. If LONG and LAT
%    are vectors, M_TEXT writes the text at all locations given. 
%    If 'string' is an array the same number of rows as the
%    length of LONG and LAT, M_TEXT marks each point with the 
%    corresponding row of the 'string' array.
%     
%    M_TEXT returns a column vector of handles to TEXT objects, one
%    handle per text object. TEXT objects are children of AXES objects.
% 
%    M_TEXT(LONG,LAT,'string',property/value pairs) can be used to
%    change fontsize, weight, color, etc using the standard TEXT
%    properties.
%
%    See also TEXT.

% Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


% 31/Jul/99 - changed to allow for X/Y vectors.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)

global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

[X,Y]=m_ll2xy(long,lat,'clip','off');
%h=text('position',[X(:) Y(:)],'tag','m_text','string',varargin{:});
% Fix to allow vectors of X/Y to work
h=text(X(:),Y(:),varargin{:});
if ~isempty(h) && isempty(get(h(1),'tag'))
 set(h,'tag','m_text');
end

if nargout==0
 clear h
end
