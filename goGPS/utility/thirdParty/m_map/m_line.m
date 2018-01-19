function h=m_line(long,lat,varargin)
% M_LINE Create a line on a map
%    M_LINE(LONG,LAT) adds the line in vectors LONG and LAT to the 
%    current axes. If LONG and LAT are matrices the same size, one 
%    line per column is added.
% 
%    LINE returns a column vector of handles to LINE objects, 
%    one handle per line. LINEs are children of AXES objects.
% 
%    The LONG,LAT pair can be followed by 
%    parameter/value pairs to specify additional properties of the lines.
%    These are standard 'line' properties.
%
%    See also LINE, M_LL2XY


% Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998
% 8/Dec/98 - changed switch test to only 3 letters (thus letting
%            "tag" properties through).
% 18/july/00 - Fixed m_line so you could do clipping through it.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)

%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

clp='on';

k=1;
while k<length(varargin)
  switch lower(varargin{k}(1:3))
    case 'cli'
      clp=varargin{k+1};
      if isempty(findstr(clp,'on'))
        varargin{k+1}='off';
      else
        varargin{k+1}='on';
%        varargin([k k+1])=[];
      end
      k=k+2;
    otherwise
      k=k+2;
  end
end

[X,Y]=m_ll2xy(long,lat,'clipping',clp);

if nargout>0
  h=line(X,Y,'tag','m_line',varargin{:});
else
  line(X,Y,'tag','m_line',varargin{:});
end


