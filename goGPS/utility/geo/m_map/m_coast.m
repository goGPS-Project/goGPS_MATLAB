function h=m_coast(varargin);
% M_COAST Add a coastline to a given map.
%         M_COAST draw a coastline as either filled patches (slow) or
%         lines (fast) on a given projection. It uses a coastline database with
%         a resolution of about 1/4 degree. 
%
%         M_COAST( (standard line option,...,...) )  or
%         M_COAST('line', (standard line option,...,...) ) draws the coastline
%         as a simple line.
%         M_COAST('patch' ( ,standard patch options,...,...) ) draws the 
%         coastline as a number of patches. 
%
%    
%         See also M_PROJ, M_GRID     

% Rich Pawlowicz (rich@ocgy.ubc.ca) 15/June/98
%
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.


% Set current projection to geographic
Currentmap=m_coord('set');
m_coord('geographic');


h=mu_coast('default',varargin{:},'tag','m_coast');

m_coord(Currentmap.name);
