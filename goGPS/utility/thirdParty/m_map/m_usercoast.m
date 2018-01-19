function m_usercoast(varargin)
% M_USERCOAST Add a coastline using a user-specified coastline file.
%         M_USERCOAST(FILENAME) uses data previously extracted and stored
%         in FILENAME to draw a coast.
%
%         M_USERCOAST(... ,(standard line option,...,...) ) draws the coastline
%         as a simple line.
%         M_USERCOAST(..., 'patch' ( ,standard patch options,...,...) ) draws the 
%         coastline as a number of patches. 
%         M_USERCOAST(..., 'speckle' ( ,standard line options,...,...) ) draws the 
%         coastline as a speckled line. 
%
%    
%         See also M_PROJ, M_GRID, M_COAST, M_GSHHS_C    

% Rich Pawlowicz (rich@ocgy.ubc.ca) 15/June/98
%
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

mu_coast('user',varargin{:},'tag','m_usercoast');

