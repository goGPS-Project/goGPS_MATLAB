function h=m_gshhs_l(varargin)
% M_GSHHS_L Add a coastline to a given map using the 'low' resolution of
%           the Global Self-consistant Hierarchical High-resolution 
%           Shorelines.
%
%         M_GSHHS_L( (standard line option,...,...) ) draws the coastline
%         as a simple line.
%         M_GSHHS_L('patch' ( ,standard patch options,...,...) ) draws the 
%         coastline as a number of patches. 
%
%         M_GSHHS_L('save',FILENAME) saves the extracted coastline data
%         for the current projection in a file FILENAME. This allows 
%         speedier replotting using M_USERCOAST(FILENAME). 
%    
%         See also M_PROJ, M_GRID, M_COAST, M_GSHHS_H, M_GSHHS_I, M_GSHHS_C 
%         M_USERCOAST    

% Rich Pawlowicz (rich@ocgy.ubc.ca) 15/June/98
%
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% Changes
%  Nov/2017 - changed this into a stub calling gshhs.m (kept
%             for backwards compatability)

h=m_gshhs('lc',varargin{:});

