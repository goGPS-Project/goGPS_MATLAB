function h=m_gshhs(resolution,varargin)
% M_GSHHS Add a coastline to a given map using 
%           the Global Self-consistant Hierarchical High-resolution 
%           Shorelines, Rivers, and Borders
%
%         M_GSHHS(RES, (standard line option,...,...) ) draws the coastline
%         river network, or borders as  simple lines.
%
%         M_GSHHS(RES,'patch' ( ,standard patch options,...,...) ) draws the 
%         coastline as a number of patches (rivers and borders are not
%         arranged so patches can be drawn).
%
%         M_GSHHS(RES,'save',FILENAME) saves the extracted coastline data
%         for the current projection in a file FILENAME. This allows 
%         speedier replotting using M_USERCOAST(FILENAME). 
%    
%         RES: A one-char string (optionally 2 or 3)
%         
%         First char: resolution - one of
%                      'c'  crude
%                      'l'  low
%                      'i'  intermediate
%                      'h'  high
%                      'f'  full
%
%         Second char: type - one of
%                      'c' GSHHS coastline (default)
%                      'b' WDB Border
%                      'r' WDB River
%  
%         Third char - if 2nd char is 'b':
%                      '1' Country borders
%                      '2' State/Province and Country borders
%                    - if 2nd char is 'r': '1','2','3','4' 
%                      add successively more tributaries
%
%         (also maintained is this optional format:
%
%         RES - selections resolution
%                  1  or 'crude'	
%                  2  or 'low'  	
%                  3  or 'intermediate'  
%                  4  or 'high' 	
%                  5  or 'full  	
%
%         but please don't use this).
%
%         See also M_PROJ, M_GRID, M_COAST, M_GSHHS_L, M_GSHHS_H, M_GSHHS_C 
%         M_USERCOAST    

% Rich Pawlowicz (rich@ocgy.ubc.ca) 15/June/98
%
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
%  16/Dec/2005
%*********************************************************************
%  Modified after code provided by Bruce Lipphardt (brucel@udel.edu) to 
%  reduce the hierarchy of M_GSHHS_* routines to a single routine with a
%  variable resolution input:
% 20/Jan/2008 - added borders and rivers from gshhs v1.10
% 4/DEc/11 - isstr to ischar
% Sep/14 - added hierarchy to borders



% Root of directories where all the gshhs_X.b, wdb_borders-X.b and wdb_rivers_X.b
% files live
FILNAME='private/';


%-------------don't change below here----------------------------

res_list = char('c','l','i','h','f') ;
typ_list=char('c','b','r');
typ_names={'gshhs_','wdb_borders_','wdb_rivers_'};

typ=1;
flaglim='9';

if ischar(resolution)
 if length(resolution)>=2
   typ = strmatch(lower(resolution(2)),typ_list);
 end
 if length(resolution)>=3
   flaglim = resolution(3);
 end  
 resolution = strmatch(lower(resolution(1)),res_list);
end
 
 
if isempty(resolution) || resolution<1 || resolution> length(res_list)
  error('**Don''t recognize the specified resolution');
end
if isempty(typ) || typ<1 || typ> length(res_list)
  error('**Don''t recognize the specified type');
end
  
res_char = res_list(resolution) ;
file     = [FILNAME,sprintf('%s%s.b',typ_names{typ},res_char)] ;
tag_name = sprintf('m_%s%s',typ_names{typ},res_char) ;


% Set current projection to geographic
Currentmap=m_coord('set');
m_coord('geographic');


if length(varargin)>1 && strcmp(varargin{1},'save')
  [ncst,Area,k]=mu_coast(res_char,file);
  save(varargin{2},'ncst','k','Area');
else
  h=mu_coast([res_char flaglim],file,varargin{:},'tag',tag_name);
end

m_coord(Currentmap.name);

if nargout==0
    clear h
end
