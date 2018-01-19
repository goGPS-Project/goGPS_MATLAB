function [values,longs,lats]=m_elev(varargin)
% M_ELEV Contour elevation onto a map using a 1-degree database
%        M_ELEV contours elevations at 1000m intervals for the map.
%        M_ELEV(OPTN (,LEVELS) (,ARGS,...) ) lets you change various options.
%        if OPTN=='contour', contour lines are drawn. For OPTN=='contourf',
%        filled contours are drawn. LEVELS are the levels used, and ARGS
%        are optional patch arguments of line types, colors, etc. 
%
%        [CS,H]=M_ELEV(...) allows access to the return arguments of the
%        contour/contourf call.
%
%        [ELEV,LONG,LAT]=M_ELEV([LONG_MIN LONG_MAX LAT_MIN LAT_MAX])
%        extracts elevation data for the given lat/long limits (without plotting).
%
%        See also M_PROJ, M_GRID, M_COAST

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% 20/5/97 - moved registration over by 1/2 degree (seems to fit better)
% 2/10/97 - a minor bug in the edge-handling for filled contours!
% 8/ 1/98 - better handling for options.
% 23/1/98 - redid everything to allow for raw bathymetry output option.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 4/DEc/11 - isstr to ischar
% 30/Aug/13 - hack edge-handling in azimuthal projections that wrap
%             around.

global MAP_PROJECTION MAP_VAR_LIST 

% Have to have initialized a map first

draw_map=1;
if nargin==1 && ~ischar(varargin{1}) && length(varargin{1})==4
  draw_map=0;
end


if draw_map==1 && isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

if MAP_PROJECTION.IsOctave
  disp('The elevation database called by ''m_elev.m'' does not exist in Octave');
  return;
end
  
% Set current projection to geographic
Currentmap=m_coord('set');
m_coord('geographic');

if draw_map

  blat=max(floor(MAP_VAR_LIST.lats(1)+.5),-89)-.5;
  tlat=min(ceil(MAP_VAR_LIST.lats(2)+.5),90)-.5;
  llong=floor(MAP_VAR_LIST.longs(1)+.5)-.5;
  rlong=ceil(MAP_VAR_LIST.longs(2)+.5)-.5;

else

  blat=max(floor(varargin{1}(3)+.5),-89)-.5;
  tlat=min(ceil(varargin{1}(4)+.5),90)-.5;
  llong=floor(varargin{1}(1)+.5)-.5;
  rlong=ceil(varargin{1}(2)+.5)-.5;

end

load topo

if rlong>360, rlong=rlong-360; llong=llong-360; end
if llong<-360, rlong=rlong+360; llong=llong+360; end

lts=(blat:tlat);
lgs=(llong:rlong);

% RP Aug/2013
% This is a hack necessary to make filled contours come out correctly
% in cases where we have to handle a 'wrap around'. Generally we
% want to get data past the map edges (so we can clip to map boundaries),
% but if there is a 'wrap' then this creates patches that overlap on
% themselves, which (when filled) show up in the background colour.
if strcmp(MAP_PROJECTION.routine,'mp_azim') && length(lgs)==362
  lgs=lgs(1:361);
  rlong=lgs(end); 
end
  

if rlong<0
  topo=topo(lts+90.5,lgs+360.5);
elseif llong<0 && rlong>=0
  topo=topo(lts+90.5,[(360.5+llong:end) (1:rlong+0.5)]);
else
  topo=topo(lts+90.5,lgs+.5);
end

 
if draw_map

  if nargin==0
   levels=[-7000:1000:-1000 1000:1000:5000];
   optn='contour';
   n_opt=1;
  else
   if ischar(varargin{1})
     optn=varargin{1};
   end
   if nargin==1
     levels=[-7000:1000:-1000 1000:1000:5000];
     n_opt=2;
   else
     if ischar(varargin{2})
       levels=[-7000:1000:-1000 1000:1000:5000];
       n_opt=2;
    else
       levels=varargin{2};
       n_opt=3;
     end
   end
  end

  [lg,lt]=meshgrid(lgs,lts);
  %[X,Y]=m_ll2xy(lg,lt,'clip','on');


  hold on;
  switch optn
   case 'contour'
      [values,longs]=m_contour(lg,lt,topo,levels);
   case 'contourf'
      [values,longs]=m_contourf(lg,lt,topo,levels);
   case 'pcolor'
      [longs]=m_pcolor(lg,lt,topo);
  end  

  set(longs,'tag','m_elev');  
  if n_opt<length(varargin)
      for l=1:length(longs) 
          set(longs(l),varargin{n_opt:end}); 
      end 
  end

else

  [longs,lats]=meshgrid(lgs,lts);
  values=topo;

end

% Reset map coords
m_coord(Currentmap.name);

if nargout==0
 clear values lats longs
end
