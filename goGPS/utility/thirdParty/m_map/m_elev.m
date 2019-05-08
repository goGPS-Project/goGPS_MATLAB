function [values,longs,lats]=m_elev(varargin)
% M_ELEV Contour elevation onto a map using a 1-degree database
%        M_ELEV contours elevations at 1000m intervals for the map.
%        M_ELEV(OPTN) where OPTN is one of the following strings lets
%        you draw a particular background:
%        OPTN: 'contour' -  contour lines are drawn.
%              'contourf' -  filled contours are drawn. 
%                               LEVELS are the levels used, and ARGS
%                               are optional patch arguments of line types, 
%                               colors, etc. 
%              'pcolor'    - pcolor call
%              'image'     - displays pixellated image
%              'shadedrelief' - shaded relief map.
%        M_ELEV(OPTN,args,...) lets you pass arguments to the OPTN call.
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
% Mar/20/2019 - changed call sequence, other stuff for shaded relief.

global MAP_PROJECTION MAP_VAR_LIST 


% First case if we just want to extract dat, otherwise we are drawing.
if nargin==1 && ~ischar(varargin{1}) && length(varargin{1})==4
  draw_map=0;
else
  draw_map=1;
  optn='contour';
  levels=[-7000:1000:-1000 1000:1000:5000];
  
  if nargin>=1 && ischar(varargin{1})
     optn=varargin{1};
     varargin(1)=[];
  end
end


   

if draw_map

  % Have to have initialized a map first
  if isempty(MAP_PROJECTION)
      disp('No Map Projection initialized - call M_PROJ first!');
      return;
  end
    
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

% Extract topographhy for the desired lat/long limits


% Replace the matlab topography with one I have derived.
currloc=mfilename('fullpath');  % Octave has problems finding the file 
load([currloc(1:end-6) 'private/m_topo.mat']);
%%load private/m_topo
%load topo

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
if exist('MAP_PROJECTION.routine') && strcmp(MAP_PROJECTION.routine,'mp_azim') && length(lgs)==362 && ...
   strcmp(optn,'contourf')
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

 
% ...and draw if required

if draw_map

  % Set current projection to geographic
  Currentmap=m_coord('set');
  m_coord('geographic');

  hold on;
  switch optn
   case 'contour'
       if ~isempty(varargin) && ~ischar(varargin{1})
           levels=varargin{1};
           varargin(1)=[];
       end
      [values,longs]=m_contour(lgs,lts,topo,levels);
      if ~isempty(varargin)
          for l=1:length(longs) 
              set(longs(l),varargin{:}); 
          end 
      end
   case 'contourf'
       if ~isempty(varargin) && ~ischar(varargin{1})
           levels=varargin{1};
           varargin(1)=[];
       end
      [values,longs]=m_contourf(lgs,lts,topo,levels);
      if ~isempty(varargin)
          for l=1:length(longs) 
              set(longs(l),varargin{:}); 
          end 
      end
   case 'pcolor'
      [longs]=m_pcolor(lgs,lts,topo,varargin{:});
   case 'image'
      [longs]=m_image(lgs,lts,topo,varargin{:});
   case 'shadedrelief'
       [Im,X,Y]=m_image(lgs,lts,topo);
       [longs]=m_shadedrelief(X,Y,Im,'coords','map','gradient',1,varargin{:});
   otherwise
      error(['m_elev:Unrecognized option: ' optn]);       
  end  

  set(longs,'tag','m_elev');  
  
  % Reset map coords
  m_coord(Currentmap.name);
 
else  % Just the data ma'am...

  [longs,lats]=meshgrid(lgs,lts);
  values=topo;

end


if nargout==0
 clear values lats longs
end
