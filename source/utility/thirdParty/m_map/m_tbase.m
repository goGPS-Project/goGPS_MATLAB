function [values,longs,lats]=m_tbase(varargin)
% M_TBASE Contour elevation onto a map using the 5-minute TerrainBase database
%        M_TBASE contours elevations at 1000m intervals for the map.
%        M_TBASE(OPTN (,LEVELS) (,ARGS,...) ) lets you change various options.
%        if OPTN=='contour', contour lines are drawn. for OPTN=='contourf',
%        filled contours are drawn. LEVELS are the levels used, and ARGS
%        are optional patch arguments of line types, colors, etc. 
%
%        [CS,H]=M_TBASE(...) allows access to the return arguments of the
%        contour/contourf call.
%
%        [ELEV,LONG,LAT]=M_TBASE([LONG_MIN LONG_MAX LAT_MIN LAT_MAX])
%        extracts elevation data for the given lat/long limits (without plotting).
%
%        See also M_PROJ, M_GRID, M_COAST


% Rich Pawlowicz (rich@ocgy.ubc.ca) 10/Sep/97
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% 17/1/98 - Allowed output of raw data, fixed small bug in selection that left
%           some things off by 1/12 deg lat.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 28/Mar/04 - defaulted to m_elev database (prevents problems with m-demo)
% 4/DEc/11 - isstr to ischar

%%% This will have to be set by YOU the USER!

PATHNAME='/ocean/rich/more/mmapbase/tbase_5/';   % Be sure to end the path with a "/" or
                                     % whatever your separator is.

%%% You probably won't want to change this...
decmax=500;
% What it does is set an upper limit (of DECMAX points) in the
% displayed resolution of the database. Why do I do this? Because the
% problem with a 5-minute database is that it is often way, way too much detail
% for a simple map of, e.g., the Pacific Ocean. And that means it takes a looong
% time to contour the data, not to mention scads of memory. But feel free
% to change that number to something else if you think I have done something
% unreasonable! (PS - I'd appreciate knowing *why* you changed it).



%%% Don't change anything below this... 

efid=fopen([PATHNAME 'tbase.int'],'r');

if efid==-1
 warning(['Cannot open ' PATHNAME 'tbase.int !! \n   Have you installed the TerrainBase database correctly?' ...
        '\n   This (optional) database must be installed separately - see the M_Map user''s guide for instructions' ...
	'\n   ----Using default elevation database instead']);
 if nargout==0
   m_elev(varargin{:});
 elseif nargout==2
   [values,longs]=m_elev(varargin{:});
 elseif nargout==3	
   [values,longs,lats]=m_elev(varargin{:});
 end
 return;
end


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

% Set current projection to geographic
Currentmap=m_coord('set');
m_coord('geographic');


if draw_map

  blat=max(floor(MAP_VAR_LIST.lats(1)*12),-90*12+1);
  tlat=ceil(MAP_VAR_LIST.lats(2)*12);
  llong=floor(MAP_VAR_LIST.longs(1)*12);
  rlong=ceil(MAP_VAR_LIST.longs(2)*12);

  lngdec=ceil((rlong-llong)/decmax);
  latdec=ceil((tlat-blat)/decmax);

else

  blat=max(floor(varargin{1}(3)*12),-90*12+1);
  tlat=ceil(varargin{1}(4)*12);
  llong=floor(varargin{1}(1)*12);
  rlong=ceil(varargin{1}(2)*12);
  lngdec=1;
  latdec=1;

end

lgs=[llong:lngdec:rlong]/12;
lts=fliplr([blat:latdec:tlat]/12);

if rlong>(360*12-1), rlong=rlong-360*12; llong=llong-360*12; end
%%if llong<-360*12, rlong=rlong+360*12; llong=llong+360*12; end;
if llong<0, rlong=rlong+360*12; llong=llong+360*12; end

eaxes=[llong rlong 90*12-blat 90*12-tlat];


if (eaxes(2)>4319 )   % Read it in in 2 pieces!


  nlat=round((eaxes(3)-eaxes(4)))+1;
  nlgr=round( eaxes(2)-4320 )+1;
  nlgl=4320-eaxes(1); %%%nlng-nlgr
  nlng=nlgr+nlgl;

  values=zeros(nlat,nlng);
  for ii=[1:nlat]
   fseek(efid,(ii-1+eaxes(4))*(360*12*2),'bof');
   values(ii,nlng+[-nlgr:-1]+1)=fread(efid,[1 nlgr],'int16');
   fseek(efid,(ii-1+eaxes(4))*(360*12*2)+eaxes(1)*2,'bof');
   values(ii,1:nlgl)=fread(efid,[1 nlgl],'int16');
  end

else  % Read it in one piece

  nlat=round((eaxes(3)-eaxes(4)))+1;
  nlng=round((eaxes(2)-eaxes(1)))+1;
  values=zeros(nlat,nlng);
  for ii=[1:nlat]
   fseek(efid,(ii-1+eaxes(4))*(360*12*2)+eaxes(1)*2,'bof');
   values(ii,:)=fread(efid,[1 nlng],'int16');
  end

end


if draw_map

  if nargin==0
   levels=[-7000:1000:-1000 000:1000:5000];
   optn='contour';
   n_opt=1;
  else
   if ischar(varargin{1})
     optn=varargin{1};
   end
   if nargin==1
     levels=[-7000:1000:-1000 000:1000:5000];
     n_opt=2;
   else
     if ischar(varargin{2})
       levels=[-7000:1000:-1000 000:1000:5000];
       n_opt=2;
    else
       levels=varargin{2};
       n_opt=3;
     end
   end
  end

  topo=values(1:latdec:end,1:lngdec:end);

  if all(levels<0)
   topo=-topo;
   levels=-levels;
  end


  hold on;
  switch optn
   case 'contour'
      [values,longs]=m_contour(lgs,lts,topo,levels);
   case 'contourf'
      [values,longs]=m_contourf(lgs,lts,topo,levels);
  end  
  set(longs,'tag','m_tbase');
  if n_opt<length(varargin)
    for l=1:length(longs), set(longs(l),varargin{n_opt:end}); end
  end

else

  [longs,lats]=meshgrid(lgs,lts);

end

m_coord(Currentmap.name);

if nargout==0
 clear values longs lats
end
