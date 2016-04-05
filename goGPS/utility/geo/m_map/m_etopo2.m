function [values,longs,lats]=m_etopo2(varargin);
% M_ETOPO2 Contour elevation onto a map using the 2-minute ETOPO2 database
%        M_ETOPO2 contours elevations at 1000m intervals for the map.
%        M_ETOPO2(OPTN (,LEVELS) (,ARGS,...) ) lets you change various options.
%        if OPTN=='contour', contour lines are drawn. for OPTN=='contourf',
%        filled contours are drawn. LEVELS are the levels used, and ARGS
%        are optional patch arguments of line types, colors, etc. 
%
%        [CS,H]=M_ETOPO2E(...) allows access to the return arguments of the
%        contour/contourf call.
%
%        [ELEV,LONG,LAT]=M_ETOPO2([LONG_MIN LONG_MAX LAT_MIN LAT_MAX])
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
% 21/Mar/06 - modified for etopo2
% 4/DEc/11 - isstr to ischar

%%% This will have to be set by YOU the USER!

PATHNAME='/ocean/rich/more/mmapbase/etopo2/';   % Be sure to end the path with a "/" or
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

efid=fopen([PATHNAME 'etopo2_2006apr.raw'],'r','b'); % in big-endian format
%efid=fopen([PATHNAME 'etopo2.i2'],'r','b'); % in big-endian format

if efid==-1,
 warning(sprintf(['Cannot open ' PATHNAME 'etopo2 !! \n   Have you installed the TerrainBase database correctly?' ...
        '\n   This (optional) database must be installed separately - see the M_Map user''s guide for instructions' ...
	'\n   ----Using default elevation database instead']));
 if nargout==0,
   m_elev(varargin{:});
  elseif nargout==2,
   [values,longs]=m_elev(varargin{:});
  elseif nargout==3,	
   [values,longs,lats]=m_elev(varargin{:});
  end;	
  return;
end;


global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

draw_map=1;
if nargin==1 & ~ischar(varargin{1}) & length(varargin{1})==4,
  draw_map=0;
end;

if draw_map==1 & isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;



if draw_map,

  blat=floor(MAP_VAR_LIST.lats(1)*30);
  tlat=ceil(MAP_VAR_LIST.lats(2)*30);
  llong=floor(MAP_VAR_LIST.longs(1)*30);
  rlong=ceil(MAP_VAR_LIST.longs(2)*30);

  lngdec=ceil((rlong-llong)/decmax);
  latdec=ceil((tlat-blat)/decmax);

else

  blat=floor(varargin{1}(3)*30);
  tlat=ceil(varargin{1}(4)*30);
  llong=floor(varargin{1}(1)*30);
  rlong=ceil(varargin{1}(2)*30);
  lngdec=1;
  latdec=1;

end;

lgs=[llong:lngdec:rlong]/30;
lts=fliplr([blat:latdec:tlat]/30);

if rlong>(180*30), rlong=rlong-360*30; llong=llong-360*30; end;
if llong<-(180*30), rlong=rlong+360*30; llong=llong+360*30; end;

eaxes=[llong+180*30 rlong+180*30 90*30-blat 90*30-tlat];


if (eaxes(2)>10800 ),   % Read it in in 2 pieces!


  nlat=round((eaxes(3)-eaxes(4)))+1;
  nlgr=round( eaxes(2)-10800 )+1
  nlgl=10800-eaxes(1); %%%nlng-nlgr
  nlng=nlgr+nlgl

  values=zeros(nlat,nlng);
  for ii=[1:nlat],
   fseek(efid,(ii-1+eaxes(4))*((360*30+1)*2),'bof');
   values(ii,nlng+[-nlgr:-1]+1)=fread(efid,[1 nlgr],'int16');
   fseek(efid,(ii-1+eaxes(4))*((360*30+1)*2)+eaxes(1)*2,'bof');
   values(ii,1:nlgl)=fread(efid,[1 nlgl],'int16');
  end;

else  % Read it in one piece

  nlat=round((eaxes(3)-eaxes(4)))+1;
  nlng=round((eaxes(2)-eaxes(1)))+1;
  values=zeros(nlat,nlng);
  for ii=[1:nlat],
   fseek(efid,(ii-1+eaxes(4))*((360*30+1)*2)+eaxes(1)*2,'bof');
   values(ii,:)=fread(efid,[1 nlng],'int16');
  end;

end;


if draw_map,

   % Set current projection to geographic
   Currentmap=m_coord('set');
   m_coord('geographic');
 
   if nargin==0,
   levels=[-7000:1000:-1000 000:1000:5000];
   optn='contour';
   n_opt=1;
  else
   if ischar(varargin{1}),
     optn=varargin{1};
   end;
   if nargin==1,
     levels=[-7000:1000:-1000 000:1000:5000];
     n_opt=2;
   else
     if ischar(varargin{2}),
       levels=[-7000:1000:-1000 000:1000:5000];
       n_opt=2;
    else
       levels=varargin{2};
       n_opt=3;
     end;
   end;
  end;

  topo=values(1:latdec:end,1:lngdec:end);

  if all(levels<0),
   topo=-topo;
   levels=-levels;
  end;
  
 hold on;
 switch optn,
   case 'contour',
      [values,longs]=m_contour(lgs,lts,topo,levels);
   case 'contourf',
      [values,longs]=m_contourf(lgs,lts,topo,levels);
  end;  
  set(longs,'tag','m_etopo2');
  if n_opt<length(varargin), 
    for l=1:length(longs), set(longs(l),varargin{n_opt:end}); end; 
  end;
  
  m_coord(Currentmap.name);

else

  [longs,lats]=meshgrid(lgs,lts);

end;


if nargout==0
 clear values longs lats
end;
