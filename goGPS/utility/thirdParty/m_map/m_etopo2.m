function [values,longs,lats]=m_etopo2(varargin)
% M_ETOPO2 Contour elevation onto a map using the  ETOPOX database
%        M_ETOPO2 contours elevations at 1000m intervals for the map.
%        M_ETOPO2(OPTN (,LEVELS) (,ARGS,...) ) lets you change various options.
%        if OPTN=='contour', contour lines are drawn. for OPTN=='contourf',
%        filled contours are drawn. LEVELS are the levels used, and ARGS
%        are optional patch arguments of line types, colors, etc. 
%
%        [CS,H]=M_ETOPO2(...) allows access to the return arguments of the
%        contour/contourf call.
%
%        [ELEV,LONG,LAT]=M_ETOPO2([LONG_MIN LONG_MAX LAT_MIN LAT_MAX])
%        extracts elevation data for the given lat/long limits (without plotting).
%
%        Can be configured for etopo2v2 and etopo1, if you want to use
%        that behemoth  (name remains etopo2 for backward compatability)
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
% May/28/2014 - modified for etopo2v2 (which has actualy ben around for years...)
% May/29/2014 - modified for etopo1  



%%% This will have to be set by YOU the USER! ---------------------

%PATHNAME='/ocean/rich/more/mmapbase/etopo2v2/';   % Be sure to end the path with a "/" or
PATHNAME='/ocean/rich/more/mmapbase/etopo1/';   % Be sure to end the path with a "/" or
                                                 % whatever your separator is.
%  Note - etopo2v2 now comes in 4 flavours - both grid-centered
% and cell-centered, with both in either big-endian or little-endian
% formats. Use the grid-centered version, 
% and right now I have used the big-endian version.
%

%%efid=fopen([PATHNAME 'etopo2_2006apr.raw'],'r','b'); % in big-endian format
%%efid=fopen([PATHNAME 'etopo2.i2'],'r','b'); % in big-endian format


%efid=fopen([PATHNAME 'ETOPO2v2g_i2_MSB.bin'],'r','b'); % in big-endian format
%efid=fopen([PATHNAME 'ETOPO2v2c_i2_MSB.bin'],'r','b'); % in big-endian format
efid=fopen([PATHNAME 'etopo1_ice_g_i2.bin'],'r','l'); % apparently little-endian format

% Now, specify whether this file is grid or cell-referenced
grid=1; % 1 for grid reference, 0 for cell - but I haven't gotten cell reference working yet

% And you have to get the resolution right as well!
resolution=1;  % 2 = 2 minute (etopo2), 1 = 1 minute (etopo1)

%---------------------------------------------------------------------


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

ptsperdeg=60/resolution;
nx=360*(ptsperdeg);  % Grid width
ny=180*(ptsperdeg);  % Grid width



if efid==-1
 warning(['Cannot open ' PATHNAME 'etopo2 !! \n   Have you installed the Etopo2 database correctly?' ...
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

% This finds the nearest cell boundaries

if grid

   if draw_map

     blat=floor(MAP_VAR_LIST.lats(1)*ptsperdeg);
     tlat=ceil(MAP_VAR_LIST.lats(2)*ptsperdeg);
     llong=floor(MAP_VAR_LIST.longs(1)*ptsperdeg);
     rlong=ceil(MAP_VAR_LIST.longs(2)*ptsperdeg);

     lngdec=ceil((rlong-llong+1)/decmax);
     latdec=ceil((tlat-blat+1)/decmax);

   else

     blat=floor(varargin{1}(3)*ptsperdeg);
     tlat=ceil(varargin{1}(4)*ptsperdeg);
     llong=floor(varargin{1}(1)*ptsperdeg);
     rlong=ceil(varargin{1}(2)*ptsperdeg);
     lngdec=1;
     latdec=1;

   end

 
   lgs=        [llong:rlong]/ptsperdeg;   
   lts=fliplr( [blat:tlat]/ptsperdeg    );  

   ptsperline=nx+1;
   
else   % Cell entering

   if draw_map

     blat=floor(MAP_VAR_LIST.lats(1)*ptsperdeg-1/2);
     tlat=ceil(MAP_VAR_LIST.lats(2)*ptsperdeg-1/2);
     llong=floor(MAP_VAR_LIST.longs(1)*ptsperdeg-1/2);
     rlong=ceil(MAP_VAR_LIST.longs(2)*ptsperdeg-1/2);

     lngdec=ceil((rlong-llong+1)/decmax);
     latdec=ceil((tlat-blat+1)/decmax);

   else

     blat=floor(varargin{1}(3)*ptsperdeg-1/2);
     tlat=ceil(varargin{1}(4)*ptsperdeg-1/2);
     llong=floor(varargin{1}(1)*ptsperdeg-1/2);
     rlong=ceil(varargin{1}(2)*ptsperdeg-1/2);
     lngdec=1;
     latdec=1;

   end

   % Cell centers are moved 1/60 = .5/30 awa from the cell boundaries
   %
   % (1,1) point is at -180W, 90N.

   lgs=([llong:rlong]+1/2)/ptsperdeg;  % move right
   lts=fliplr( ([blat:tlat]-1/2)/ptsperdeg );  % move down

   ptsperline=nx;
end   

eaxes=[llong+nx/2 rlong+nx/2 ny/2-blat ny/2-tlat];  % indexes of edges (start with 0)

% Get it inside, or just off the right edge if edge-crossing   
if eaxes(2)>nx,  eaxes([1 2])=eaxes([1 2])-nx; end
if eaxes(1)<0,   eaxes([1 2])=eaxes([1 2])+nx; end



if (eaxes(2)>nx  )   % Read it in in 2 pieces!


  nlat=round((eaxes(3)-eaxes(4)))+1;
  nlgr=round( eaxes(2)-nx )+1;
  nlgl=round(nx-eaxes(1)); %%%nlng-nlgr
  nlng=nlgr+nlgl;

  values=zeros(nlat,nlng);
  for ii=[1:nlat]
   fseek(efid,( (ii-1+eaxes(4))*ptsperline )*2,'bof');
   values(ii,nlng+[-nlgr:-1]+1)=fread(efid,[1 nlgr],'int16');
   fseek(efid,( (ii-1+eaxes(4))*ptsperline+eaxes(1) )*2,'bof');
   values(ii,1:nlgl)=fread(efid,[1 nlgl],'int16');
  end

else  % Read it in one piece

  nlat=round((eaxes(3)-eaxes(4)))+1;
  nlng=round((eaxes(2)-eaxes(1)))+1;
  values=zeros(nlat,nlng);
  for ii=[1:nlat]
   fseek(efid,( (ii-1 +eaxes(4))*ptsperline +eaxes(1) )*2,'bof');
   values(ii,:)=fread(efid,[1 nlng],'int16');
  end

end
 
 

if draw_map

  % Set current projection to geographic
  Currentmap=m_coord('set');
  m_coord('geographic');
 
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
  lts=lts(1:latdec:end);
  lgs=lgs(1:lngdec:end);
   
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
 set(longs,'tag','m_etopo2');
 if n_opt<length(varargin) 
    for l=1:length(longs), set(longs(l),varargin{n_opt:end}); end 
 end
  
 m_coord(Currentmap.name);

else

  [longs,lats]=meshgrid(lgs,lts);

end


if nargout==0
 clear values longs lats
end
