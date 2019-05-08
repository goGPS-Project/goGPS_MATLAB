function  [X,Y,vals,labI]=mp_cyl(optn,varargin)
% MP_CYL   Cylindrical projections
%           This function should not be used directly; instead it is
%           is accessed by various high-level functions named M_*.


% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% 21/5/97 - allow unequal lat limits for mercator and miler projections
%           (to match with user's guide).
%
% Mathematical formulas for the projections and thie inverses are taken from
%
%      Snyder, John P., Map Projections used by the US Geological Survey, 
%      Geol. Surv. Bull. 1532, 2nd Edition, USGPO, Washington D.C., 1983.
%
% These are cylindrical projections wrapped around the equator
% Mercator - conformal
% Miller - "looks" nice.
% Equidistant - basically plotting by lat/long, with distances stretched.
%
%  Jan/18 - added option to put map into units of meters by specifying spheroid.
%  May/18 - rescaled equidistant projection so it does distances properly
%           (multiplied x by cos(lat) instead of dividing y by that value).

global MAP_PROJECTION MAP_VAR_LIST

MAP_ELLIP=mc_ellips;

name={'Mercator','Miller Cylindrical','Equidistant Cylindrical'};

pi180=pi/180;


switch optn

  case 'name'

     X=name;

  case {'usage','set'}

     m_names=fieldnames(MAP_ELLIP);

     X=char({['     ''' varargin{1} ''''],...
              '     <,''lon<gitude>'',( [min max] | center)>',...
              '     <,''lat<itude>'',( maxlat | [min max]>',...
              '     <,''sph<ere>'', one of',...
           reshape(sprintf('         %6s',m_names{:}),15,length(m_names))'} );

  case 'get'

     X=char([' Projection: ' MAP_PROJECTION.name '  (function: ' MAP_PROJECTION.routine ')'],...
            [' longitudes: ' num2str(MAP_VAR_LIST.ulongs) ],...
            [' Latitudes: ' num2str(MAP_VAR_LIST.ulats) ],...  
            [' sphere: ' MAP_VAR_LIST.ellipsoid ]);

  case 'initialize'

    MAP_VAR_LIST=[];
    MAP_PROJECTION.name=varargin{1};
    MAP_VAR_LIST.ulongs=[-180 180];
    MAP_VAR_LIST.ulats=[-85 85];
    MAP_VAR_LIST.clong=0;
    MAP_VAR_LIST.rectbox='off'; %always...(this is because we actually want lat/long grids
                                % to include the edges; this is turned off in rectboxes).
    MAP_VAR_LIST.ellipsoid = 'normal';
    MAP_VAR_LIST.aussiemode=false;
    
    k=2;
    while k<length(varargin)
      switch varargin{k}(1:3)
         case 'lon'
           if length(varargin{k+1})>1
             MAP_VAR_LIST.ulongs=varargin{k+1}(:)';
             MAP_VAR_LIST.clong=mean(MAP_VAR_LIST.ulongs);
           else
             MAP_VAR_LIST.clong=varargin{k+1};
             MAP_VAR_LIST.ulongs=MAP_VAR_LIST.clong+[-180 180];
           end
         case 'lat'
           MAP_VAR_LIST.ulats=varargin{k+1}(:)';
         case 'sph'
           MAP_VAR_LIST.ellipsoid=varargin{k+1};
	 case 'aus' % aussiemode - my joke
	   if strcmp(varargin{k+1},'on')
	     MAP_VAR_LIST.aussiemode=true;
	   end   
         otherwise
           disp(['Unknown option: ' varargin{k}]);
      end
       k=k+2;
    end
   
     if length(MAP_VAR_LIST.ulats)==1
       MAP_VAR_LIST.ulats=[-1 1]*abs(MAP_VAR_LIST.ulats(1));
     end
     
    if ~isfield(MAP_ELLIP,MAP_VAR_LIST.ellipsoid)
       MAP_VAR_LIST.ellipsoid = 'normal';
    end
    MAP_VAR_LIST.ellip=getfield(MAP_ELLIP,MAP_VAR_LIST.ellipsoid);

     mu_util('xylimits');

  case 'll2xy'


    long=varargin{1};
    lat=varargin{2};
    vals=zeros(size(long));

    % Clip out-of-range values
    
    if ~strcmp(varargin{4},'off')
        vals=vals | long<=MAP_VAR_LIST.longs(1)+eps*10 | long>=MAP_VAR_LIST.longs(2)-eps*10 | ...
	              lat<=MAP_VAR_LIST.lats(1)+eps*10 |   lat>=MAP_VAR_LIST.lats(2)-eps*10;
        [long,lat]=mu_util('clip',varargin{4},long,MAP_VAR_LIST.longs(1),long<MAP_VAR_LIST.longs(1),lat);
        [long,lat]=mu_util('clip',varargin{4},long,MAP_VAR_LIST.longs(2),long>MAP_VAR_LIST.longs(2),lat);
        [lat,long]=mu_util('clip',varargin{4},lat,MAP_VAR_LIST.lats(1),lat<MAP_VAR_LIST.lats(1),long);
        [lat,long]=mu_util('clip',varargin{4},lat,MAP_VAR_LIST.lats(2),lat>MAP_VAR_LIST.lats(2),long);
    end


    switch MAP_PROJECTION.name
      case name(1)
        X=MAP_VAR_LIST.ellip(1)*(long-MAP_VAR_LIST.clong)*pi180;
        Y=MAP_VAR_LIST.ellip(1)*atanh(sin(lat*pi180));
      case name(2)
        X=MAP_VAR_LIST.ellip(1)*(long-MAP_VAR_LIST.clong)*pi180;
        Y=MAP_VAR_LIST.ellip(1)*atanh(sin(lat*pi180*0.8))/0.8;
      case name(3)
        X=MAP_VAR_LIST.ellip(1)*(long-MAP_VAR_LIST.clong)*pi180*cos(mean(MAP_VAR_LIST.lats)*pi180);
        Y=MAP_VAR_LIST.ellip(1)*lat*pi180;
    end
    if MAP_VAR_LIST.aussiemode, Y=-Y; X=-X; end;
    
  case 'xy2ll'
  
    if MAP_VAR_LIST.aussiemode, varargin{2}=-varargin{2};varargin{1}=-varargin{1};  end;
    

    switch MAP_PROJECTION.name
      case name(1)
        X=varargin{1}/MAP_VAR_LIST.ellip(1)/pi180+mean(MAP_VAR_LIST.longs);
        Y=90-2/pi180*atan(exp(-varargin{2}/MAP_VAR_LIST.ellip(1)));
      case name(2)
        X=varargin{1}/MAP_VAR_LIST.ellip(1)/pi180+mean(MAP_VAR_LIST.longs);
        Y=(2/pi180*atan(exp(varargin{2}/MAP_VAR_LIST.ellip(1)*0.8))-90)/0.8;
      case name(3)
        X=varargin{1}/MAP_VAR_LIST.ellip(1)/pi180/cos(mean(MAP_VAR_LIST.lats)*pi180)+mean(MAP_VAR_LIST.longs);
        Y=varargin{2}/MAP_VAR_LIST.ellip(1)/pi180;
    end
        
  case 'xgrid'

    [X,Y,vals,labI]=mu_util('xgrid',MAP_VAR_LIST.longs,MAP_VAR_LIST.lats,varargin{1},3,varargin{2:3});

  case 'ygrid'
   
    [X,Y,vals,labI]=mu_util('ygrid',MAP_VAR_LIST.lats,MAP_VAR_LIST.longs,varargin{1},3,varargin{2:3});

  case 'box'

    [X,Y]=mu_util('box',2);

end


