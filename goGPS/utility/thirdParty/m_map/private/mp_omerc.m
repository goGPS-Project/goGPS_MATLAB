function  [X,Y,vals,labI]=mp_omerc(optn,varargin)
% MP_TMERC   Oblique Mercator projection
%           This function should not be used directly; instead it is
%           is accessed by various high-level functions named M_*.


% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% Mathematical formulas for the projections and thie inverses are taken from
%
%      Snyder, John P., Map Projections used by the US Geological Survey, 
%      Geol. Surv. Bull. 1532, 2nd Edition, USGPO, Washington D.C., 1983.
%
% The oblique mercator is a cylindrical projection around an arbitrary
% great-circle route. This is a handy projection for coastlines and so
% forth.
% Oblique Mercator - cylindrical conformal


global MAP_PROJECTION MAP_VAR_LIST

name={'Oblique Mercator'};

pi180=pi/180;

switch optn

  case 'name'

     X=name;

  case {'usage','set'}

     X=char({['     ''' varargin{1} ''''],...
              '     <,''lon<gitude>'',[value1 value2]>',...
              '     <,''lat<itude>'',[value1 value2]>',...
              '     <,''asp<ect>'',value>',...
              '     <,''dir<ection>'',( ''horizontal'' | ''vertical'' )'});

  case 'get'

     X=char([' Projection: ' MAP_PROJECTION.name '  (function: ' MAP_PROJECTION.routine ')'],...
            [' longitudes: ' num2str(MAP_VAR_LIST.ulongs)],...
            [' latitudes: ' num2str(MAP_VAR_LIST.ulats) ],...
            [' Aspect ratio: ' num2str(MAP_VAR_LIST.aspect)],...
            [' Baseline direction ' MAP_VAR_LIST.direc]); 

  case 'initialize'

    MAP_VAR_LIST=[];
    MAP_PROJECTION.name=varargin{1};
    MAP_VAR_LIST.ulongs=[-132 -125];
    MAP_VAR_LIST.ulats=[56 40];
    MAP_VAR_LIST.aspect=0.5;
    MAP_VAR_LIST.direc='ver';
    MAP_VAR_LIST.rectbox='on';  % THis is always the case for this projection; it's just
                                % too difficult to comtemplate the other possibility
    k=2;
    while k<length(varargin)   
      switch varargin{k}(1:3)
         case 'lon'
           MAP_VAR_LIST.ulongs=varargin{k+1};
         case 'asp'
           MAP_VAR_LIST.aspect=varargin{k+1};
         case 'lat'
           MAP_VAR_LIST.ulats=varargin{k+1};
         case 'dir'
           MAP_VAR_LIST.direc=varargin{k+1};
         otherwise
           disp(['Unknown option: ' varargin{k}]);
      end
      k=k+2;
    end
    rulats=MAP_VAR_LIST.ulats*pi180;
    rulongs=MAP_VAR_LIST.ulongs*pi180;

    % Compute the poles of the oblique transformation

    MAP_VAR_LIST.rpolelong=atan2( ...
             cos(rulats(1))*sin(rulats(2))*cos(rulongs(1)) ...
            -sin(rulats(1))*cos(rulats(2))*cos(rulongs(2)),...
             sin(rulats(1))*cos(rulats(2))*sin(rulongs(2)) ...
            -cos(rulats(1))*sin(rulats(2))*sin(rulongs(1)) );
    MAP_VAR_LIST.rpolelat=atan(-cos(MAP_VAR_LIST.rpolelong-rulongs(1))/ ...
                                 tan(rulats(1)));

    % Now get the map X/Y limits in the transformed plane
    switch MAP_VAR_LIST.direc(1:3)
      case 'hor'
        [MAP_VAR_LIST.xlims,bY]=mp_omerc('ll2xy',MAP_VAR_LIST.ulongs,MAP_VAR_LIST.ulats,'clip','off');
        MAP_VAR_LIST.xlims=[min(MAP_VAR_LIST.xlims) max(MAP_VAR_LIST.xlims)];
        MAP_VAR_LIST.ylims=diff(MAP_VAR_LIST.xlims)*[-0.5 0.5]*MAP_VAR_LIST.aspect;
      case 'ver'
        [bX,MAP_VAR_LIST.ylims]=mp_omerc('ll2xy',MAP_VAR_LIST.ulongs,MAP_VAR_LIST.ulats,'clip','off');
        MAP_VAR_LIST.ylims=[min(MAP_VAR_LIST.ylims) max(MAP_VAR_LIST.ylims)];
        MAP_VAR_LIST.xlims=diff(MAP_VAR_LIST.ylims)*[-0.5 0.5]*MAP_VAR_LIST.aspect;
    end
  
    % For further use, it is useful to have the max/min lat/longs in the visible area

    mu_util('lllimits');
    
  case 'll2xy'

    l_0=MAP_VAR_LIST.rpolelong+pi/2;
    long=varargin{1}*pi180-l_0;
    lat=varargin{2}*pi180;
    vals=zeros(size(long));

    A=sin(MAP_VAR_LIST.rpolelat)*sin(lat)- ...
      cos(MAP_VAR_LIST.rpolelat)*cos(lat).*sin(long);

    switch MAP_VAR_LIST.direc(1:3)
      case 'ver'
        Y=atan2( tan(lat)*cos(MAP_VAR_LIST.rpolelat)+ ...
                 sin(MAP_VAR_LIST.rpolelat)*sin(long), cos(long) );
        X=-atanh(A);
      case 'hor'
        X=atan2( tan(lat)*cos(MAP_VAR_LIST.rpolelat)+ ...
                 sin(MAP_VAR_LIST.rpolelat)*sin(long), cos(long) );
        Y=atanh(A);
    end

    if ~strcmp(varargin{4},'off')
        vals= vals | X<=MAP_VAR_LIST.xlims(1)+eps*10 | X>=MAP_VAR_LIST.xlims(2)-eps*10 | ...
                     Y<=MAP_VAR_LIST.ylims(1)+eps*10 | Y>=MAP_VAR_LIST.ylims(2)-eps*10;
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(1),X<MAP_VAR_LIST.xlims(1),Y);
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(2),X>MAP_VAR_LIST.xlims(2),Y);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(1),Y<MAP_VAR_LIST.ylims(1),X);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(2),Y>MAP_VAR_LIST.ylims(2),X);
    end

  case 'xy2ll'
    
     l_0=MAP_VAR_LIST.rpolelong+pi/2;
     
     switch MAP_VAR_LIST.direc(1:3)
       case 'hor'
         Y=asin( sin(MAP_VAR_LIST.rpolelat)*tanh(varargin{2}) ...
                +cos(MAP_VAR_LIST.rpolelat)*sin(varargin{1})./cosh(varargin{2}) )/pi180;
         X=(l_0+atan2( sin(MAP_VAR_LIST.rpolelat)*sin(varargin{1}) ...
                     -cos(MAP_VAR_LIST.rpolelat)*sinh(varargin{2}), cos(varargin{1}) ) )/pi180;
       case 'ver'
         Y=asin( sin(MAP_VAR_LIST.rpolelat)*tanh(-varargin{1}) ...
                +cos(MAP_VAR_LIST.rpolelat)*sin(varargin{2})./cosh(-varargin{1}) )/pi180;
         X=(l_0+atan2( sin(MAP_VAR_LIST.rpolelat)*sin(varargin{2}) ...
                     -cos(MAP_VAR_LIST.rpolelat)*sinh(-varargin{1}), cos(varargin{2}) ) )/pi180;
     end  
        
  case 'xgrid'
   
    [X,Y,vals,labI]=mu_util('xgrid',MAP_VAR_LIST.longs,MAP_VAR_LIST.lats,varargin{1},31,varargin{2:3});
    

  case 'ygrid'
   
    [X,Y,vals,labI]=mu_util('ygrid',MAP_VAR_LIST.lats,MAP_VAR_LIST.longs,varargin{1},31,varargin{2:3});

  case 'box'

    [X,Y]=mu_util('box',2);

end


