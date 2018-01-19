function [X,Y,vals,labI]=mp_azim(optn,varargin)
% MP_AZIM Azimuthal projections
%           This function should not be used directly; instead it is
%           is accessed by various high-level functions named M_*.


% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
%         13/5/97 - Added satellite perspective
%          1/6/97 - Another stab at removing some /0 errors.
%         10/8/00 - Rotation for projections?
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% Mathematical formulas for the projections and their inverses are taken from
%
%      Snyder, John P., Map Projections used by the US Geological Survey, 
%      Geol. Surv. Bull. 1532, 2nd Edition, USGPO, Washington D.C., 1983.
%
% These are azimuthal projections, best suited for circular areas. The
% stereographic is commonlu used for polar regions.
%   Stereographic  - conformal
%   Orthographic   - neither conformal nor equal-area, but looks like globe
%                    with viewpoint at infinity.
%   Az Equal-area  - equal area, but not conformal (by Lambert)
%   Az Equidistant - distance and direction from center are true 
%   Gnomonic       - all great circles are straight lines.
%   Satellite      - a perspective view from a finite distance

global MAP_PROJECTION MAP_VAR_LIST

name={'Stereographic','Orthographic ','Azimuthal Equal-area','Azimuthal Equidistant','Gnomonic','Satellite'};

pi180=pi/180;


switch optn

  case 'name'
     X=name;

  case {'usage','set'}
     X=char({['     ''' varargin{1} ''''],...
              '     <,''lon<gitude>'',center_long>',...
              '     <,''lat<itude>'', center_lat>',...
              '     <,''rad<ius>'', ( degrees | [longitude latitude] ) | ''alt<itude>'', alt_frac >',...
              '     <,''rec<tbox>'', ( ''on'' | ''off'' | ''circle'' )>',...
	      '     <,''rot<angle>'', degrees CCW>'});

  case 'get'

     X=char([' Projection: ' MAP_PROJECTION.name '  (function: ' MAP_PROJECTION.routine ')'],...
            [' center longitude: ' num2str(MAP_VAR_LIST.ulong) ],...
            [' center latitude: ' num2str(MAP_VAR_LIST.ulat) ],...
            [' radius/altitude : ' num2str(MAP_VAR_LIST.uradius) ],...
            [' Rectangular border: ' MAP_VAR_LIST.rectbox ],...
	    [' Rotation angle: ' num2str(MAP_VAR_LIST.rotang) ]);

  case 'initialize'

    MAP_VAR_LIST=[];
    MAP_PROJECTION.name=varargin{1};
    MAP_VAR_LIST.ulong=0;
    MAP_VAR_LIST.ulat=60;
    MAP_VAR_LIST.rectbox='circle';
    MAP_VAR_LIST.uradius=90;
    MAP_VAR_LIST.rotang=0;
    k=2;
    while k<length(varargin)
      switch varargin{k}(1:3)
         case 'lon'
           MAP_VAR_LIST.ulong=varargin{k+1};
         case 'lat'
           MAP_VAR_LIST.ulat=varargin{k+1};
         case {'rad','alt'}
           MAP_VAR_LIST.uradius=varargin{k+1};
         case 'rec'
           MAP_VAR_LIST.rectbox=varargin{k+1};
	 case 'rot'
	   MAP_VAR_LIST.rotang=varargin{k+1};
         otherwise
           disp(['Unknown option: ' varargin{k}]);
      end
      k=k+2;    
    end
    if strcmp(MAP_VAR_LIST.rectbox,'off'), MAP_VAR_LIST.rectbox='circle'; end

    MAP_VAR_LIST.rlong=MAP_VAR_LIST.ulong*pi180;
    MAP_VAR_LIST.rlat=MAP_VAR_LIST.ulat*pi180;

    % Compute the various limits
    % For the perspective viewpoints, we can't go *quite* to the horizon because this causes
    % problems in the clipping later on - but we can get pretty close (within .995) of it.
    % This is a fudge factor that can probably be changed if we increase the number of points
    % in grid lines...

    if length(MAP_VAR_LIST.uradius)==1
      if strcmp(MAP_PROJECTION.name,name{2}) && abs(MAP_VAR_LIST.uradius-90)<.2
        MAP_VAR_LIST.radius=89.8;
      elseif strcmp(MAP_PROJECTION.name,name{5})
        MAP_VAR_LIST.radius=min(80,MAP_VAR_LIST.uradius);
      elseif strcmp(MAP_PROJECTION.name,name{6})
        MAP_VAR_LIST.radius=acos(1/(1+MAP_VAR_LIST.uradius))/pi180*.98; % uradius is the height fraction here
      else
        MAP_VAR_LIST.radius=MAP_VAR_LIST.uradius;
      end
      rradius=MAP_VAR_LIST.radius*pi180;
    else
      % do some sperical trig
      edge=MAP_VAR_LIST.uradius*pi180 - [MAP_VAR_LIST.rlong 0];
      cosc=sin(MAP_VAR_LIST.rlat)*sin(edge(2))+cos(MAP_VAR_LIST.rlat)*cos(edge(2))*cos(edge(1));
      sinc=sqrt( ( cos(edge(2))*sin(edge(1)))^2 + ... 
                 (cos(MAP_VAR_LIST.rlat)*sin(edge(2))-sin(MAP_VAR_LIST.rlat)*cos(edge(2))*cos(edge(1)))^2);
      rradius=atan2(sinc,cosc);
      MAP_VAR_LIST.radius=rradius/pi180;
    end
    MAP_VAR_LIST.cosradius=cos(rradius);  
 
    switch MAP_PROJECTION.name
      case name(1)
        MAP_VAR_LIST.rhomax=2*tan(rradius/2);
      case name(2)
        MAP_VAR_LIST.rhomax=sin(rradius);
      case name(3)
        MAP_VAR_LIST.rhomax=2*sin(rradius/2);
      case name(4)
        MAP_VAR_LIST.rhomax=rradius;
      case name(5)
        MAP_VAR_LIST.rhomax=tan(rradius);
      case name(6)
        MAP_VAR_LIST.rhomax=sin(rradius)/(1+(1-cos(rradius))/MAP_VAR_LIST.uradius);
    end

    if strcmp(MAP_VAR_LIST.rectbox,'on')
      if length(MAP_VAR_LIST.uradius)==1
        MAP_VAR_LIST.xlims=[-1 1]/sqrt(2)*MAP_VAR_LIST.rhomax;
        MAP_VAR_LIST.ylims=[-1 1]/sqrt(2)*MAP_VAR_LIST.rhomax;        
      else
        [X,Y]=mp_azim('ll2xy',MAP_VAR_LIST.uradius(1),MAP_VAR_LIST.uradius(2),'clip','off');
        MAP_VAR_LIST.xlims=[-abs(X) abs(X)];
        MAP_VAR_LIST.ylims=[-abs(Y) abs(Y)];
      end
    else
      MAP_VAR_LIST.xlims=[-MAP_VAR_LIST.rhomax MAP_VAR_LIST.rhomax];
      MAP_VAR_LIST.ylims=[-MAP_VAR_LIST.rhomax MAP_VAR_LIST.rhomax];
    end

    mu_util('lllimits');
 
    

  case 'll2xy'

    long=varargin{1}*pi180-MAP_VAR_LIST.rlong;
    lat=varargin{2}*pi180;
    vals=zeros(size(long));

    pi180=pi/180;     
    cosc     =sin(MAP_VAR_LIST.rlat)*sin(lat)+cos(MAP_VAR_LIST.rlat)*(cos(lat).*cos(long));
    sinAzsinc=sin(long).*cos(lat);
    cosAzsinc=cos(MAP_VAR_LIST.rlat)*sin(lat)-sin(MAP_VAR_LIST.rlat)*(cos(lat).*cos(long));
    sinc=sqrt(sinAzsinc.^2+cosAzsinc.^2);
  
    switch MAP_PROJECTION.name
      case name(1)
        cosc(cosc==-1)=-1+eps;
        rho=2*sinc./(1+cosc);  % = 2*tan(c/2)
      case name(2)
        rho=sinc;   % = sinc
      case name(3)
        cosc(cosc==-1)=-1+eps;
        rho=sqrt(2)*sinc./sqrt(1+cosc);   % = 2*sin(c/2)
      case name(4)
        rho=atan2(sinc,cosc); % = c
      case name(5)
        rho=sinc./cosc; % = tan(c)
      case name(6)
        rho=sinc./(1+(1-cosc)/MAP_VAR_LIST.uradius); % 
    end

    sinc(sinc==0)=eps;
    Az=(sinAzsinc+sqrt(-1)*cosAzsinc)./sinc;
    Az(abs(Az)==0)=-1;

    % Clip out-of-range values. We test against cos(c) (where c is the angular
    % distance from map center) rather than directly against rhomax, because
    % in the orthographic map rho->0 for points on the other side of the
    % globe whereas c does not!
    
    % Also, we clip on rho even if we later clip on X/Y because in some projections (e.g. the 
    % orthographic) the X/Y locations wrap back. 
    if ~strcmp(varargin{4},'off')
        vals = vals | cosc<=MAP_VAR_LIST.cosradius+eps*10;
        [rho,Az]=mu_util('clip',varargin{4},rho,MAP_VAR_LIST.rhomax,cosc<MAP_VAR_LIST.cosradius,Az);
        Az=Az./abs(Az);
    end

     X=rho.*real(Az*exp(i*pi180*MAP_VAR_LIST.rotang));
     Y=rho.*imag(Az*exp(i*pi180*MAP_VAR_LIST.rotang));

    if strcmp(MAP_VAR_LIST.rectbox,'on')  && ~strcmp(varargin{4},'off')
        vals= vals | X<=MAP_VAR_LIST.xlims(1)+eps*10 | X>=MAP_VAR_LIST.xlims(2)-eps*10 | ...
                     Y<=MAP_VAR_LIST.ylims(1)+eps*10 | Y>=MAP_VAR_LIST.ylims(2)-eps*10;
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(1),X<MAP_VAR_LIST.xlims(1) | isnan(X),Y);
        [X,Y]=mu_util('clip',varargin{4},X,MAP_VAR_LIST.xlims(2),X>MAP_VAR_LIST.xlims(2) | isnan(X),Y);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(1),Y<MAP_VAR_LIST.ylims(1) | isnan(Y),X);
        [Y,X]=mu_util('clip',varargin{4},Y,MAP_VAR_LIST.ylims(2),Y>MAP_VAR_LIST.ylims(2) | isnan(Y),X);
    end

  case 'xy2ll'


    rho=sqrt(varargin{1}.^2+varargin{2}.^2);
    Z=exp(i*(atan2(varargin{2},varargin{1})-MAP_VAR_LIST.rotang*pi180));
    V1=rho.*real(Z);
    V2=rho.*imag(Z);
    
    ir=rho==0;  % To prevent /0 warnings when rho is 0
    rho(ir)=eps;

    switch MAP_PROJECTION.name
      case name(1)
        c=2*atan(rho/2);
      case name(2)
        c=asin(rho);
      case name(3)
        c=2*asin(rho/2);
      case name(4)
        c=rho;
      case name(5)
        c=atan(rho);
      case name(6)
        c=asin((MAP_VAR_LIST.uradius+1)./sqrt(1+(MAP_VAR_LIST.uradius./rho).^2)) - atan(rho/MAP_VAR_LIST.uradius);
    end
    c(ir)=eps; % we offset this slightly so that the correct limit is achieved in the
               % division below:

%    Y=(asin(cos(c)*sin(MAP_VAR_LIST.rlat) + ...
%            cos(MAP_VAR_LIST.rlat)*sin(c).*varargin{2}./rho))/pi180;
%
%    switch MAP_VAR_LIST.ulat,
%      case 90,
%        X=(MAP_VAR_LIST.rlong+atan2(varargin{1},-varargin{2}))/pi180;
%      case -90,
%        X=(MAP_VAR_LIST.rlong+atan2(varargin{1},varargin{2}))/pi180;
%      otherwise
%        X=(MAP_VAR_LIST.rlong+atan2( varargin{1}.*sin(c), ...
%          cos(MAP_VAR_LIST.rlat)*cos(c).*rho - sin(MAP_VAR_LIST.rlat)*varargin{2}.*sin(c) ) )/pi180; 
%     end;

   % Can be problem if the argument is slightly larger than 1 - then the asin
   % returns a complex number.
   arg=cos(c)*sin(MAP_VAR_LIST.rlat) + ...
            cos(MAP_VAR_LIST.rlat)*sin(c).*V2./rho;
    arg=min(max(arg,-1),1);
    
    Y=(asin(arg))/pi180;

    switch MAP_VAR_LIST.ulat
      case 90
        X=(MAP_VAR_LIST.rlong+atan2(V1,-V2))/pi180;
      case -90
        X=(MAP_VAR_LIST.rlong+atan2(V1,V2))/pi180;
      otherwise
        X=(MAP_VAR_LIST.rlong+atan2( V1.*sin(c), ...
          cos(MAP_VAR_LIST.rlat)*cos(c).*rho - sin(MAP_VAR_LIST.rlat)*V2.*sin(c) ) )/pi180; 
    end

  case 'xgrid'
   
    [X,Y,vals,labI]=mu_util('xgrid',MAP_VAR_LIST.longs,MAP_VAR_LIST.lats,varargin{1},31,varargin{2:3});

  case 'ygrid'

    [X,Y,vals,labI]=mu_util('ygrid',MAP_VAR_LIST.lats,MAP_VAR_LIST.longs,varargin{1},91,varargin{2:3});

  case 'box'

     [X,Y]=mu_util('box',31);

end



