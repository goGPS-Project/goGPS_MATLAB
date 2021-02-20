function r=geodistance(ci,cf,m)

%GEODISTANCE: Calculates the distance in meters between two points on earth surface.
%
% Usage:  r = geodistance( coordinates1 , coordinates2 , method ) ; 
%         
%	  Where coordinates1 = [longitude1,latitude1] defines the
%	  initial position and coordinates2 = [longitude2,latitude2]
%	  defines the final position.
%	  Coordinates values should be specified in decimal degrees.
%	  Method can be an integer between 1 and 23, default is m = 6. 
%         Methods 1 and 2 are based on spherical trigonometry and a 
%         spheroidal model for the earth, respectively.  
%	  Methods 3 to 24 use Vincenty's formulae, based on ellipsoid 
%         parameters. 
%         Here it follows the correspondence between m and the type of 
%         ellipsoid:
%
%         m =  3 -> ANS ,        m =  4 -> GRS80,    m = 5 -> WGS72, 
%         m =  6 -> WGS84,       m =  7 -> NSWC-9Z2, 
%         m =  8 -> Clarke 1866, m =  9 -> Clarke 1880,
%         m = 10 -> Airy 1830,
%         m = 11 -> Bessel 1841 (Ethiopia,Indonesia,Japan,Korea),
%         m = 12 -> Bessel 1841 (Namibia),
%         m = 13 -> Sabah and Sarawak (Everest,Brunei,E.Malaysia),
%         m = 14 -> India 1830, m = 15 -> India 1956, 
%         m = 16 -> W. Malaysia and Singapore 1948, 
%         m = 17 -> W. Malaysia 1969, 
%         m = 18 -> Helmert 1906, m = 19 -> Helmert 1960,
%         m = 20 -> Hayford International 1924, 
%         m = 21 -> Hough 1960, m = 22 -> Krassovsky 1940,
%         m = 23 -> Modified Fischer 1960, 
%         m = 24 -> South American 1969. 
%
%	  Important notes:
%
%	 1)South latitudes are negative.
%	 2)East longitudes are positive.
%	 3)Great circle distance is the shortest distance between two points 
%          on a sphere. This coincides with the circumference of a circle which 
%          passes through both points and the centre of the sphere.
%	 4)Geodesic distance is the shortest distance between two points on a spheroid.
%	 5)Normal section distance is formed by a plane on a spheroid containing a 
%          point at one end of the line and the normal of the point at the other end. 
%          For all practical purposes, the difference between a normal section and a 
%          geodesic distance is insignificant.
%	 6)The method m=2 assumes a spheroidal model for the earth with an average 
%          radius of 6364.963 km. It has been derived for use within Australia. 
%          The formula is estimated to have an accuracy of about 200 metres over 50 km, 
%          but may deteriorate with longer distances. 
%          However, it is not symmetric when the points are exchanged. 
%  
%  Examples: A = [150 -30]; B = [150 -31]; L = [151 -80];
%            [geodistance(A,B,1) geodistance(A,B,2) geodistance(A,B,3)]
%            [geodistance(A,L,1) geodistance(A,L,2) geodistance(A,L,3)]
%            geodistance([0 0],[2 3])
%            geodistance([2 3],[0 0])
%            geodistance([0 0],[2 3],1)
%            geodistance([2 3],[0 0],1)
%            geodistance([0 0],[2 3],2)
%            geodistance([2 3],[0 0],2)
%            for m = 1:24
%            r(m) = geodistance([150 -30],[151 -80],m);
%            end
%            plot([1:m],r), box on, grid on
%
%***************************************************************************************
% Second version: 07/11/2007
% Third  version: 03/08/2010
% 
% Contact: orodrig@ualg.pt
% 
% Any suggestions to improve the performance of this 
% code will be greatly appreciated. 
% 
% Reference: Geodetic Calculations Methods
%            Geoscience Australia
%            (http://www.ga.gov.au/geodesy/calcs/)
%
%***************************************************************************************
% Changed a little bit by Mao,
% So, it works in vector computaiton now.
%
% Contact: argansos@hotmail.com
% Ganquan Mao, 15/09/2011
%
%***************************************************************************************

if size(ci,2)~=2 || size(cf,2)~=2
    error('ci cf must be a n*2 matrix!')
end

if size(ci,1)~=size(cf,1)
    error('ci cf must have same size!')
end

r=[];%for convergence

if nargin==2
    m=6;
end 

lambda1=ci(:,1)*pi/180; 
   phi1=ci(:,2)*pi/180;

lambda2=cf(:,1)*pi/180; 
   phi2=cf(:,2)*pi/180; 

L=lambda2-lambda1;

if m==1%great circle distance, based on spherical trigonometry
    r=180*1.852*60*acos(sin(phi1).*sin(phi2)...
      +cos(phi1).*cos(phi2).*cos(lambda2-lambda1))./pi;
    r=1000*abs(r);

elseif m==2%spheroidal model for the earth
    term1=111.08956*(ci(:,2)-cf(:,2)+0.000001);
    term2=cos(phi1+((phi2-phi1)/2));
    term3=(cf(:,1)-ci(:,1)+0.000001)/(cf(:,2)-ci(:,2)+0.000001);
    r=1000*abs(term1./cos(atan(term2.*term3)));

else%apply Vincenty's formulae (as long as the points are not coincident)
    alla=[0,0,6378160,6378137,6378135,6378137,6378145,6378206.4,...
          6378249.145,6377563.396,6377397.155,6377483.865,6377298.556,...
          6377276.345,6377301.243,6377304.063,6377295.664,6378200,...
          6378270,6378388,6378270,6378245,6378155,6378160];

    allf=[0,0,1/298.25,1/298.257222101,1/298.26,1/298.257223563,...
          1/298.25,1/294.9786982,1/293.465,1/299.3249646,1/299.1528128,...
          1/299.1528128,1/300.8017,1/300.8017,1/300.8017,1/300.8017,...
          1/300.8017,1/298.3,1/297,1/297,1/297,1/298.3,1/298.3,1/298.25];

    a=alla(m);
    f=allf(m);

    b=a*(1-f);

    axa=a^2;
    bxb=b^2;

    U1=atan((1-f)*tan(phi1));
    U2=atan((1-f)*tan(phi2));

    lambda=L;
    lambda_old=sqrt(-1);%there is no way a complex number is going to coincide with a real number!

    ntrials=0;%just in case...

    while any(abs(lambda-lambda_old)>1e-9)
        
        ntrials=ntrials+1;

        lambda_old=lambda;
        sin_sigma=sqrt((cos(U2).*sin(lambda)).^2+(cos(U1).*sin(U2)...
                  -sin(U1).*cos(U2).*cos(lambda)).^2);
        cos_sigma=sin(U1).*sin(U2)+cos(U1).*cos(U2).*cos(lambda);
        sigma=atan2(sin_sigma,cos_sigma);		 	
        sin_alpha=cos(U1).*cos(U2).*sin(lambda)./sin_sigma;
        cos2_alpha=1-sin_alpha.^2;
        cos_2sigmam=cos_sigma-2*sin(U1).*sin(U2)./cos2_alpha;

        C=(f/16)*cos2_alpha.*(4+f*(4-3*cos2_alpha));

        lambda=L+f*(1-C).*sin_alpha.*(sigma+C.*sin_sigma.*(cos_2sigmam...
               +C.*cos_sigma.*(-1+2*(cos_2sigmam).^2)));

       % stop the function if convergence is not achieved
        if ntrials>1000
            disp('Convergence failure...')
            return
        end

    end

    % convergence achieved? get the distance
    uxu=cos2_alpha*(axa-bxb)/bxb;
    A=1+(uxu/16384).*(4096+uxu.*(-768+uxu.*(320-175*uxu)));
    B=(uxu/1024).*(256+uxu.*(-128+uxu.*(74-47*uxu)));
    delta_sigma=B.*sin_sigma.*(cos_2sigmam+(B/4).*(cos_sigma.*(-1+...
                2*cos_2sigmam.^2)-(B/6).*cos_2sigmam.*(-3+...
                4*sin_sigma.^2).*(-3+4*cos_2sigmam.^2)));
    r=b*A.*(sigma-delta_sigma);

end

r(isnan(r))=0;

end