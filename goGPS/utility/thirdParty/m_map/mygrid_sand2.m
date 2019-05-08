function  [image_data,vlat,vlon] = mygrid_sand2(region,ssfname)
%  MYGRID_SAND2  Read bathymetry data from Sandwell Database
%    [Z,LAT,LON] = MYGRID_SAND2(REGION) extracts data from
%    the Sandwell and Smith bathymetry, which is now at 1-minute
%    resolution.
%
%
% WARNING: change ssfname and ssversion to the correct one for 
%          your machine
%
%                                               Catherine de Groot-Hedlin
%                                               modified Rich Pawlowicz
%
% latitudes must be between -80.738 and 80.738;
%       input:
%               REGION =[west east south north];
%       output:
%               Z - matrix of sandwell bathymetry/topography
%               LAT - vector of latitudes associated with image_data
%               LON - vector of longitudes
%
%  An odd Z value of say -2001m signifies that this pixel was constrained
%  by a real depth sounding while an even depth of say -2000m is
%  a predicted depth.  Plot rem(Z,2) to see this.
%
%  Changes  Jan/18 - modified to work with post v9 database versions
%                    also order of arguments, and various other
%                    things that bothered me.

if nargin<1
    region=[-140 -121 46 52];
end

if nargin<2
  %ssfname='/ocean/rich/more/mmapbase/s_and_s/topo_8.2.img';
  %ssversion=8.2;
  ssfname='/ocean/rich/more/mmapbase/s_and_s/topo_18.1.img';
  ssversion=18.1;
end


%------no changes below here-------------

% determine the requested region
wlon = region(1);
elon = region(2);
blat = region(3);
tlat = region(4);

% Setup the parameters for reading Sandwell data version 
if ssversion<9
   db_res         = 2/60;          % 2 minute resolution
   db_loc         = [-72.006 72.006 0.0 360-db_res];
   db_size        = [6336 10800];
else
   db_res         = 1/60;          % 1 minute resolution
   db_loc         = [-80.738  80.738 0.0 360-db_res];
   db_size        = [17280 21600];    
end   
nbytes_per_lat = db_size(2)*2;  % 2-byte integers

% Check ranges

if blat<db_loc(1)
    blat=db_loc(1);
end
if tlat>db_loc(2)
    tlat=db_loc(2);
end


% Determine if the database needs to be read twice (overlapping prime meridian)
if ((wlon<0)&&(elon>=0))
      wlon      = [wlon           0];
      elon      = [360-db_res  elon];
end

% Calculate number of "records" down to start (latitude) (0 to db_size(1)-1)
% (mercator projection)
rad=pi/180;

arg1=log(tand(45+db_loc(1)/2));
arg2=log(tand(45+blat/2));
iblat = fix(db_size(1) +1 - (arg2-arg1)/(db_res*rad));

arg2=log(tand(45+tlat/2));
itlat = fix(db_size(1) +1 - (arg2-arg1)/(db_res*rad));

if (iblat < 0 ) || (itlat > db_size(1)-1)
        errordlg([' Requested latitude is out of file coverage ']);
end

% Go ahead and read the database
% Open the data file
fid = fopen( ssfname, 'r','b');
if (fid < 0)
        error(['Could not open database: '  ssfname],'Error');
end

image_data     = [];
for i = 1:length(wlon)


        % Make sure the longitude data goes from 0 to 360
        if wlon(i) < 0
                wlon(i) = 360 + wlon(i);
        end

        if elon(i) < 0
                elon(i) = 360 + elon(i);
        end

        % Calculate the longitude indices into the matrix (0 to db_size(1)-1)
        iwlon(i) = fix((wlon(i)-db_loc(3))/db_res);
        ielon(i) = fix((elon(i)-db_loc(3))/db_res);
        if (iwlon(i) < 0 ) || (ielon(i) > db_size(2)-1)
               error([' Requested longitude is out of file coverage ']);
        end

        % allocate memory for the data
        data = zeros(iblat-itlat+1,ielon(i)-iwlon(i)+1);

        % Skip into the appropriate spot in the file, and read in the data
        disp('Reading in bathymetry data');
        for ilat = itlat:iblat
                offset = ilat*nbytes_per_lat + iwlon(i)*2;
                status = fseek(fid, offset, 'bof');
                data(iblat-ilat+1,:)=fread(fid,[1,ielon(i)-iwlon(i)+1],'integer*2');
        end


        % put the two files together if necessary
        if (i>1)
                image_data = [image_data data];
        else
                image_data = data;
        end
end
% close the file
fclose(fid);


% Determine the coordinates of the image_data
vlat=zeros(1,iblat-itlat+1);
arg2 = log(tand(45+db_loc(1)/2.));
for ilat=itlat+1:iblat+1
        arg1 = rad*db_res*(db_size(1)-ilat+0.5);
        term=exp(arg1+arg2);
        vlat(iblat-ilat+2)=2*atand(term)-90;
end
vlon=db_res*((iwlon+1:ielon+1)-0.5);

% to plot it up

if nargout==0
  [xx,yy]=meshgrid(vlon,vlat);
  pcolor(xx,yy,image_data),shading flat,colormap(jet),colorbar('vert')
  xlabel('longitude'),ylabel('latitude'),title('Smith and Sandwell bathymetry Example')
  clear image_data vlon vlat
end
