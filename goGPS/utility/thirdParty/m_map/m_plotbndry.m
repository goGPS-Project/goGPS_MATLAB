function [bndry_lon,bndry_lat] = m_plotbndry(name,varargin)
% M_PLOTBNDRY plots text files of Lat,Lon for political boundaries.
% Text files (derived from the DCW) are obtained from 
%   http://www.maproom.psu.edu/cgi-bin/ian/points/index.cgi
%
%     M_PLOTBNDRY(NAME) plots the state or country specified in the
%     string NAME, which may include a path name. The routine will
%     then a) search the specified path for a mat-file of that name
%          b) search the specified path for an ascii *2pts.txt file
%             of that name, and, if found convert it to a mat-file.
%          c) failing that, it will open a file dialog box.
%
%     M_PLOTBNDRY(NAME, ...line properties) will use the specified 
%     line properties in drawing the boundary.
%
%     [LON,LAT]=M_PLOTBNDRY(...) returns vectors of the boundary
%     points.
%
% Note: If errors occur when a file is first plotted, check that 
% the entire file was downloaded.  It should end with two consecutive
% END lines.


% Original Author: Michael W. Mann
%
% Changes: R. Pawlowicz 21/12/98 - changed interface to read from
%          given directory, allow output, allow various line 
%          properties to be specified.
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 19/Mar/04 - .mat files not being created because of a bug found by James Connor.


% Set current projection to geographic
Currentmap=m_coord('set');
m_coord('geographic');


targetfile = [name,'2pts.mat']; %try to find binary file
targetpath = dir(targetfile);

if ( length(targetpath)~=0 ) %then load file.
   load(targetfile);
   m_line(bndry_lon,bndry_lat,'tag','m_plotbndry',varargin{:});
   
else %Can't find binary file, load and process text file.
   %  Create .mat file from data read.
   
   targetfile = [name,'2pts.txt'];
   targetpath = dir(targetfile);
   if ( length(targetpath)==0 ) %then
      [filename,pathname]=uigetfile('*2pts.txt',['Select ',targetfile,' file']);
      if ( filename == 0 ), error(['Could not find ',targetfile]), end
      targetpath = [pathname,filename];
   end %if
   [fid,message] = fopen(targetfile,'r');
   if ( fid < 0 ) %then
      fclose(fid);
      error(message);
   end %if
   
   namein = fgetl(fid);
   if ( ~findstr(lower(name),lower(namein)) ) %then
      fclose(fid);
      error(['File contains wrong state! ','Desired: ',name,' Found: ',namein]);
   end %if
   
   done = 0;
   polynum = 0;
   numpts = [];
   while ( ~done )
      %read in polygon index
      line = fgetl(fid);
      field = sscanf(line,'%3c');
      if ( strcmp(upper(field),'END') ) %then found EOF
         done = 1;
      else %start reading in polygon
         polynum = polynum + 1;
         polynumin = sscanf(line,'%u');
         if ( polynum ~= polynumin ) %then
            fclose(fid);
            error(['Didn''t find polynum ',num2str(polynum),...
                  ', found ',num2str(polynumin),' .']);
         end %if
         line = fgetl(fid); %skip polygon centroid
         line = fgetl(fid); %first point
         ptcount = 0;
         while ( ~strcmp(upper(line(1:3)),'END') )
            ptcount = ptcount + 1;
            line = fgetl(fid);
         end %while
         numpts = [numpts,ptcount];
      end %if
   end %while
   
   Npolygon = length(numpts);
   fclose(fid); %Close file.
   
   % Vectors for composite outline
   bndry_lat = zeros(1,sum(numpts)+Npolygon);
   bndry_lon = zeros(1,sum(numpts)+Npolygon);
   i1 = 0; %index to last filled entry in vectors

   fid = fopen(targetfile,'r'); %reopen file
   line = fgetl(fid); %ignore namefield
   
   for i = 1:Npolygon
      line = fgetl(fid); %ignore polygon index
      line = fgetl(fid); %ignore centroid lat,lon
      [poly,count] = fscanf(fid,'%g',[2,inf]); %read data
      line = fgetl(fid); %ignore END
     
      % Save data for reuse as .mat file
      i0 = i1 + 1;  i1 = i0 + numpts(i);
      bndry_lon(i0:i1) = [poly(1,:),NaN]; %add NaN to lift pen
      bndry_lat(i0:i1) = [poly(2,:),NaN]; %add NaN to lift pen
      
   end %for
   
   fclose(fid);
   
   m_line(bndry_lon,bndry_lat,'tag','m_plotbndry',varargin{:});
   nchar = length(targetfile);   % Bug fix thanks to James Connor
   matfile = [targetfile(1:(nchar-4)),'.mat'];
   save(matfile,'bndry_lat','bndry_lon');
   
end %if

m_coord(Currentmap.name);

if nargout==0
 clear bndry_lon bndry_lat
end
