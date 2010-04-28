function KML_write (lamR,phiR,hR,lamM,phiM,hM,nsat,date)

% SYNTAX:
%   KML_write (lamR,phiR,hR,lamM,phiM,hM,nsat,date);
%
% INPUT:
%   lamR = rover longitude [degrees]
%   phiR = rover latitude [degrees]
%   hR   = rover orthometric height [m]
%   lamM = master longitude [degrees]
%   phiM = master latitude [degrees]
%   hM   = master orthometric height [m]
%   nsat = number of visible satellites
%   date = date expressed as [year,month,day,hour,minutes,seconds)
%
% DESCRIPTION:
%   Write a KML file (Goole Earth) displaying rover and master position.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 beta
%
% Copyright (C) 2009-2010 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Graduate School for Creative Cities, Osaka City University, Japan
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

global kml_filename

%-------------------------------------------------------------------------------
% PARAMETERS FOR THE KML FILE
%-------------------------------------------------------------------------------

%"clampedToGround" plots the points attached to the ground
%"absolute" uses the height defined in the tag <coordinates>;
%N.B. Google Earth uses orthometric heights
z_pos = 'clampedToGround';
%z_pos = 'absolute';

%URL to load the icon for the rover position
iconR = 'http://maps.google.com/mapfiles/kml/pal2/icon26.png';

%URL to load the icon for the master station position
iconM = 'http://maps.google.com/mapfiles/kml/shapes/square.png';

%string representing the ARGB color of the points
if (nsat >= 4)
    point_colorR = 'FFF5005A';
else
    point_colorR = 'FF0000FF';
end

%point size
scaleR = 0.4;

%label color
label_colorR = point_colorR;

%label size
label_scaleR = 0.7;

%string representing the ARGB color of the points
point_colorM = 'FF00FFFF';

%point size
scaleM = 0.8;

%label color
label_colorM = point_colorM;

%label size
label_scaleM = 0.7;

%-------------------------------------------------------------------------------
% FORMAT DATE FOR KML
%-------------------------------------------------------------------------------

year = num2str(date(1,1));
month = num2str(date(1,2));
day = num2str(date(1,3));
hour = num2str(date(1,4));
minute = num2str(date(1,5));
second = num2str(floor(date(1,6)));

[null, ncifre] = size(year); %#ok<ASGLU>
if (ncifre == 1)
    [year] = sprintf('200%s',year);
elseif (ncifre == 2)
    [year] = sprintf('20%s',year);
end

[null, ncifre] = size(month); %#ok<ASGLU>
if (ncifre == 1)
    [month] = sprintf('0%s',month);
end

[null, ncifre] = size(day); %#ok<ASGLU>
if (ncifre == 1)
    [day] = sprintf('0%s',day);
end

[null, ncifre] = size(hour); %#ok<ASGLU>
if (ncifre == 1)
    [hour] = sprintf('0%s',hour);
end

[null, ncifre] = size(minute); %#ok<ASGLU>
if (ncifre == 1)
    [minute] = sprintf('0%s',minute);
end

[null, ncifre] = size(second); %#ok<ASGLU>
if (ncifre == 1)
    [second] = sprintf('0%s',second);
end

%-------------------------------------------------------------------------------
% FILE WRITING
%-------------------------------------------------------------------------------

fkml=fopen(kml_filename,'wt');

while (fkml == -1)
    fkml=fopen(kml_filename,'wt');
end

fprintf(fkml, '<?xml version="1.0" standalone="yes"?>\n');
fprintf(fkml, '<kml creator="goGPS" xmlns="http://earth.google.com/kml/2.2">\n');
fprintf(fkml, '  <Document>\n');
fprintf(fkml, '    <name><![CDATA[%s]]></name>\n', kml_filename);
fprintf(fkml, '    <Snippet><![CDATA[created by goGPS]]></Snippet>\n');
if (lamM ~= 0 | phiM ~= 0 | hM ~= 0)
    fprintf(fkml, '      <Placemark>\n');
    fprintf(fkml, '        <name>Master station</name>\n');
    fprintf(fkml, '        <Point>\n');
    fprintf(fkml, '          <altitudeMode>%s</altitudeMode>\n',z_pos);
    fprintf(fkml, '          <coordinates>%.8f,%.8f,%.3f</coordinates>\n',lamM,phiM,hM);
    fprintf(fkml, '        </Point>\n');
    fprintf(fkml, '        <Snippet></Snippet>\n');
    fprintf(fkml, '        <Style>\n');
    fprintf(fkml, '          <IconStyle>\n');
    fprintf(fkml, '            <Icon>\n');
    fprintf(fkml, '              <href>%s</href>\n',iconM);
    fprintf(fkml, '            </Icon>\n');
    fprintf(fkml, '            <color>%s</color>\n',point_colorM);
    fprintf(fkml, '            <colorMode>normal</colorMode>\n');
    fprintf(fkml, '            <scale>%.2f</scale>\n',scaleM);
    fprintf(fkml, '          </IconStyle>\n');
    fprintf(fkml, '          <LabelStyle>\n');
    fprintf(fkml, '            <color>%s</color>\n',label_colorM);
    fprintf(fkml, '            <scale>%s</scale>\n',label_scaleM);
    fprintf(fkml, '          </LabelStyle>\n');
    fprintf(fkml, '        </Style>\n');
    fprintf(fkml, '        <description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation:</i> %.1f m<br/> <i>Time:</i> %s-%s-%s %s:%s:%s]]></description>\n',phiM,lamM,hM,year,month,day,hour,minute,second);
    fprintf(fkml, '      </Placemark>\n');
end
fprintf(fkml, '      <Placemark>\n');
fprintf(fkml, '        <name>%d</name>\n', nsat);
fprintf(fkml, '        <Point>\n');
fprintf(fkml, '          <altitudeMode>%s</altitudeMode>\n',z_pos);
fprintf(fkml, '          <coordinates>%.8f,%.8f,%.3f</coordinates>\n',lamR,phiR,hR);
fprintf(fkml, '        </Point>\n');
fprintf(fkml, '        <Snippet></Snippet>\n');
fprintf(fkml, '        <Style>\n');
fprintf(fkml, '          <IconStyle>\n');
fprintf(fkml, '            <Icon>\n');
fprintf(fkml, '              <href>%s</href>\n',iconR);
fprintf(fkml, '            </Icon>\n');
fprintf(fkml, '            <color>%s</color>\n',point_colorR);
fprintf(fkml, '            <colorMode>normal</colorMode>\n');
fprintf(fkml, '            <scale>%.2f</scale>\n',scaleR);
fprintf(fkml, '          </IconStyle>\n');
fprintf(fkml, '          <LabelStyle>\n');
fprintf(fkml, '            <color>%s</color>\n',label_colorR);
fprintf(fkml, '            <scale>%s</scale>\n',label_scaleR);
fprintf(fkml, '          </LabelStyle>\n');
fprintf(fkml, '        </Style>\n');
fprintf(fkml, '        <description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation:</i> %.1f m<br/> <i>Satellites:</i> %d <br/> <i>Time:</i> %s-%s-%s %s:%s:%s]]></description>\n',phiR,lamR,hR,nsat,year,month,day,hour,minute,second);
fprintf(fkml, '      </Placemark>\n');
fprintf(fkml, '  </Document>\n</kml>');
fclose(fkml);
