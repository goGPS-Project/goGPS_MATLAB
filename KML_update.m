function KML_update (filename,lam,phi,h,nsat,date)

% SYNTAX:
%   KML_update (filename,lam,phi,h,nsat,date);
%
% INPUT:
%   filename = name of the file with extension
%   lam = longitude [degrees]
%   phi = latitude [degrees]
%   h = orthometric height [m]
%   nsat = number of visible satellites
%   date = date expressed as [year,month,day,hour,minutes,seconds)
%
% DESCRIPTION:
%   Update a KML file (Goole Earth).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini*
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
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

%-------------------------------------------------------------------------------
% PARAMETERS FOR THE KML FILE
%-------------------------------------------------------------------------------

%"clampedToGround" plots the points attached to the ground
%"absolute" uses the height defined in the tag <coordinates>;
%N.B. Google Earth uses orthometric heights
z_pos = 'clampedToGround';
%z_pos = 'absolute';

%URL to load the icon for the points
icon = 'http://maps.google.com/mapfiles/kml/pal2/icon26.png';

%string representing the ARGB color of the points
if (nsat >= 4)
    point_color = 'FFF5005A';
else
    point_color = 'FF0000FF';
end

%point size
scale = 0.4;

%label color
label_color = point_color;

%label size
label_scale = 0.7;

%-------------------------------------------------------------------------------
% FORMAT DATE FOR KML
%-------------------------------------------------------------------------------

year = num2str(date(1,1));
month = num2str(date(1,2));
day = num2str(date(1,3));
hour = num2str(date(1,4));
minute = num2str(date(1,5));
second = num2str(floor(date(1,6)));

%[second] = num2str(floor(str2double(second)));
[null, ncifre] = size(year);
if (ncifre == 1)
    [year] = sprintf('200%s',year);
elseif (ncifre == 2)
    [year] = sprintf('20%s',year);
end

[null, ncifre] = size(month);
if (ncifre == 1)
    [month] = sprintf('0%s',month);
end

[null, ncifre] = size(day);
if (ncifre == 1)
    [day] = sprintf('0%s',day);
end

[null, ncifre] = size(hour);
if (ncifre == 1)
    [hour] = sprintf('0%s',hour);
end

[null, ncifre] = size(minute);
if (ncifre == 1)
    [minute] = sprintf('0%s',minute);
end

[null, ncifre] = size(second);
if (ncifre == 1)
    [second] = sprintf('0%s',second);
end

[datekml] = sprintf('%s-%s-%sT%s:%s:%sZ',year,month,day,hour,minute,second);

%-------------------------------------------------------------------------------
% FILE UPDATE
%-------------------------------------------------------------------------------

fkml=fopen(filename,'r+t');

while (fkml == -1)
    fkml=fopen(filename,'r+t');
end

fseek(fkml,-20,1);

%-------------------------------------------------------------------------------
% INSERTION OF A POINT
%-------------------------------------------------------------------------------

fprintf(fkml, '      <Placemark>\n');
fprintf(fkml, '        <name>%d</name>\n', nsat);
fprintf(fkml, '        <Point>\n');
fprintf(fkml, '          <altitudeMode>%s</altitudeMode>\n',z_pos);
fprintf(fkml, '          <coordinates>%.6f,%.6f,%.6f</coordinates>\n',lam,phi,h);
fprintf(fkml, '        </Point>\n');
fprintf(fkml, '        <Snippet></Snippet>\n');
fprintf(fkml, '        <Style>\n');
fprintf(fkml, '          <IconStyle>\n');
fprintf(fkml, '            <Icon>\n');
fprintf(fkml, '              <href>%s</href>\n',icon);
fprintf(fkml, '            </Icon>\n');
fprintf(fkml, '            <color>%s</color>\n',point_color);
fprintf(fkml, '            <colorMode>normal</colorMode>\n');
fprintf(fkml, '            <scale>%.2f</scale>\n',scale);
fprintf(fkml, '          </IconStyle>\n');
fprintf(fkml, '          <LabelStyle>\n');
fprintf(fkml, '            <color>%s</color>\n',label_color);
fprintf(fkml, '            <scale>%s</scale>\n',label_scale);
fprintf(fkml, '          </LabelStyle>\n');
fprintf(fkml, '        </Style>\n');
fprintf(fkml, '        <TimeStamp>\n');
fprintf(fkml, '          <when>%s</when>\n',datekml);
fprintf(fkml, '        </TimeStamp>\n');
fprintf(fkml, '        <description><![CDATA[ <i>Latitude:</i> %.6f &#176;<br/> <i>Longitude:</i> %.6f &#176;<br/> <i>Elevation:</i> %.1f m<br/> <i>Satellites:</i> %d <br/> <i>Time:</i> %s-%s-%s %s:%s:%s]]></description>\n',phi,lam,h,nsat,year,month,day,hour,minute,second);
fprintf(fkml, '      </Placemark>\n');
fprintf(fkml, '  </Document>\n</kml>');
fclose(fkml);
