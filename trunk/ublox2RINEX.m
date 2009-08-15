function ublox2RINEX(msg, filename)

% SYNTAX:
%   ublox2RINEX(msg, filename);
%
% INPUT:
%
% OUTPUT:
%
% DESCRIPTION:
%   Conversion from ublox binary message to RINEX format.

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

global XM YM ZM

%----------------------------------------------------------------------------------------------

%message decoding
[TOW, WEEK, NSV, RES, CPM, PRM, DOM, MQI, CNO, LLI] = decode_ublox(msg);

%date decoding
DAT = datevec(TOW/(3600*24) + 7*WEEK + datenum([1980,1,6,0,0,0]));

%----------------------------------------------------------------------------------------------

%create RINEX observation file
fid = fopen(filename,'w');	

%write header
fprintf(fid,'     2.10           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE\n');
fprintf(fid,'goGPS               Geomatics Lab.      %11s         PGM / RUN BY / DATE \n', date);
fprintf(fid,'Antenna marker                                              MARKER NAME         \n'); 
fprintf(fid,'Geomatics Lab.      Politecnico Milano                      OBSERVER / AGENCY   \n');
fprintf(fid,'                    ublox AEK-4T                            REC # / TYPE / VERS \n');
fprintf(fid,'                    ANN-MS                                  ANT # / TYPE        \n');
fprintf(fid,'%14.4f%14.4f%14.4f                  APPROX POSITION XYZ \n', XM, YM, ZM);
fprintf(fid,'        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N\n');
fprintf(fid,'     2     0                                                WAVELENGTH FACT L1/2\n');
fprintf(fid,'     4    C1    L1    S1    D1                              # / TYPES OF OBSERV \n');
fprintf(fid,'     1                                                      INTERVAL            \n');
fprintf(fid,'%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS   \n', ...
        DAT(1,1), DAT(1,2), DAT(1,3), DAT(1,4), DAT(1,5), DAT(1,6));
fprintf(fid,'                                                            END OF HEADER       \n');

%-------------------------------------------------------------------------------

%number of records
N = length(TOW);

%write data
for i = 1 : N
    sat = find(PRM(:,i) ~= 0);
    n = length(sat);
    fprintf(fid,' %02d %2d %2d %2d %2d %10.7f  0 %2d', ...
            DAT(i,1)-2000, DAT(i,2), DAT(i,3), DAT(i,4), DAT(i,5), round(DAT(i,6)), n);
    for j = 1 : n
        fprintf(fid,'G%02d',sat(j));
    end
    fprintf(fid,'\n');
    for j = 1 : n
        fprintf(fid,'%14.3f %1d',PRM(sat(j),i),floor(CNO(sat(j),i)/6));
        if (CPM(sat(j),i) > 1e-100)
            fprintf(fid,'%14.3f%1d%1d',CPM(sat(j),i),LLI(sat(j),i),floor(CNO(sat(j),i)/6));
        else
            fprintf(fid,'                ');
        end
        fprintf(fid,'%14.3f %1d',CNO(sat(j),i),floor(CNO(sat(j),i)/6));
        fprintf(fid,'%14.3f %1d',DOM(sat(j),i),floor(CNO(sat(j),i)/6));
        fprintf(fid,'\n');
    end    
end

%-------------------------------------------------------------------------------

%close RINEX file
fclose(fid);
